"""
Receptor conversion to PDBQT, offline.

Why this exists rather than a single call to Meeko.

Meeko prepares a receptor from built-in templates for the twenty
standard amino acids, and that works with no connection. For anything
else -- a magnesium, a heme, bound ADP -- it goes to the internet to
fetch a chemical template from the Chemical Component Dictionary. That
was discovered only when the whole pipeline was run end to end.

For this software that is disqualifying. The specification requires it
to open, teach and run with no connection, and a workshop room with no
wifi is the normal case rather than the exception. Worse, the failure
lands exactly on structures that contain a cofactor -- which are the
structures Stage 3's cofactor-preservation design exists to serve.

So the receptor is built in three parts, all offline:

    protein    Meeko's receptor path, which needs no network
    cofactors  Meeko's ligand path, using templates shipped below
    metals     written directly; a single atom needs no template

and the parts are concatenated. The result is one receptor PDBQT with
the catalytic metals and cofactors still in it, which is the point.
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

# ----------------------------------------------------------------------
# Metals: a single atom needs no template, only a type and a charge.
# AutoDock type, then the usual formal charge in a protein.
# ----------------------------------------------------------------------

METAL_PARAMS: dict[str, tuple[str, float]] = {
    "ZN": ("Zn", 2.0), "MG": ("Mg", 2.0), "MN": ("Mn", 2.0),
    "FE": ("Fe", 3.0), "FE2": ("Fe", 2.0), "CA": ("Ca", 2.0),
    "CU": ("Cu", 2.0), "CU1": ("Cu", 1.0), "NI": ("Ni", 2.0),
    "CO": ("Co", 2.0), "CD": ("Cd", 2.0), "K": ("K", 1.0),
    "NA": ("Na", 1.0),
}

# ----------------------------------------------------------------------
# Cofactor templates, shipped rather than fetched.
#
# This is the list Meeko would otherwise download. It covers what turns
# up in practice; anything else falls back to element typing and the
# user is told so plainly rather than quietly given something wrong.
# ----------------------------------------------------------------------

COFACTOR_SMILES: dict[str, str] = {
    "ADP": "Nc1ncnc2n(cnc12)C1OC(COP(O)(=O)OP(O)(O)=O)C(O)C1O",
    "ATP": "Nc1ncnc2n(cnc12)C1OC(COP(O)(=O)OP(O)(=O)OP(O)(O)=O)C(O)C1O",
    "AMP": "Nc1ncnc2n(cnc12)C1OC(COP(O)(O)=O)C(O)C1O",
    "ANP": "Nc1ncnc2n(cnc12)C1OC(COP(O)(=O)OP(O)(=O)NP(O)(O)=O)C(O)C1O",
    "GTP": "Nc1nc2n(cnc2c(=O)[nH]1)C1OC(COP(O)(=O)OP(O)(=O)OP(O)(O)=O)C(O)C1O",
    "GDP": "Nc1nc2n(cnc2c(=O)[nH]1)C1OC(COP(O)(=O)OP(O)(O)=O)C(O)C1O",
    "NAD": ("NC(=O)c1ccc[n+](c1)C1OC(COP(O)(=O)OP(O)(=O)OCC2OC(n3cnc4c(N)"
            "ncnc34)C(O)C2O)C(O)C1O"),
    "NAP": ("NC(=O)c1ccc[n+](c1)C1OC(COP(O)(=O)OP(O)(=O)OCC2OC(n3cnc4c(N)"
            "ncnc34)C(OP(O)(O)=O)C2O)C(O)C1O"),
    "FMN": ("Cc1cc2nc3c(=O)[nH]c(=O)nc3n(CC(O)C(O)C(O)COP(O)(O)=O)c2cc1C"),
    "SAM": ("C[S+](CCC(N)C(O)=O)CC1OC(n2cnc3c(N)ncnc32)C(O)C1O"),
    "SAH": ("NC(CCSCC1OC(n2cnc3c(N)ncnc32)C(O)C1O)C(O)=O"),
    "PLP": "Cc1ncc(COP(O)(O)=O)c(C=O)c1O",
    "UDP": "OC1C(O)C(COP(O)(=O)OP(O)(O)=O)OC1n1ccc(=O)[nH]c1=O",
    "COA": ("CC(C)(COP(O)(=O)OP(O)(=O)OCC1OC(n2cnc3c(N)ncnc32)C(OP(O)(O)=O)"
            "C1O)C(O)C(=O)NCCC(=O)NCCS"),
    "MGF": "F[Mg-](F)F",
    "BEF": "F[Be-](F)F",
    "ALF": "F[Al-](F)F",
}


@dataclass
class ReceptorResult:
    ok: bool
    path: Path | None = None
    protein_atoms: int = 0
    extras: list[str] = field(default_factory=list)   # what was added
    skipped: list[str] = field(default_factory=list)  # what could not be
    problem: str = ""
    suggestion: str = ""


# ----------------------------------------------------------------------

def _split(pdb_text: str) -> tuple[list[str], dict[str, list[str]]]:
    """Protein ATOM records, and the het groups keyed by label."""
    protein: list[str] = []
    hets: dict[str, list[str]] = {}
    for line in pdb_text.splitlines():
        if line.startswith("ATOM  "):
            protein.append(line)
        elif line.startswith("HETATM"):
            name = line[17:20].strip().upper()
            if name in ("HOH", "WAT", "DOD"):
                continue
            label = f"{name} {line[21]}{line[22:27].strip()}"
            hets.setdefault(label, []).append(line)
    return protein, hets


def _metal_record(line: str, serial: int) -> str | None:
    """A metal is one atom: copy the coordinates, set type and charge."""
    resname = line[17:20].strip().upper()
    params = METAL_PARAMS.get(resname)
    if params is None:
        return None
    ad_type, charge = params
    return (
        f"ATOM  {serial:5d} {line[12:16]}{line[16]}{line[17:20]}"
        f" {line[21]}{line[22:27]}   {line[30:54]}"
        f"  1.00  0.00    {charge:+6.3f} {ad_type:<2}"
    )


def _cofactor_records(label: str, lines: list[str],
                      serial_start: int) -> tuple[list[str], str]:
    """
    Turn a bound cofactor into receptor PDBQT atoms, offline.

    The crystal coordinates carry no bond orders, so a shipped SMILES
    template supplies them. With bond orders right, Meeko's ligand path
    -- which needs no network -- gives correct atom types and charges.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem

    RDLogger.DisableLog("rdApp.*")
    resname = label.split()[0]
    smiles = COFACTOR_SMILES.get(resname)
    if smiles is None:
        return [], f"{label}: no template is bundled for {resname}"

    # Hydrogens off here too, for the same reason as the protein: the
    # template below carries implicit hydrogens, so a fully protonated
    # fragment has more atoms than the template and the match fails.
    # Meeko puts the polar ones back.
    heavy = [ln for ln in lines
             if (ln[76:78].strip() or ln[12:16].strip()[:1]).upper() != "H"]
    block = "\n".join(heavy) + "\nEND\n"
    raw = Chem.MolFromPDBBlock(block, removeHs=True, sanitize=False)
    if raw is None:
        return [], f"{label}: the coordinates could not be read"

    template = Chem.MolFromSmiles(smiles)
    if template is None:
        return [], f"{label}: the bundled template is invalid"

    try:
        mol = AllChem.AssignBondOrdersFromTemplate(template, raw)
        Chem.SanitizeMol(mol)
    except Exception as exc:  # noqa: BLE001
        # Say what actually went wrong. A bare exception name sends the
        # reader looking in the wrong place, which cost real time on
        # this project more than once.
        return [], (f"{label}: does not match the bundled {resname} "
                    f"template - it may be a variant or incomplete "
                    f"({type(exc).__name__}: {str(exc)[:90]})")

    # Hydrogens off for the template match, then back on for Meeko.
    # Both are needed and in this order: the template carries implicit
    # hydrogens so the match fails with explicit ones, and Meeko refuses
    # a molecule that has implicit ones. addCoords places them, since a
    # receptor atom without coordinates is of no use.
    try:
        mol = Chem.AddHs(mol, addCoords=True)
    except Exception as exc:  # noqa: BLE001
        return [], (f"{label}: hydrogens could not be added "
                    f"({type(exc).__name__}: {str(exc)[:90]})")

    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy

        setups = MoleculePreparation().prepare(mol)
        text, ok, err = PDBQTWriterLegacy.write_string(setups[0])
        if not ok:
            return [], f"{label}: {err}"
    except Exception as exc:  # noqa: BLE001
        return [], (f"{label}: Meeko could not prepare it "
                    f"({type(exc).__name__}: {str(exc)[:110]})")

    # Keep only the atom records, renumbered, and mark them as part of
    # the receptor rather than as a separate molecule.
    out: list[str] = []
    serial = serial_start
    for ln in text.splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            out.append("ATOM  " + f"{serial:5d}" + ln[11:])
            serial += 1
    return out, ""


def prepare_receptor_pdbqt(
    prepared_pdb: Path,
    out_path: Path,
    box: dict[str, Any] | None = None,
    timeout: float = 900.0,
) -> ReceptorResult:
    """
    Build a complete receptor PDBQT with its cofactors, offline.
    """
    text = Path(prepared_pdb).read_text(errors="ignore")
    protein, hets = _split(text)
    if not protein:
        return ReceptorResult(
            ok=False,
            problem="There are no protein atoms in this file.",
            suggestion="Check that the protein preparation step finished.",
        )

    work = Path(tempfile.mkdtemp(prefix="recpdbqt_"))
    try:
        # --- protein, through Meeko's receptor path -------------------
        #
        # Hydrogens are stripped first, and this is not optional.
        #
        # Stage 3 adds a full set of hydrogens at pH 7.4, which is right
        # for minimisation and for showing the user what they have. But
        # an AutoDock receptor carries polar hydrogens only, and Meeko
        # adds those itself from its own residue templates. Handing it a
        # fully protonated structure makes those templates fail, with an
        # error that names an arbitrary residue -- a different one on
        # each run -- and says nothing about hydrogens at all.
        #
        # Isolated by comparison rather than by guesswork: with
        # hydrogens it fails whether or not the terminal OXT is present;
        # without them it succeeds. Two earlier theories about this
        # failure, both about the C-terminus, were wrong.
        heavy = [
            ln for ln in protein
            if (ln[76:78].strip() or ln[12:16].strip()[:1]).upper() != "H"
        ]
        protein_pdb = work / "protein.pdb"
        protein_pdb.write_text("\n".join(heavy) + "\nTER\nEND\n")
        stem = work / "protein_out"

        args = [sys.executable, "-m", "meeko.cli.mk_prepare_receptor",
                "--read_pdb", str(protein_pdb), "-o", str(stem),
                "-p", "--allow_bad_res"]
        if box:
            args += ["--box_center", str(box["center_x"]),
                     str(box["center_y"]), str(box["center_z"]),
                     "--box_size", str(box["size_x"]), str(box["size_y"]),
                     str(box["size_z"])]
        try:
            run = subprocess.run(args, capture_output=True, text=True,
                                 timeout=timeout)
        except subprocess.TimeoutExpired:
            return ReceptorResult(
                ok=False,
                problem="Converting the protein took too long.",
                suggestion="The structure may be unusually large.",
            )

        produced = list(work.glob("protein_out*.pdbqt"))
        if not produced:
            msg = (run.stderr or run.stdout or "").strip().splitlines()
            real = [m for m in msg if "Error" in m or "error" in m]
            return ReceptorResult(
                ok=False,
                problem="The protein could not be converted to PDBQT.",
                suggestion=(real[-1][:200] if real else
                            (msg[-1][:200] if msg else "no message")),
            )

        lines = [ln for ln in produced[0].read_text().splitlines()
                 if ln.startswith(("ATOM", "HETATM"))]
        n_protein = len(lines)

        # --- metals and cofactors, appended --------------------------
        extras: list[str] = []
        skipped: list[str] = []
        serial = n_protein + 1

        for label, het_lines in sorted(hets.items()):
            resname = label.split()[0]
            if resname in METAL_PARAMS and len(het_lines) == 1:
                rec = _metal_record(het_lines[0], serial)
                if rec:
                    lines.append(rec)
                    serial += 1
                    extras.append(f"{label} ({METAL_PARAMS[resname][0]})")
                else:
                    skipped.append(f"{label}: unrecognised metal")
                continue

            records, why = _cofactor_records(label, het_lines, serial)
            if records:
                lines.extend(records)
                serial += len(records)
                extras.append(f"{label} ({len(records)} atoms)")
            else:
                skipped.append(why)

        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text("\n".join(lines) + "\nTER\nEND\n")

        return ReceptorResult(
            ok=True, path=out_path, protein_atoms=n_protein,
            extras=extras, skipped=skipped,
        )
    finally:
        import shutil as _sh
        _sh.rmtree(work, ignore_errors=True)

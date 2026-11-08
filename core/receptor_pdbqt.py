"""
Direct PDB to PDBQT conversion for a rigid receptor.

Why this exists instead of calling Meeko.

Meeko builds a receptor by perceiving the whole chain as a polymer and
matching residue templates. That is more than a rigid receptor needs,
and it failed on both test structures in different ways: on 4JSV it
could not template residues once hydrogens were present, and on 1M17 it
raised "Expected 4 paddings for (A:997, A:998)". Two structures, two
unrelated failures, each needing its own investigation. A teaching tool
cannot have a step that fails differently on every protein.

What makes a direct converter safe here is a fact about Vina.

AutoDock Vina's scoring function has five terms -- two gaussians, a
repulsion, a hydrophobic term and a hydrogen-bond term. All of them
depend on atom TYPES and distances. There is no electrostatic term.
Partial charges appear in a PDBQT file because AutoDock 4 used them;
Vina reads them and ignores them. (If this software ever adds AutoDock
4, charges start to matter and this file must be revisited -- which is
why that is written down here rather than assumed.)

So a Vina receptor needs correct atom types and nothing else, and for
the twenty standard amino acids an atom type follows from the residue
name and the atom name. That is a lookup, not a perception problem:
nothing to infer, nothing to fetch, nothing to fail per structure.

Types used, from the AutoDock set:

    C   aliphatic carbon          A   aromatic carbon
    N   nitrogen, donor or none   NA  nitrogen, acceptor
    OA  oxygen, acceptor          SA  sulfur, acceptor
    HD  polar hydrogen (donor)    plus element symbols for metals
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

# ----------------------------------------------------------------------
# Aromatic carbons, by residue and atom name. Everything else carbon in
# a standard residue is aliphatic.
# ----------------------------------------------------------------------

AROMATIC_CARBONS: dict[str, set[str]] = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "HIS": {"CG", "CD2", "CE1"},
    "HID": {"CG", "CD2", "CE1"},
    "HIE": {"CG", "CD2", "CE1"},
    "HIP": {"CG", "CD2", "CE1"},
}

# Nitrogens that can accept a hydrogen bond -- but only when they are
# not carrying a hydrogen themselves.
#
# A neutral histidine has one NH and one bare ring nitrogen. Which is
# which is the tautomer, and it cannot be read off the heavy atoms.
#
# An earlier version typed both as acceptors, reasoning that this was
# the cautious choice. It is not: it invents a second acceptor that
# does not exist, and an invented acceptor can manufacture a hydrogen
# bond and with it a docking score. Comparing against Meeko atom by
# atom showed this as the only disagreement in 8,608 atoms -- 42 of
# them, every one a histidine NE2.
#
# The real answer needs no guessing. The protein-preparation step adds
# hydrogens at a stated pH, so the protonation is already decided and
# recorded in the file: a ring nitrogen with a hydrogen is a donor, one
# without is an acceptor. Determined, not assumed.
#
# Only these nitrogens are candidates. Backbone nitrogen and the amide
# nitrogens of asparagine and glutamine are poor acceptors because the
# lone pair is delocalised into the carbonyl, so they are never typed
# NA however they are protonated.
ACCEPTOR_NITROGENS: dict[str, set[str]] = {
    "HIS": {"ND1", "NE2"},
    "HID": {"ND1", "NE2"},
    "HIE": {"ND1", "NE2"},
    "HIP": set(),            # doubly protonated: both carry an H
}


def _nitrogen_has_hydrogen(atom_name: str, residue_atoms: set[str]) -> bool:
    """
    Is there a hydrogen on this nitrogen, by PDB naming?

    A hydrogen is named after the heavy atom it sits on: HD1 on ND1,
    HE2 on NE2, H on N. That is the PDB convention, so this is a
    lookup rather than a geometric guess.
    """
    suffix = atom_name[1:]                     # "ND1" -> "D1"
    candidates = {f"H{suffix}", f"H{suffix}1", f"H{suffix}2"}
    if atom_name == "N":
        candidates |= {"H", "H1", "H2", "H3", "HN"}
    return bool(candidates & residue_atoms)

# Sulfurs. Both cysteine and methionine sulfur can accept weakly.
ACCEPTOR_SULFURS = {"SG", "SD"}

METALS = {
    "ZN": "Zn", "MG": "Mg", "MN": "Mn", "FE": "Fe", "CA": "Ca",
    "CU": "Cu", "NI": "Ni", "CO": "Co", "CD": "Cd", "K": "K",
    "NA": "Na", "MO": "Mo", "W": "W", "V": "V",
}

WATER = {"HOH", "WAT", "DOD", "H2O"}

STANDARD = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "HID", "HIE", "HIP", "CYX", "MSE", "SEC", "PYL",
}


@dataclass
class Converted:
    ok: bool
    path: Path | None = None
    n_atoms: int = 0
    n_protein: int = 0
    n_hetero: int = 0
    kept_hetero: list[str] = field(default_factory=list)
    type_counts: dict[str, int] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    problem: str = ""
    suggestion: str = ""


# ----------------------------------------------------------------------

def _element(line: str) -> str:
    el = line[76:78].strip().upper()
    if el:
        return el
    # Fall back to the atom name, which is where older files keep it.
    name = line[12:16]
    return (name[0] if name[0] != " " else name[1]).upper()


def _is_polar_hydrogen(line: str, by_residue: dict[str, list[str]]) -> bool:
    """
    A hydrogen counts as a donor when it hangs off nitrogen or oxygen.

    Decided from the atom name rather than from geometry: in every
    standard residue a hydrogen is named after the heavy atom it is
    attached to -- HG on OG, HE1 on NE1, H on N. That naming is part of
    the PDB specification, so it is a lookup rather than a guess.
    """
    name = line[12:16].strip().upper()
    if not name.startswith("H") and not (len(name) > 1 and name[1] == "H"):
        return False
    stem = name.lstrip("0123456789").lstrip("H")
    if not stem:
        return True          # plain backbone "H", which sits on N
    for heavy in by_residue.get("names", []):
        if heavy.startswith(("N", "O")) and heavy[1:] == stem:
            return True
        if heavy in ("N", "O") and stem in ("", heavy):
            return True
    return stem[0] in ("N", "O")


def atom_type(line: str, resname: str, residue_atoms: set[str],
              protonated: bool = True) -> str | None:
    """
    The AutoDock type for one atom, or None if it should be dropped.
    """
    name = line[12:16].strip().upper()
    el = _element(line)

    if el == "H":
        # Only polar hydrogens belong in a Vina receptor.
        stem = name.lstrip("0123456789").lstrip("H")
        if not stem:
            return "HD"                      # backbone amide H, on N
        if stem[0] in ("N", "O"):
            return "HD"
        if stem in residue_atoms and stem[0] in ("N", "O"):
            return "HD"
        return None                          # non-polar, dropped

    if el == "C":
        if name in AROMATIC_CARBONS.get(resname, ()):
            return "A"
        return "C"

    if el == "N":
        # Acceptor only if it is a candidate ring nitrogen AND carries
        # no hydrogen of its own. See ACCEPTOR_NITROGENS above.
        if name in ACCEPTOR_NITROGENS.get(resname, ()):
            if protonated:
                # The file says which nitrogen carries the hydrogen,
                # so the tautomer is known rather than assumed.
                if not _nitrogen_has_hydrogen(name, residue_atoms):
                    return "NA"
                return "N"
            # No hydrogens anywhere in the file, so the tautomer cannot
            # be read. Fall back to the common neutral form, N-epsilon-H:
            # the hydrogen sits on NE2, leaving ND1 free to accept. This
            # is what other tools assume, and getting it backwards -- as
            # a first attempt here did -- swaps which nitrogen can form
            # a hydrogen bond. The caller warns the user, because it is
            # a default rather than a reading.
            return "NA" if name == "ND1" else "N"
        return "N"

    if el == "O":
        return "OA"                          # protein oxygens all accept

    if el == "S":
        return "SA" if name in ACCEPTOR_SULFURS else "S"

    if el in ("SE",):
        return "SA"
    if el in ("F", "CL", "BR", "I"):
        return {"F": "F", "CL": "Cl", "BR": "Br", "I": "I"}[el]
    if el == "P":
        return "P"
    if el in METALS:
        return METALS[el]

    return None


def convert(
    pdb_path: Path,
    out_path: Path,
    keep_hetero: list[str] | None = None,
    drop_waters: bool = True,
) -> Converted:
    """
    Write a Vina-ready receptor PDBQT.

    `keep_hetero` holds labels such as "MG A2602" -- normally whatever
    the protein-preparation step decided to keep. Anything not listed
    is left out, and what was left out is reported.
    """
    text = Path(pdb_path).read_text(errors="ignore")

    # None means "no filter, keep every hetero group".
    # An empty list means "keep none of them" -- which is a real and
    # important request: when re-docking, the ligand must come out of
    # the receptor, and after removing it the keep list is empty.
    #
    # Treating [] as "no filter" left the ligand in the receptor while
    # it was being docked into that same pocket. Vina ran, returned
    # nine ranked poses and a score, and nothing anywhere said the
    # pocket was already full.
    filter_hetero = keep_hetero is not None
    keep = set(keep_hetero or [])

    # Atom names present in each residue, so a hydrogen can be matched
    # to the heavy atom it belongs to.
    residue_atoms: dict[tuple, set[str]] = {}
    for line in text.splitlines():
        if line[:6] in ("ATOM  ", "HETATM"):
            key = (line[21], line[22:27].strip(), line[17:20].strip())
            residue_atoms.setdefault(key, set()).add(line[12:16].strip())

    # Does this file carry hydrogens at all? Needed before any atom is
    # typed, because histidine typing depends on it.
    protonated = any(
        _element(ln) == "H"
        for ln in text.splitlines() if ln[:6] in ("ATOM  ", "HETATM")
    )

    out: list[str] = []
    counts: dict[str, int] = {}
    warnings: list[str] = []
    n_protein = n_hetero = 0
    kept_labels: list[str] = []
    unknown: set[str] = set()
    serial = 0

    for line in text.splitlines():
        rec = line[:6]
        if rec not in ("ATOM  ", "HETATM"):
            continue
        resname = line[17:20].strip().upper()
        chain = line[21]
        resseq = line[22:27].strip()
        label = f"{resname} {chain}{resseq}"

        if resname in WATER:
            if drop_waters:
                continue

        is_std = resname in STANDARD
        if not is_std:
            if filter_hetero and label not in keep:
                continue
            if resname not in METALS and resname not in WATER:
                unknown.add(resname)

        key = (chain, resseq, resname)
        ad = atom_type(line, resname, residue_atoms.get(key, set()),
                       protonated=protonated)
        if ad is None:
            continue

        serial += 1
        counts[ad] = counts.get(ad, 0) + 1
        if is_std:
            n_protein += 1
        else:
            n_hetero += 1
            if label not in kept_labels:
                kept_labels.append(label)

        # Charges are written as zero. Vina has no electrostatic term
        # and ignores them; see this module's docstring.
        out.append(
            f"ATOM  {serial:5d} {line[12:16]}{line[16]}{line[17:20]} "
            f"{chain}{line[22:27]}   {line[30:54]}  1.00  0.00    "
            f"{0.0:+6.3f} {ad:<2}"
        )

    if not out:
        return Converted(
            ok=False,
            problem="No usable atoms were found in this structure.",
            suggestion="Check that the protein preparation step finished.",
        )

    if not protonated:
        warnings.append(
            "This structure has no hydrogens, so the histidine tautomers "
            "could not be read and a default was assumed. If a histidine "
            "lines your binding site, prepare the protein first - that "
            "step decides the protonation properly and this conversion "
            "will then follow it."
        )

    if unknown:
        warnings.append(
            "Typed by element only, because no residue template applies: "
            + ", ".join(sorted(unknown))
            + ". For a cofactor this is usually right; check it if the "
            "group sits in your binding site."
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(out) + "\nTER\nEND\n")

    return Converted(
        ok=True, path=out_path, n_atoms=serial, n_protein=n_protein,
        n_hetero=n_hetero, kept_hetero=kept_labels,
        type_counts=dict(sorted(counts.items())), warnings=warnings,
    )

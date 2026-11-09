"""
Stage 6 -- Dock.

This is the stage everything else exists to protect, and the one whose
output is most often believed without reason. A docking score is a
number with units of energy that is not a measurement of energy, is not
comparable between proteins, and correlates only loosely with activity
even within one target. It looks exactly like a result.

Two controls run before any user result is displayed, and neither can
be turned off:

  RE-DOCKING
      Take the ligand that was crystallised in this pocket, throw away
      its pose, and dock it back. If the software cannot reproduce a
      pose that was measured experimentally, it will not get an unknown
      one right either. Root-mean-square deviation under 2 A is the
      conventional pass mark.

  DECOY BASELINE
      Dock a set of molecules matched to the user's compound on size,
      lipophilicity and flexibility, but otherwise arbitrary. If the
      user's compound does not score better than molecules chosen for
      having nothing to do with the target, the score is telling you
      about its size, not its fit. Most reported "hits" fail this and
      almost nobody runs it.

NOTE ON VERIFICATION
    Everything here except the call to Vina itself was built and tested
    where this was written. AutoDock Vina could not be installed on
    that machine -- no SWIG, no Boost headers, no package route -- so
    the subprocess call and the parsing of its output have not been run
    against the real program. check_docking.py tests exactly that on a
    machine where Vina exists. Run it before trusting this stage.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from core.contract import (
    CheckResult,
    ControlOutcome,
    ProvenanceEntry,
    Readiness,
    Stage,
    StageResult,
    Teaching,
    Tier,
    Verdict,
)
from core.install import how_to_install

RMSD_PASS = 2.0          # angstroms; the conventional threshold
DEFAULT_EXHAUSTIVENESS = 8
DEFAULT_SEED = 0xD0CC          # fixed, and recorded, so a run repeats
DEFAULT_POSES = 9


@dataclass
class Pose:
    rank: int
    score: float                      # kcal/mol, as Vina reports it
    rmsd_lb: float = 0.0
    rmsd_ub: float = 0.0
    block: str = ""                   # the PDBQT text of this pose


@dataclass
class DockResult:
    ok: bool
    poses: list[Pose] = field(default_factory=list)
    log: str = ""
    error: str = ""

    @property
    def best(self) -> Pose | None:
        return self.poses[0] if self.poses else None


# ----------------------------------------------------------------------
# Finding and running Vina
# ----------------------------------------------------------------------

def find_vina() -> str | None:
    """Use the private native binary in a bundle, or a source install's PATH."""
    from core.runtime import resource_root
    names = (("vina.exe", "vina", "AutoDock-Vina", "vina_1.2.5") if sys.platform == "win32"
             else ("vina", "AutoDock-Vina", "vina_1.2.5"))

    def bundled():
        here = resource_root()
        for folder in (here / "bin", here, here / "vendor", here / "Scripts", here / "Library" / "bin"):
            for name in names:
                candidate = folder / name
                if candidate.is_file() and (sys.platform == "win32" or os.access(candidate, os.X_OK)):
                    return str(candidate)
        return None

    if getattr(sys, "frozen", False):
        path = bundled()
        if path:
            return path
    for name in names:
        path = shutil.which(name)
        if path:
            return path
    if not getattr(sys, "frozen", False):
        return bundled()
    return None


def parse_vina_pdbqt(text: str) -> list[Pose]:
    """
    Read the poses out of Vina's output file.

    Each model carries a REMARK line of the form

        REMARK VINA RESULT:    -9.5      0.000      0.000

    which is affinity, then lower- and upper-bound RMSD from the best
    pose. Written defensively: a changed line format should produce no
    poses rather than wrong numbers.
    """
    poses: list[Pose] = []
    current: list[str] = []
    score = rmsd_lb = rmsd_ub = None
    rank = 0

    for line in text.splitlines():
        if line.startswith("MODEL"):
            current, score, rmsd_lb, rmsd_ub = [], None, None, None
        elif line.startswith("REMARK VINA RESULT:"):
            nums = re.findall(r"-?\d+\.?\d*", line[len("REMARK VINA RESULT:"):])
            if len(nums) >= 3:
                score, rmsd_lb, rmsd_ub = (float(nums[0]), float(nums[1]),
                                           float(nums[2]))
        elif line.startswith("ENDMDL"):
            if score is not None:
                rank += 1
                poses.append(Pose(rank=rank, score=score,
                                  rmsd_lb=rmsd_lb or 0.0,
                                  rmsd_ub=rmsd_ub or 0.0,
                                  block="\n".join(current)))
        else:
            current.append(line)

    # A single-model file may have no MODEL/ENDMDL wrapper at all.
    if not poses and "REMARK VINA RESULT:" in text:
        nums = re.findall(r"REMARK VINA RESULT:\s*(-?\d+\.?\d*)", text)
        for i, n in enumerate(nums, start=1):
            poses.append(Pose(rank=i, score=float(n), block=text))
    return poses


def run_vina(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    box: dict,
    out_pdbqt: Path,
    exhaustiveness: int = DEFAULT_EXHAUSTIVENESS,
    n_poses: int = DEFAULT_POSES,
    seed: int = DEFAULT_SEED,
    timeout: float = 900.0,
) -> DockResult:
    """One docking run. Never raises; a failure comes back as text."""
    exe = find_vina()
    if exe is None:
        return DockResult(ok=False, error="AutoDock Vina is not installed.")

    cmd = [
        exe,
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--out", str(out_pdbqt),
        "--center_x", str(box["center_x"]),
        "--center_y", str(box["center_y"]),
        "--center_z", str(box["center_z"]),
        "--size_x", str(box["size_x"]),
        "--size_y", str(box["size_y"]),
        "--size_z", str(box["size_z"]),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(n_poses),
        "--seed", str(seed),
    ]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=timeout)
    except subprocess.TimeoutExpired:
        return DockResult(
            ok=False,
            error=(f"Docking did not finish within {timeout / 60:.0f} "
                   "minutes. The box may be far too large."),
        )
    except Exception as exc:  # noqa: BLE001
        return DockResult(ok=False, error=f"Could not run Vina: {exc}")

    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()
        return DockResult(
            ok=False,
            log=proc.stdout,
            error=("Vina stopped with an error: "
                   + (tail[-1][:200] if tail else "no message")),
        )

    if not out_pdbqt.exists():
        return DockResult(ok=False, log=proc.stdout,
                          error="Vina finished but wrote no poses.")

    poses = parse_vina_pdbqt(out_pdbqt.read_text(errors="ignore"))
    if not poses:
        return DockResult(ok=False, log=proc.stdout,
                          error="Vina's output could not be read.")
    return DockResult(ok=True, poses=poses, log=proc.stdout)


# ----------------------------------------------------------------------
# Symmetry-corrected RMSD
# ----------------------------------------------------------------------

def best_rmsd(reference_sdf: Path, probe_sdf: Path) -> float | None:
    """
    RMSD between two poses of the same molecule, allowing for symmetry.

    This has to be symmetry-aware or it is wrong. A benzene ring turned
    through 180 degrees is the same pose; atom-by-atom RMSD calls it a
    large deviation because atom 3 is now where atom 5 was. Scoring a
    correct re-docking as a failure would make the validation control
    useless in exactly the aromatic, symmetric molecules that most
    drugs are. RDKit's GetBestRMS tries every symmetry-equivalent
    mapping and returns the lowest.
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    ref = Chem.MolFromMolFile(str(reference_sdf), removeHs=True)
    probe = Chem.MolFromMolFile(str(probe_sdf), removeHs=True)
    if ref is None or probe is None:
        return None
    try:
        # CalcRMS, not GetBestRMS. The difference is the whole control.
        #
        # GetBestRMS superimposes the two molecules before measuring, so
        # it reports how similar their SHAPES are. CalcRMS measures them
        # where they actually sit. Re-docking asks whether the ligand
        # landed in the right PLACE, so only the second one answers it.
        #
        # Measured on erlotinib displaced from its crystal pose:
        #
        #     displacement   CalcRMS   GetBestRMS
        #           5 A        5.000        0.000
        #          25 A       25.000        0.000
        #
        # The first version of this function used GetBestRMS and
        # reported a comfortable 1.796 A pass for a run whose pocket was
        # still occupied by the ligand being docked. A pose on the far
        # side of the protein would have passed too.
        return round(rdMolAlign.CalcRMS(probe, ref), 3)
    except Exception:  # noqa: BLE001
        # Different molecules, or no common substructure.
        return None


# ----------------------------------------------------------------------
# Property-matched decoys
# ----------------------------------------------------------------------

DECOY_LIBRARY: list[tuple[str, str]] = [
    # Ordinary, commercially common molecules with no relationship to
    # any particular target, spanning drug-like size and lipophilicity.
    #
    # The span is the point. An earlier version of this list held only
    # small molecules, 9 to 20 heavy atoms. Matched against erlotinib
    # (29 heavy atoms, 10 rotatable bonds) the closest decoy had 18
    # atoms and 4 rotatable bonds, so erlotinib won on contact area
    # alone -- which is exactly the illusion this control exists to
    # expose. A baseline that cannot be matched to the query passes
    # everything, and a control that passes everything is worse than
    # no control, because it manufactures confidence.
    ("c1ccc2[nH]ccc2c1", "indole"),
    ("CCOC(=O)c1ccccc1", "ethyl benzoate"),
    ("CCN(CC)C(=O)c1ccccc1", "a diethylbenzamide"),
    ("O=C(Nc1ccccc1)c1ccccc1", "benzanilide"),
    ("c1ccc(cc1)S(=O)(=O)N1CCCCC1", "a phenylsulfonyl piperidine"),
    ("O=C(OCc1ccccc1)c1ccccc1", "benzyl benzoate"),
    ("CCCCOc1ccc(cc1)C(=O)OCC", "a butoxy ethylbenzoate"),
    ("CC(=O)Nc1ccc(cc1)OCc1ccccc1", "a benzyloxy acetanilide"),
    ("COc1ccc(cc1)C(=O)NCc1ccccc1", "a methoxy benzylbenzamide"),
    ("O=C(NCCc1ccccc1)OCc1ccccc1", "a benzyl carbamate"),
    ("CC(C)(C)OC(=O)N1CCC(CC1)Oc1ccc(cc1)C#N", "a Boc piperidinyl ether"),
    ("COc1ccc(cc1OC)CCNC(=O)CCc1ccccc1", "a dimethoxyphenethyl amide"),
    ("CCOc1ccc(cc1)C(=O)N1CCN(CC1)c1ccccc1", "an ethoxybenzoyl piperazine"),
    ("COc1cc(cc(c1OC)OC)C(=O)NCCc1ccccc1", "a trimethoxybenzamide"),
    ("O=S(=O)(N1CCOCC1)c1ccc(cc1)NC(=O)c1ccccc1", "a morpholinosulfonyl anilide"),
    ("CCOC(=O)c1ccc(cc1)N1CCN(CC1)CCOc1ccccc1", "a piperazinyl benzoate"),
    ("COc1ccccc1OCCN1CCC(CC1)NC(=O)c1ccccc1", "a methoxyphenoxyethyl piperidine"),
    ("COCCOc1cc2c(cc1OCCOC)ncnc2Nc1ccccc1", "a dimethoxyquinazoline aniline"),
    ("CN(C)CCOc1ccc(cc1)C(=O)Nc1ccc(cc1)OCc1ccccc1", "a dimethylaminoethoxy amide"),
    ("O=C(NCc1ccccc1)N1CCC(CC1)c1ccc(cc1)OCc1ccccc1", "a benzyl urea piperidine"),
    ("O=C(Nc1ccc(cc1)S(=O)(=O)N1CCOCC1)c1ccc(cc1)OCc1ccccc1", "a larger sulfonamide"),
    ("CCN1CCN(CC1)c1ccc(cc1)NC(=O)c1ccc(cc1)OCCc1ccccc1", "an ethylpiperazinyl amide"),
    # Larger still. Kinase inhibitors and similar run to 35-45 heavy
    # atoms and are much of what this audience will dock, so the
    # baseline has to reach them or it refuses on every one.
    ("COc1cc2nccc(Oc3ccc(NC(=O)C4(C(=O)Nc5ccccc5)CC4)cc3)c2cc1OC",
     "a large diaryl amide"),
    ("CN1CCN(CC1)CC(=O)Nc1ccc(cc1)Oc1ccnc2cc(OC)c(OC)cc12",
     "a methylpiperazinyl quinoline amide"),
    ("O=C(Nc1ccc(cc1)N1CCOCC1)c1ccc(cc1)Oc1ccc(cc1)C(F)(F)F",
     "a morpholino diaryl ether amide"),
    ("CCOc1cc2ncnc(Nc3cccc(c3)C(=O)NCCN3CCOCC3)c2cc1OCC",
     "a diethoxyquinazoline morpholine"),
    ("CN(C)CCCNC(=O)c1ccc(cc1)Nc1nccc(n1)-c1cccnc1",
     "an aminopyrimidine benzamide"),
    # Polar and acidic. Without these the library is all neutral and
    # greasy, and any polar or zwitterionic drug -- lisinopril,
    # methotrexate, most antibiotics -- gets refused for want of a
    # comparable baseline rather than judged.
    ("NC(CCC(=O)NC(CS)C(=O)NCC(O)=O)C(O)=O", "glutathione"),
    ("OC(=O)CN(CC(O)=O)CCN(CC(O)=O)CC(O)=O", "EDTA"),
    ("NC(Cc1ccc(O)cc1)C(=O)NC(CC(O)=O)C(=O)NCC(O)=O", "a short peptide"),
    ("OCC1OC(O)C(O)C(O)C1O", "glucose"),
    ("NC(=O)c1ccc(cn1)C(=O)NCCCC(N)C(O)=O", "a pyridine dicarboxamide"),
    ("OC(=O)C(N)CCCNC(N)=N", "arginine"),
    ("NC(Cc1c[nH]c2ccccc12)C(=O)NC(CCC(O)=O)C(O)=O", "a tryptophan dipeptide"),
]

# Known coverage limits, stated rather than discovered later: the
# library spans roughly 9-36 heavy atoms. Molecules much larger than
# that (ritonavir and most peptidomimetics) have no comparable decoy,
# and the control refuses rather than returning a flattering verdict.
# Extending it is a matter of adding rows, not changing code.


def pick_decoys(smiles: str, n: int = 8,
                tolerance: float = 0.30) -> list[tuple[str, str, dict]]:
    """
    Choose decoys matched to the user's molecule on the three
    properties that drive a docking score regardless of chemistry:
    heavy-atom count, lipophilicity and rotatable bonds.

    The matching matters. Compare a 40-atom compound against 15-atom
    decoys and it will win on size alone, which is precisely the
    illusion this control exists to dispel.
    """
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors

    query = Chem.MolFromSmiles(smiles)
    if query is None:
        return []
    q = (query.GetNumHeavyAtoms(), Crippen.MolLogP(query),
         Descriptors.NumRotatableBonds(query))

    scored = []
    for smi, name in DECOY_LIBRARY:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        d = (m.GetNumHeavyAtoms(), Crippen.MolLogP(m),
             Descriptors.NumRotatableBonds(m))
        # Relative distance across the three properties, with size
        # weighted heavily. Size is what produces the artefact this
        # control exists to catch -- a bigger molecule makes more
        # contacts and scores better whatever it is -- so a decoy of
        # the wrong size is useless however well its greasiness
        # matches. Unweighted, a rigid query such as celecoxib pulled
        # in small flexible-matched decoys and the comparison
        # collapsed.
        gap = (
            3.0 * abs(d[0] - q[0]) / max(q[0], 1)
            + abs(d[1] - q[1]) / max(abs(q[1]), 1.0)
            + 0.5 * abs(d[2] - q[2]) / max(q[2], 1)
        )
        scored.append((gap, smi, name,
                       {"heavy_atoms": d[0], "logp": round(d[1], 2),
                        "rotatable": d[2]}))

    scored.sort(key=lambda t: t[0])
    chosen = [(smi, name, props) for gap, smi, name, props in scored[:n]]
    return chosen


def match_quality(smiles: str,
                  decoys: list[tuple[str, str, dict]]) -> dict[str, Any]:
    """
    How well do these decoys actually match the query?

    This is a control on a control, and it exists because the first
    version of the decoy library did not have it. The library held only
    small molecules; matched against a real drug the closest decoy was
    a third smaller, the drug beat it on contact area alone, and the
    control reported a comfortable pass. A baseline that cannot be
    matched does not merely fail to help -- it certifies whatever it is
    given.
    """
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors

    q = Chem.MolFromSmiles(smiles)
    if q is None or not decoys:
        return {"usable": False, "why": "no decoys could be matched"}

    q_heavy = q.GetNumHeavyAtoms()
    q_logp = Crippen.MolLogP(q)
    q_rot = Descriptors.NumRotatableBonds(q)

    heavy = [d[2]["heavy_atoms"] for d in decoys]
    logp = [d[2]["logp"] for d in decoys]

    heavy_gap = abs(sum(heavy) / len(heavy) - q_heavy) / max(q_heavy, 1)
    logp_gap = abs(sum(logp) / len(logp) - q_logp)

    # 15% on size is about where contact area starts to decide the
    # comparison on its own.
    usable = heavy_gap <= 0.15 and logp_gap <= 1.5
    return {
        "usable": usable,
        "query_heavy_atoms": q_heavy,
        "query_logp": round(q_logp, 2),
        "query_rotatable": q_rot,
        "decoy_mean_heavy_atoms": round(sum(heavy) / len(heavy), 1),
        "decoy_mean_logp": round(sum(logp) / len(logp), 2),
        "size_mismatch": f"{heavy_gap:.0%}",
        "logp_mismatch": round(logp_gap, 2),
        "why": "" if usable else (
            f"the decoys average {sum(heavy) / len(heavy):.0f} heavy atoms "
            f"against your molecule's {q_heavy} - too different for the "
            "comparison to mean anything"
        ),
    }


def size_confound(results: list[dict]) -> dict[str, Any]:
    """
    Is this ranking a ranking of fit, or a ranking of molecular size?

    Vina's score becomes more negative as a molecule gets bigger and
    greasier, whether or not it fits anything. So when several
    compounds are docked together, the league table that comes out
    often reproduces their heavy-atom counts almost exactly -- and it
    looks, to anybody reading it, like structure-activity data.

    This is the check for that, and it is the one that most often
    changes what a result means. If the correlation between score and
    size is strong, the ordering carries little information about the
    pocket, and ligand efficiency (score divided by heavy atoms) is
    the honest way to compare.

    It is a different question from the decoy baseline. That asks
    whether one compound beats arbitrary molecules. This asks whether
    the ordering among the user's own compounds says anything at all.
    """
    try:
        from rdkit import Chem
    except ImportError:
        return {"usable": False, "why": "RDKit is not installed"}

    rows: list[tuple[str, float, int, float]] = []
    for r in results:
        if "score" not in r or not r.get("smiles"):
            continue
        mol = Chem.MolFromSmiles(r["smiles"])
        if mol is None:
            continue
        heavy = mol.GetNumHeavyAtoms()
        if heavy < 1:
            continue
        rows.append((r.get("name") or "molecule", float(r["score"]),
                     heavy, -float(r["score"]) / heavy))

    if len(rows) < 4:
        return {"usable": False,
                "why": "fewer than four molecules were docked"}

    n = len(rows)
    xs = [row[2] for row in rows]
    ys = [row[1] for row in rows]
    if len(set(xs)) < 2:
        return {"usable": False,
                "why": "every molecule is the same size"}

    mx, my = sum(xs) / n, sum(ys) / n
    sx = sum((x - mx) ** 2 for x in xs) ** 0.5
    sy = sum((y - my) ** 2 for y in ys) ** 0.5
    if sx == 0 or sy == 0:
        return {"usable": False, "why": "the scores are all the same"}
    r = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (sx * sy)

    by_score = sorted(rows, key=lambda t: t[1])
    by_efficiency = sorted(rows, key=lambda t: -t[3])
    spread = max(ys) - min(ys)

    return {
        "usable": True,
        "r": round(r, 2),
        "r_squared": round(r * r, 2),
        "molecules": n,
        "spread": round(spread, 2),
        "by_score": [f"{nm} {sc:.2f} ({h} atoms)"
                     for nm, sc, h, _ in by_score],
        "by_efficiency": [f"{nm} {le:.3f}" for nm, _, _, le in by_efficiency],
        "order_changes": [nm for nm, _, _, _ in by_score]
                         != [nm for nm, _, _, _ in by_efficiency],
    }


def hard_for_docking(results: list[dict]) -> str:
    """
    Is this the kind of molecule docking is worst at?

    Scoring functions were fitted to drug-sized, moderately greasy
    molecules in enclosed pockets. Two things take a compound outside
    where they work, and both are invisible in the output:

    SMALL. Below about fifteen heavy atoms a molecule is a fragment.
    It makes few contacts, the score has almost no dynamic range, and
    almost anything that fits scores about the same. Fragment work is
    done by crystallography and NMR precisely because docking cannot
    rank fragments.

    POLAR. Vina has no explicit desolvation term. Stripping the water
    off a polar molecule and off a polar pocket costs real energy, and
    that cost is simply not in the equation. For a compound with a
    negative logP the omission can be larger than every term that is
    there.

    This is said before the scores are read, not after, because a user
    who knows their compound is a small polar fragment will interpret
    the table correctly, and one who does not will publish it.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Crippen
    except ImportError:
        return ""

    small: list[str] = []
    polar: list[str] = []
    for item in results:
        if "score" not in item or not item.get("smiles"):
            continue
        mol = Chem.MolFromSmiles(item["smiles"])
        if mol is None:
            continue
        name = item.get("name") or "a molecule"
        if mol.GetNumHeavyAtoms() < 15:
            small.append(name)
        if Crippen.MolLogP(mol) < 0.0:
            polar.append(name)

    parts: list[str] = []
    if small:
        parts.append(
            f"{len(small)} of these are fragment-sized (under 15 heavy "
            "atoms). Docking scores have very little dynamic range for "
            "fragments: almost anything that fits the pocket scores "
            "about the same, which is why fragment screening is done "
            "by crystallography rather than by docking"
        )
    if polar:
        parts.append(
            f"{len(polar)} of these are polar (logP below zero). Vina "
            "has no desolvation term at all, so the energy it costs to "
            "strip water off a polar molecule and off a polar pocket "
            "is missing from the score entirely"
        )
    if not parts:
        return ""
    return (" These compounds are also near the edge of what any "
            "scoring function handles: " + "; and ".join(parts) + ".")


def size_verdict(confound: dict[str, Any],
                 difficulty: str = "") -> ControlOutcome:
    """Say what the size correlation means, in words that act."""
    if not confound.get("usable"):
        return ControlOutcome(
            name="Is this ranking about fit, or about size?",
            verdict=Verdict.NOT_APPLICABLE,
            plain_verdict=(
                "This check needs at least four molecules docked "
                f"together ({confound.get('why', '')})."
            ),
            what_to_do=(
                "Dock four or more compounds in one run if you want to "
                "compare them. Comparing scores from separate runs is "
                "not safe either way."
            ),
            detail={},
        )

    detail = {
        "molecules compared": confound["molecules"],
        "correlation of score with molecular size": confound["r"],
        "score-size linear association (R-squared)":
            f"{confound['r_squared'] * 100:.0f}%",
        "total spread in scores": f"{confound['spread']} kcal/mol",
        "score-spread caution threshold used here": "2.5 kcal/mol (heuristic)",
        "ranked by score": "; ".join(confound["by_score"]),
        "ranked by ligand efficiency (score per heavy atom)":
            "; ".join(confound["by_efficiency"]),
    }

    strong = confound["r_squared"] >= 0.5
    narrow = confound["spread"] < 2.5

    if strong:
        return ControlOutcome(
            name="Is this ranking about fit, or about size?",
            verdict=Verdict.FAIL if narrow else Verdict.WARN,
            plain_verdict=(
                "The scores are strongly associated with molecular size "
                f"in this set (R-squared {confound['r_squared']:.2f}). "
                "Size can influence contact-based scores. This "
                "correlation does not establish the cause of the ranking "
                "or distinguish true binding from a scoring bias."
                + (f" The whole table spans only {confound['spread']} "
                   "kcal/mol, below the 2.5 kcal/mol caution threshold "
                   "used here. Treat this small spread cautiously."
                   if narrow else "")
                + difficulty
            ),
            what_to_do=(
                "Do not present this ordering as structure-activity "
                "data. Two things are safe to say instead. First, give "
                "ligand efficiency alongside the score - it is the "
                "score divided by heavy atoms, and it is listed below; "
                "check whether it changes their order. "
                "Second, if you need to compare compounds, compare ones "
                "of similar size and review their poses. MD is an "
                "optional later exercise and does not resolve affinity "
                "ranking by itself."
            ),
            detail=detail,
            blocks_progress=True,
        )

    return ControlOutcome(
        name="Is this ranking about fit, or about size?",
        verdict=Verdict.PASS,
        plain_verdict=(
            "This set did not cross the size-association warning "
            f"threshold (R-squared {confound['r_squared']:.2f}). This "
            "does not show that the ordering measures pocket fit or "
            "binding affinity."
            + (f" Be careful even so: the whole table spans "
               f"{confound['spread']} kcal/mol, below the 2.5 kcal/mol "
               "caution threshold used here."
               if narrow else "")
            + difficulty
        ),
        what_to_do=(
            "Report the scores and ligand efficiencies descriptively. "
            "Review pose plausibility and the other checks; this "
            "diagnostic does not validate a potency ranking."
        ),
        detail=detail,
    )


def decoy_verdict(user_score: float, decoy_scores: list[float],
                  user_name: str = "your compound",
                  match: dict[str, Any] | None = None) -> ControlOutcome:
    """
    Did the user's molecule beat arbitrary molecules of the same size?

    Deliberately blunt. A compound that does not clear the best decoy
    has not been shown to do anything, and the honest way to say that
    is to say it.
    """
    if not decoy_scores:
        return ControlOutcome(
            name="Does it beat an arbitrary molecule?",
            verdict=Verdict.NOT_APPLICABLE,
            plain_verdict=(
                "No baseline could be run, so there is nothing to compare "
                "this score against. A docking score on its own means very "
                "little."
            ),
            what_to_do=(
                "Check the calculation log and comparison-molecule "
                "availability. If the issue can be fixed, rerun; "
                "otherwise record that no baseline was available."
            ),
            detail={},
        )

    if match is not None and not match.get("usable", True):
        return ControlOutcome(
            name="Does it beat an arbitrary molecule?",
            verdict=Verdict.NOT_APPLICABLE,
            plain_verdict=(
                "This comparison cannot be made: "
                + match.get("why", "the decoys do not match your molecule")
                + ". Poorly matched molecules can confound this "
                "comparison, so no verdict is given."
            ),
            what_to_do=(
                "This usually happens when your molecule is very large or "
                "very small. Review the available baseline library. "
                "Adding your own comparison compounds can be a separate "
                "exercise, but does not change this automatic baseline."
            ),
            detail=match,
            blocks_progress=True,
        )

    best_decoy = min(decoy_scores)          # more negative is better
    mean_decoy = sum(decoy_scores) / len(decoy_scores)
    margin = round(best_decoy - user_score, 2)      # positive = user better

    detail = {
        f"{user_name} score": user_score,
        "best decoy": round(best_decoy, 2),
        "average decoy": round(mean_decoy, 2),
        "decoys run": len(decoy_scores),
        "margin over best decoy": margin,
    }

    if user_score > best_decoy:
        return ControlOutcome(
            name="Does it beat an arbitrary molecule?",
            verdict=Verdict.FAIL,
            plain_verdict=(
                f"Score comparison needs review. {user_name} scored {user_score}, "
                f"while the best matched comparison scored {best_decoy:.2f}. "
                "This run did not favour the leading drug over that comparison. "
                "This does not overturn experimental evidence for a known drug and target. "
                "The comparison molecules are not proven nonbinders."
            ),
            what_to_do=(
                "Report the comparison and inspect both poses, preparation and search box. "
                "For a known reference drug, use this result to review the method. "
                "A new hit claim needs additional evidence. The score alone does not establish binding strength."
            ),
            detail=detail,
            blocks_progress=True,
        )

    if margin < 1.0:
        return ControlOutcome(
            name="Does it beat an arbitrary molecule?",
            verdict=Verdict.WARN,
            plain_verdict=(
                f"Only just. {user_name} scored {user_score} against "
                f"{best_decoy:.2f} for an arbitrary molecule of similar "
                f"size - a margin of {margin} kcal/mol. This is below "
                "the 1.0 kcal/mol separation threshold used in this "
                "exercise and is not a reliable potency difference."
            ),
            what_to_do=(
                "Treat this as no clear separation, not as a hit. In your "
                "write-up say the compound did not separate from the "
                "baseline. If you want to take this molecule further, "
                "review the poses, preparation and comparison set. "
                "Crossing a score threshold alone does not establish "
                "improved binding."
            ),
            detail=detail,
            blocks_progress=True,
        )

    return ControlOutcome(
        name="Does it beat an arbitrary molecule?",
        verdict=Verdict.PASS,
        plain_verdict=(
            f"{user_name} scored {user_score}, beating the best arbitrary "
            f"molecule of similar size by {margin} kcal/mol. This "
            "crosses the exercise's separation threshold for this small "
            "comparison set. It does not establish binding or exclude "
            "size-related scoring effects."
        ),
        what_to_do=(
            "This is a good result and you may report it - but report it "
            "in the right words. Correct: 'docked favourably, and "
            f"outperformed size-matched baseline molecules by {margin} "
            "kcal/mol'. Not correct: 'binds strongly' or 'is a potent "
            "inhibitor'. Docking cannot show either of those. "
            "Experimental testing is needed to assess binding or "
            "activity; MD is an optional additional modelling exercise."
        ),
        detail=detail,
    )


def crystal_pose(prepared_pdb: Path, label: str,
                 template_smiles: str = "") -> "Any":
    """
    The bound molecule, where the experiment actually put it.

    Pulled straight out of the prepared receptor rather than
    downloaded, so it needs no network and is guaranteed to be the
    same coordinates the docking box was built from.

    Bond orders are not in a PDB file, so when the SMILES of the same
    compound is available it is used as a template. Without that step
    the two molecules being compared can end up with different bond
    graphs, and the RMSD either fails outright or is computed over a
    partial match, which is worse because it returns a number.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    try:
        text = prepared_pdb.read_text(errors="ignore")
    except OSError:
        return None

    resname, _, chain_seq = label.partition(" ")
    want_chain = chain_seq[:1] if chain_seq else ""
    want_seq = chain_seq[1:].strip() if chain_seq else ""

    lines = []
    for line in text.splitlines():
        if not line.startswith("HETATM"):
            continue
        if line[17:20].strip().upper() != resname.upper():
            continue
        if want_chain and line[21] != want_chain:
            continue
        if want_seq and line[22:27].strip() != want_seq:
            continue
        if (line[76:78].strip().upper() or "C") == "H":
            continue
        lines.append(line)
    if len(lines) < 4:
        return None

    raw = Chem.MolFromPDBBlock("\n".join(lines) + "\nEND\n",
                               removeHs=True, sanitize=False)
    if raw is None:
        return None
    if not template_smiles:
        return raw
    template = Chem.MolFromSmiles(template_smiles)
    if template is None:
        return raw
    try:
        return AllChem.AssignBondOrdersFromTemplate(template, raw)
    except Exception:  # noqa: BLE001
        return raw


def same_compound(mol_a: "Any", mol_b: "Any") -> bool:
    """Do two molecules have the same heavy-atom composition?"""
    if mol_a is None or mol_b is None:
        return False
    def formula(m) -> dict:
        counts: dict[str, int] = {}
        for atom in m.GetAtoms():
            if atom.GetAtomicNum() == 1:
                continue
            sym = atom.GetSymbol()
            counts[sym] = counts.get(sym, 0) + 1
        return counts
    return formula(mol_a) == formula(mol_b)


def redock_rmsd(prepared_pdb: Path, label: str,
                results: list[dict]) -> tuple[float | None, str, str]:
    """
    Did the run put the known molecule back where it really sits?

    No extra docking is needed when the user has already docked the
    compound that was bound in the structure -- which, for anybody
    working on a drug whose complex has been solved, is the usual
    case. The pose is there; it only has to be compared.

    Returns the RMSD, the name of the molecule compared, and a
    sentence explaining any failure to compare.
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    if not label:
        return None, "", ("there is no bound molecule in this structure "
                          "to compare against")

    bare = crystal_pose(prepared_pdb, label)
    if bare is None:
        return None, "", (f"the bound molecule {label} could not be read "
                          "back out of the prepared structure")

    # Why each candidate was rejected, so a failure is explained
    # rather than silently swallowed. An earlier version of this loop
    # used bare `continue` throughout, and when the comparison failed
    # the user was told no matching molecule had been docked -- which
    # was not true and sent them looking in the wrong place.
    trouble: list[str] = []

    for item in results:
        if "score" not in item or not item.get("file"):
            continue
        smiles = item.get("smiles", "")
        docked_mol = Chem.MolFromSmiles(smiles) if smiles else None
        if not same_compound(bare, docked_mol):
            continue

        name = item.get("name") or label
        reference = crystal_pose(prepared_pdb, label, smiles) or bare
        try:
            from stages.md import pdbqt_pose_to_pdb
            block = pdbqt_pose_to_pdb(
                Path(item["file"]).read_text(errors="ignore"), model=1)
        except Exception as exc:  # noqa: BLE001
            trouble.append(f"{name}: the docked pose file could not be "
                           f"read ({exc})")
            continue

        pose = Chem.MolFromPDBBlock(block, removeHs=True, sanitize=False)
        if pose is None:
            trouble.append(f"{name}: the docked pose could not be turned "
                           "back into a molecule")
            continue
        if docked_mol is not None:
            try:
                from rdkit.Chem import AllChem
                pose = AllChem.AssignBondOrdersFromTemplate(docked_mol, pose)
            except Exception:  # noqa: BLE001
                pass    # the comparison below may still succeed
        try:
            # CalcRMS, never GetBestRMS. See best_rmsd above: the
            # second one superimposes first and would pass a pose on
            # the far side of the protein.
            return (round(rdMolAlign.CalcRMS(pose, reference), 3),
                    name, "")
        except Exception as exc:  # noqa: BLE001
            trouble.append(f"{name}: the docked pose and the bound one "
                           f"could not be matched up ({exc})")

    if trouble:
        return None, "", "; ".join(trouble[:3])

    return None, "", (
        f"none of the molecules you docked is {label.split()[0]}, the one "
        "bound in this structure, so there is no known pose to compare "
        "against"
    )


def redock_verdict(rmsd: float | None, ligand_name: str = "") -> ControlOutcome:
    """
    Can this setup reproduce a pose that is already known?

    This tests pose reproduction for one reference ligand under the
    recorded settings; passing does not validate affinity predictions.
    """
    if rmsd is None:
        # WARN, not NOT_APPLICABLE, and this distinction matters more
        # than it looks. NOT_APPLICABLE ranks as harmless, so a run
        # that could not be validated at all would show whichever
        # other check happened to pass -- and a user would read a
        # green verdict on a docking run with nothing anchoring it.
        # "We could not check this" is a warning about the result, not
        # an absence of one.
        return ControlOutcome(
            name="Can this setup reproduce a known pose?",
            verdict=Verdict.WARN,
            plain_verdict=(
                "Reference-pose reproduction could not be checked. A "
                "usable bound ligand and a matching docked pose are "
                "needed for this test. Interpret the scores with this "
                "missing check clearly recorded."
            ),
            what_to_do=(
                "Go back to the Structure step and look through the "
                "list for an entry with a molecule bound inside it. "
                "A suitable complex provides a reference pose for this "
                "test; also review its local structure quality.\n\n"
                "Check the Binding site step too. An automatically "
                "detected cavity is not proof of a biological binding "
                "site; a manually chosen box also needs evidence.\n\n"
                "If this protein genuinely has no complex in the "
                "database, you can continue - but say plainly in your "
                "write-up that the docking protocol was not validated."
            ),
            detail={},
            blocks_progress=True,
        )

    detail = {"RMSD to the crystal pose": f"{rmsd} A",
              "pass mark": f"{RMSD_PASS} A",
              "ligand re-docked": ligand_name}

    if rmsd <= RMSD_PASS:
        return ControlOutcome(
            name="Can this setup reproduce a known pose?",
            verdict=Verdict.PASS,
            plain_verdict=(
                f"Yes. The known ligand was re-docked to within {rmsd} "
                "angstroms of its experimental reference pose. This "
                f"meets this exercise's {RMSD_PASS} angstrom threshold "
                "for pose reproduction. It supports this setup for "
                "this ligand, but does not validate other poses, "
                "binding affinities or molecule rankings."
            ),
            what_to_do=(
                "Record this result and the settings in your write-up: "
                "'the protocol reproduced the crystallographic "
                f"pose of {ligand_name or 'the bound ligand'} with an RMSD "
                f"of {rmsd} angstroms'. Review the comparison and other "
                "checks before interpreting new compounds."
            ),
            detail=detail,
        )

    if rmsd <= 3.0:
        return ControlOutcome(
            name="Can this setup reproduce a known pose?",
            verdict=Verdict.WARN,
            plain_verdict=(
                f"Not quite. The known ligand came back {rmsd} angstroms "
                f"from its real position - over the usual {RMSD_PASS} "
                "angstrom threshold used in this exercise. The "
                "reference-pose test has not passed."
            ),
            what_to_do=(
                "Review the reference and predicted poses, preparation "
                "and search box before interpreting other results. "
                "Report the RMSD and this warning. The calculation "
                "alone does not establish binding in this region or "
                "specific experimental contacts."
            ),
            detail=detail,
            blocks_progress=True,
        )

    return ControlOutcome(
        name="Can this setup reproduce a known pose?",
        verdict=Verdict.FAIL,
        plain_verdict=(
            f"No. The known ligand was placed {rmsd} angstroms from "
            "the experimental reference pose. It failed this "
            "pose-reproduction test, so confidence in other poses from "
            "this setup is reduced. Investigate the cause."
        ),
        what_to_do=(
            "Stop and fix this before reading any score below. Check in "
            "the following possible causes: "
            "(1) the box - is it centred on the real pocket and not too "
            "large? (2) the protein - was a metal or cofactor removed at "
            "the Prepare step? (3) holes - is there a gap in the chain "
            "beside the pocket? Fix the cause and run again. If nothing "
            "helps, this structure is not suitable; choose another one."
        ),
        detail=detail,
        blocks_progress=True,
    )


# ----------------------------------------------------------------------
# The stage
# ----------------------------------------------------------------------

class DockStage(Stage):
    """
    Docking, with both controls run before any result is shown.

    The order matters and is not configurable. Re-docking first: if the
    setup cannot reproduce a pose that is already known, nothing below
    it is worth reading. Then the user's molecules. Then the decoy
    baseline, which is what turns a score into a comparison.
    """

    stage_id = "06_dock"
    title = "Dock"
    tier = Tier.BUNDLED
    requires = ["03_receptor", "04_site", "05_ligands"]

    def check(self) -> CheckResult:
        if find_vina() is None:
            return CheckResult(
                readiness=Readiness.BLOCKED,
                reason="AutoDock Vina is not installed.",
                fix=how_to_install("vina"),
            )
        for needed, what in (("03_receptor", "the prepared protein"),
                             ("04_site", "the binding site"),
                             ("05_ligands", "your molecules")):
            if not self.project.is_complete(needed):
                return CheckResult(
                    readiness=Readiness.NEEDS_INPUT,
                    reason=f"{what.capitalize()} is not ready yet.",
                    fix=f"Complete that step first.",
                )
        return CheckResult(readiness=Readiness.READY)

    def explain(self) -> Teaching:
        return Teaching(
            button="Dock your molecules, and check whether the result means "
                   "anything.",
            panel=(
                "Docking puts your molecule into the pocket in many "
                "orientations, scores each one, and reports the best. The "
                "score looks like an energy and is quoted in kilocalories "
                "per mole, which makes it feel like a measurement. It is "
                "not one.\n\n"
                "It is a number from a function fitted to known "
                "structures, and it is useful for ranking poses of one "
                "molecule in one pocket. It is much weaker at ranking "
                "different molecules against each other, and it is not "
                "comparable at all between different proteins. A score of "
                "-9 against your target and -9 against another tells you "
                "nothing about selectivity.\n\n"
                "Because the number alone is so easy to over-read, "
                "checks are attempted automatically and their outcomes "
                "are shown with the results. The first compares a "
                "re-docked ligand with an experimental reference whose "
                "position is known, to find out whether this setup can "
                "reproduce an answer we already have. The second scores a "
                "set of arbitrary molecules matched to yours for size and "
                "greasiness, to find out whether your compound is doing "
                "better than something chosen for having nothing to do "
                "with the target. This small baseline is a diagnostic, "
                "not a validated set of confirmed non-binders. A "
                "size-related ranking check is also reported when there "
                "are enough molecules."
            ),
            simple=(
                "This is the main step. The computer takes your "
                "molecule, puts it into the protein pocket in thousands "
                "of different positions and angles, gives each one a "
                "score, and shows you the best.\n\n"
                "The score looks like this: -8.5 kcal/mol. More negative "
                "means more favourable under this scoring method. It "
                "does not directly measure physical fit or potency.\n\n"
                "IMPORTANT: this score is NOT a measurement of binding "
                "strength. It is a rough estimate from a formula. You "
                "cannot say 'my drug binds at -8.5 so it is potent'. "
                "That sentence appears in many published papers and it "
                "is wrong.\n\n"
                "The software checks the reference pose and a small "
                "comparison baseline. Where enough molecules are "
                "available, it also flags size-related rankings.\n\n"
                "CHECK 1 - Can we trust this setup?\n"
                "We take a drug whose real position inside this pocket "
                "is already known from the laboratory. We hide that "
                "position and ask the computer to find it again. If the "
                "computer cannot find an answer we already know, it "
                "has not passed this validation test. Investigate the "
                "setup before interpreting other poses.\n\n"
                "CHECK 2 - Is your score actually good?\n"
                "We also dock a few ordinary molecules that have nothing "
                "to do with your protein, chosen to be the same size and "
                "greasiness as the leading compound. Failure to separate "
                "from this baseline weakens the score-based argument. "
                "A few comparison molecules do not establish binding "
                "or prove that size is the only cause of the score."
            ),
            not_mean=[
                "A good score is not evidence of binding. It is a "
                "hypothesis for an experiment, not confirmation of activity.",
                "Scores from different proteins cannot be compared. The "
                "function is not calibrated across targets, so -9 here and "
                "-7 there says nothing about selectivity.",
                "A difference of less than about one kilocalorie per mole "
                "should not be treated as a reliable potency difference. "
                "Small score differences may not support a stable ranking.",
                "Docking assumes a rigid protein. If your pocket opens or "
                "reshapes to accept a ligand - which many do - the result "
                "can be wrong in a way no score will reveal.",
            ],
            learn_more=[
                ("AutoDock Vina - what the score is",
                 "https://autodock-vina.readthedocs.io/en/latest/"),
                ("Protein preparation and docking protocols",
                 "https://doi.org/10.1038/nprot.2016.051"),
            ],
        )

    # -- helpers -----------------------------------------------------------

    def _receptor(self, folder: Path, exclude: str) -> tuple[Path | None, str]:
        from core.receptor_pdbqt import convert

        step = self.project.get_step("03_receptor")
        if not step:
            return None, "The protein has not been prepared."
        values = json.loads(step["values_json"] or "{}")
        prepared = Path(values.get("prepared_file", ""))
        if not prepared.exists():
            return None, "The prepared protein file is missing."

        # The molecule being docked must not be in the receptor.
        keep = [k["label"] for k in values.get("kept", [])
                if k["label"] != exclude]
        out = folder / "receptor.pdbqt"
        result = convert(prepared, out, keep_hetero=keep)
        if not result.ok:
            return None, f"{result.problem} {result.suggestion}"
        return out, ""

    @staticmethod
    def _ligand_pdbqt(smiles: str, folder: Path, name: str,
                      seed: int = DEFAULT_SEED, prepared_file: str = "") -> tuple[Path | None, str]:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = (next(iter(Chem.SDMolSupplier(prepared_file, removeHs=False)), None)
               if prepared_file else Chem.MolFromSmiles(smiles))
        if mol is None:
            return None, f"{name}: the structure could not be read."
        if prepared_file:
            expected = Chem.MolFromSmiles(smiles)
            if expected is None or Chem.MolToSmiles(Chem.RemoveHs(mol)) != Chem.MolToSmiles(expected):
                return None, f"{name}: the prepared structure and recorded identity differ. Prepare the drug again."
            if not mol.GetNumConformers():
                return None, f"{name}: the prepared drug has no 3D coordinates."
            mol = Chem.AddHs(mol, addCoords=True)
        else:
            mol = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = seed
            if AllChem.EmbedMolecule(mol, params) != 0:
                return None, f"{name}: no 3D shape could be built."
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            except Exception:  # noqa: BLE001
                pass
        setups = MoleculePreparation().prepare(mol)
        text, ok, err = PDBQTWriterLegacy.write_string(setups[0])
        if not ok:
            return None, f"{name}: {err}"
        safe = "".join(c if c.isalnum() else "_" for c in name) or "ligand"
        path = folder / f"{safe}.pdbqt"
        path.write_text(text)
        return path, ""

    # -- run ---------------------------------------------------------------

    def run(self, exhaustiveness: int = DEFAULT_EXHAUSTIVENESS,
            n_decoys: int = 3, **kwargs: Any) -> StageResult:
        folder = self.project.step_folder(self.stage_id)

        site = self.project.get_step("04_site")
        site_values = json.loads(site["values_json"] or "{}") if site else {}
        box = site_values.get("box")
        box_source = (box or {}).get("source", "")
        if not box:
            return StageResult(success=False,
                               message="No docking box has been defined.")

        ligands_step = self.project.get_step("05_ligands")
        prepared = json.loads(
            ligands_step["values_json"] or "{}").get("prepared", [])
        if not prepared:
            return StageResult(success=False,
                               message="No prepared molecules were found.")

        receptor_ligand = json.loads(
            self.project.get_step("03_receptor")["values_json"] or "{}"
        ).get("ligand", "")

        receptor, why = self._receptor(folder, exclude=receptor_ligand)
        if receptor is None:
            return StageResult(success=False, message=why)

        results: list[dict[str, Any]] = []
        for item in prepared:
            lig, why = self._ligand_pdbqt(item["smiles"], folder,
                                          item.get("name", "molecule"), prepared_file=item.get("file", ""))
            if lig is None:
                results.append({"name": item.get("name", ""), "error": why})
                continue
            out = folder / f"{lig.stem}_docked.pdbqt"
            dr = run_vina(receptor, lig, box, out,
                          exhaustiveness=exhaustiveness)
            if not dr.ok or dr.best is None:
                results.append({"name": item.get("name", ""),
                                "error": dr.error})
                continue
            results.append({
                "name": item.get("name", ""),
                "smiles": item["smiles"],
                "score": dr.best.score,
                "poses": [{"rank": p.rank, "score": p.score}
                          for p in dr.poses[:5]],
                "file": str(out),
                "input_file": str(lig),
            })

        scored = [r for r in results if "score" in r]
        if not scored:
            return StageResult(
                success=False,
                message="None of the molecules could be docked. "
                        + (results[0].get("error", "") if results else ""),
            )

        best = min(scored, key=lambda r: r["score"])
        control = self._baseline(best, folder, receptor, box,
                                 exhaustiveness, n_decoys)

        # The check that matters most, and the one this stage went
        # without for far too long: can the setup reproduce a pose
        # that is already known? Everything else is unverified until
        # this has been asked.
        prepared_step = self.project.get_step("03_receptor")
        prepared_file = Path(json.loads(
            prepared_step["values_json"] or "{}").get("prepared_file", "")
        ) if prepared_step else Path("")
        rmsd, redocked_name, why_not = redock_rmsd(
            prepared_file, receptor_ligand, scored)
        redock = redock_verdict(rmsd, redocked_name)
        if rmsd is None and why_not:
            redock.plain_verdict = (
                "The docking setup could not be validated, because "
                + why_not + "."
            )
            redock.what_to_do = (
                "Dock the molecule that is bound in this structure "
                f"({receptor_ligand.split()[0]}) alongside yours. The "
                "software will then remove it, see whether the docking "
                "puts it back where it really sits, and tell you how "
                "close it got. Until then, reference-pose reproduction "
                "has not been checked."
                if receptor_ligand else
                "Go back to the Structure step and choose an entry with "
                "a molecule bound in it if a suitable complex is "
                "available. Without one, this reference-pose test "
                "cannot be performed."
            )

        # The third question, asked only when there is a league table
        # to ask it about: is the ordering among the user's own
        # compounds about fit, or about how big they are?
        confound = size_confound(scored)
        difficulty = hard_for_docking(scored)
        size_check = (size_verdict(confound, difficulty)
                      if confound.get("usable") else None)

        rank = {Verdict.PASS: 0, Verdict.NOT_APPLICABLE: 0,
                Verdict.WARN: 1, Verdict.FAIL: 2}

        # Three checks now run. The user sees the most serious one in
        # full, and every one of them in the detail underneath -- a
        # check whose result is not shown may as well not have run.
        # Three checks run. Whichever is most serious leads, and all
        # three always appear underneath.
        #
        # This used to be done by merging one outcome into another in
        # sequence, and it silently lost results: each merge rebuilt
        # the detail from the incoming check, so the re-docking line
        # added by the first merge was thrown away by the second. The
        # user replaced the file, re-ran, and still saw no RMSD. A
        # check whose result does not reach the screen has not run as
        # far as anybody is concerned, so they are collected once,
        # here, and nothing is threaded through anything.
        checks = [("baseline", control), ("re-docking", redock)]
        if size_check is not None:
            checks.append(("size", size_check))

        leader_name, leader = max(
            checks, key=lambda pair: rank[pair[1].verdict])

        summaries = {
            f"{name} check": outcome.plain_verdict
            for name, outcome in checks if name != leader_name
        }
        control = ControlOutcome(
            name=leader.name,
            verdict=leader.verdict,
            plain_verdict=leader.plain_verdict,
            what_to_do=leader.what_to_do,
            detail={**leader.detail, **summaries},
            blocks_progress=any(o.blocks_progress for _, o in checks),
        )

        self.record(
            ProvenanceEntry(
                stage=self.stage_id,
                tool="AutoDock Vina",
                version=self._vina_version(),
                parameters={"exhaustiveness": exhaustiveness,
                            "box": box, "molecules": len(scored),
                            "receptor_ligand_excluded": receptor_ligand},
                seed=DEFAULT_SEED,
                citation="Eberhardt J et al. AutoDock Vina 1.2.0. "
                         "J Chem Inf Model. 2021;61:3891-3898.",
            )
        )
        for e in self.provenance:
            self.project.add_provenance(e)

        return StageResult(
            success=True,
            summary=(f"Docked {len(scored)} molecule(s); best "
                     f"{best['name']} at {best['score']} kcal/mol"),
            values={"results": results, "box": box,
                    "receptor": str(receptor),
                    "settings": {"exhaustiveness": exhaustiveness, "seed": DEFAULT_SEED,
                                 "n_poses": DEFAULT_POSES, "vina_version": self._vina_version()},
                    "box_source": box_source,
                    "redock_rmsd": rmsd,
                    "redock_molecule": redocked_name,
                    "checks": {name: asdict(outcome) for name, outcome in checks},
                    "redock_note": why_not,
                    "size_confound": confound},
            artefacts=[Path(r["file"]) for r in scored if r.get("file")],
            control=control,
        )

    def _baseline(self, best: dict, folder: Path, receptor: Path,
                  box: dict, exhaustiveness: int,
                  n_decoys: int) -> ControlOutcome:
        """The decoy comparison -- what turns a score into a finding."""
        decoys = pick_decoys(best["smiles"], n=n_decoys)
        quality = match_quality(best["smiles"], decoys)
        scores: list[float] = []
        for smi, name, _ in decoys:
            lig, why = self._ligand_pdbqt(smi, folder, f"decoy_{name}")
            if lig is None:
                continue
            out = folder / f"{lig.stem}_docked.pdbqt"
            dr = run_vina(receptor, lig, box, out,
                          exhaustiveness=exhaustiveness, n_poses=3)
            if dr.ok and dr.best:
                scores.append(dr.best.score)
        return decoy_verdict(best["score"], scores, best["name"],
                             match=quality)

    @staticmethod
    def _vina_version() -> str:
        exe = find_vina()
        if not exe:
            return "unknown"
        try:
            out = subprocess.run([exe, "--version"], capture_output=True,
                                 text=True, timeout=20)
            return (out.stdout or out.stderr).strip().splitlines()[0][:60]
        except Exception:  # noqa: BLE001
            return "installed"

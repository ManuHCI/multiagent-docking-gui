"""
Stage 5 -- Prepare the ligands.

Takes what the user typed or drew and turns it into the 3D structure
that docking will actually use. The gap between those two things is
wider than most newcomers expect, and it is a common source of silent
error: the molecule that gets docked is not the molecule the user
thought they drew.

So the control for this stage is a change list. Every modification is
named, in plain English, and the user has to look at it before going
on. Specification Stage 5.

Runs entirely offline. RDKit only.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, Draw
from rdkit.Chem.MolStandardize import rdMolStandardize

from core.install import how_to_install
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

RDLogger.DisableLog("rdApp.*")   # RDKit's warnings are not for this audience


# ----------------------------------------------------------------------
# What happened to one molecule
# ----------------------------------------------------------------------

@dataclass
class Prepared:
    name: str
    input_smiles: str
    smiles: str = ""              # after preparation
    changes: list[str] = field(default_factory=list)
    properties: dict[str, float] = field(default_factory=dict)
    sdf_path: Path | None = None
    error: str = ""

    @property
    def ok(self) -> bool:
        return not self.error

    @property
    def changed(self) -> bool:
        return bool(self.changes)


# ----------------------------------------------------------------------
# Preparation, one molecule at a time
# ----------------------------------------------------------------------

def prepare_one(
    smiles: str,
    name: str = "",
    ph_neutral: bool = True,
    canonical_tautomer: bool = False,
    seed: int = 0xF00D,
) -> Prepared:
    """
    Standardise, desalt, neutralise, then build a 3D conformer.

    Every step that alters the molecule appends a sentence to
    `changes`, worded for a reader who does not know what a
    "fragment parent" is.
    """
    out = Prepared(name=name or "Molecule", input_smiles=smiles)

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        out.error = (
            "That is not a structure this software can read. Check the "
            "SMILES text, or draw the molecule instead."
        )
        return out

    before = Chem.MolToSmiles(mol)
    n_frags_before = len(Chem.GetMolFrags(mol))

    # 1. clean up -- sanitise, fix common drawing problems
    try:
        mol = rdMolStandardize.Cleanup(mol)
    except Exception:  # noqa: BLE001
        out.error = "The structure could not be tidied up. It may be invalid."
        return out
    if Chem.MolToSmiles(mol) != before:
        out.changes.append(
            "Tidied the structure (standard bond and charge conventions)."
        )

    # 2. keep only the main molecule -- drop salts and counter-ions
    parent = rdMolStandardize.FragmentParent(mol)
    if len(Chem.GetMolFrags(mol)) > 1 or n_frags_before > 1:
        out.changes.append(
            "Removed the salt or counter-ion and kept only the main "
            "molecule. Check that this is the chemical entity you intend "
            "to dock."
        )
    mol = parent

    # 3. neutral form
    if ph_neutral:
        uncharged = rdMolStandardize.Uncharger().uncharge(mol)
        if Chem.MolToSmiles(uncharged) != Chem.MolToSmiles(mol):
            out.changes.append(
                "Neutralised charges where possible under this workflow's "
                "standardisation convention. This does not predict the "
                "ionisation state at pH 7.4, which can matter for binding."
            )
        mol = uncharged

    # 4. tautomer -- OFF by default, deliberately.
    #
    # RDKit's canonical tautomer is chosen by a general scoring function
    # that knows nothing about pharmacology, and on drug-like molecules
    # it can move a hydrogen that matters. Amlodipine is the worked
    # example: canonicalising it shifts the hydrogen off the ring
    # nitrogen, turning the 1,4-dihydropyridine into a pyridine and
    # destroying the NH donor that defines the entire calcium-channel-
    # blocker class. The docking result would look perfectly normal and
    # be meaningless.
    #
    # The user drew or quoted a particular form, usually the documented
    # one. Keep it. Anyone who wants the canonical form can ask.
    if canonical_tautomer:
        try:
            canon = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)
            if Chem.MolToSmiles(canon) != Chem.MolToSmiles(mol):
                out.changes.append(
                    "Chose a different tautomer - a hydrogen was moved to "
                    "another position. Check that the molecule is still "
                    "the form you meant; this can matter a great deal."
                )
            mol = canon
        except Exception:  # noqa: BLE001
            pass   # tautomer search can fail on exotic structures; not fatal

    # RDKit's standardisation operations return molecules whose ring
    # information has not been rebuilt, and the descriptor functions
    # below need it. Without this the stage crashes on any molecule
    # that was desalted or neutralised -- which is most real drugs.
    try:
        Chem.SanitizeMol(mol)
    except Exception:  # noqa: BLE001
        out.error = (
            "The prepared structure is not chemically valid. Check what "
            "you entered."
        )
        return out

    out.smiles = Chem.MolToSmiles(mol)

    # 5. properties, for the applicability check later
    out.properties = {
        "Molecular weight": round(Descriptors.MolWt(mol), 2),
        "LogP": round(Descriptors.MolLogP(mol), 2),
        "Polar surface area": round(Descriptors.TPSA(mol), 2),
        "H-bond donors": Descriptors.NumHDonors(mol),
        "H-bond acceptors": Descriptors.NumHAcceptors(mol),
        "Rotatable bonds": Descriptors.NumRotatableBonds(mol),
        "Heavy atoms": mol.GetNumHeavyAtoms(),
    }

    # 6. 3D structure
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed          # recorded, so the run repeats exactly
    if AllChem.EmbedMolecule(mol_h, params) != 0:
        out.error = (
            "A three-dimensional shape could not be built for this "
            "molecule. It may be too strained or unusual. Try a simpler "
            "related structure."
        )
        return out
    try:
        AllChem.MMFFOptimizeMolecule(mol_h, maxIters=500)
    except Exception:  # noqa: BLE001
        AllChem.UFFOptimizeMolecule(mol_h, maxIters=500)

    out._mol3d = mol_h  # type: ignore[attr-defined]
    return out


# ----------------------------------------------------------------------
# The stage
# ----------------------------------------------------------------------

class LigandStage(Stage):
    stage_id = "05_ligands"
    title = "Prepare the ligands"
    tier = Tier.BUNDLED
    requires: list[str] = []      # can be done before or after the protein

    SEED = 0xF00D

    # -- contract ---------------------------------------------------------

    def check(self) -> CheckResult:
        try:
            import rdkit  # noqa: F401
        except ImportError:
            return CheckResult(
                readiness=Readiness.BLOCKED,
                reason="The molecule tools are not installed on this "
                       "computer.",
                fix="Open Anaconda Prompt and run:  "
                    + how_to_install("rdkit"),
            )
        return CheckResult(readiness=Readiness.READY)

    def explain(self) -> Teaching:
        return Teaching(
            rule=(
                "Copy the SMILES from PubChem or ChEMBL - never type it. Then read what preparation changed, and check the molecule is big enough to be scored."
            ),
            rule_detail=[
                "Search the drug on PubChem, open its page, copy the Canonical SMILES. Typing it by hand is how an ethoxy becomes a methoxyethoxy and you dock a compound that does not exist.",
                "Small molecules can be difficult to rank with docking scores. Levetiracetam has 12 heavy atoms, paracetamol 11 and aspirin 13. Review their poses and avoid treating score differences as potency differences.",
                "Read the change list. Desalting, neutralising and tautomer choice can alter the molecule you thought you were docking.",
                'Give every molecule a name, or the results all come back called "Molecule".',
            ],
            button="Turn your molecules into the 3D structures docking needs.",
            panel=(
                "A molecule you draw on paper is a flat diagram. Docking "
                "needs a real three-dimensional object with every hydrogen "
                "in place, so this step builds one.\n\n"
                "On the way it makes several decisions that change the "
                "molecule. If you type the formula for a drug as it is "
                "supplied - amlodipine besylate, say - it contains the "
                "besylate counter-ion, which is not what binds the target. "
                "That gets removed. Charges are neutralised where possible "
                "under this workflow's standardisation convention; this "
                "does not determine the form at physiological pH.\n\n"
                "One thing this software deliberately does not do is move "
                "hydrogens around to pick a 'preferred' tautomer. Automatic "
                "tautomer tools are scored on general chemistry and know "
                "nothing about your binding site. A different tautomer "
                "can alter hydrogen donors, acceptors and bonding. The "
                "form you entered is kept unless standardisation needs "
                "to correct its representation.\n\n"
                "Each of those is reasonable and each of them changes the "
                "molecule. The list of changes is shown to you because the "
                "commonest silent error at this stage is docking something "
                "subtly different from what you meant - and a docking score "
                "for the wrong molecule looks exactly like a docking score "
                "for the right one."
            ),
            simple=(
                "You draw a molecule on paper as a flat picture. But a "
                "real molecule is a 3D object. Docking needs the proper "
                "3D shape, with all hydrogens in place. This step builds "
                "that.\n\n"
                "You give the molecule as a SMILES - a short text code "
                "for a structure. You can copy SMILES from PubChem or "
                "ChEMBL. For example, aspirin is:\n"
                "CC(=O)Oc1ccccc1C(=O)O\n\n"
                "While building the 3D shape, the software makes two "
                "changes you should know about:\n\n"
                "1. SALT REMOVED. Many drugs are sold as salts - "
                "amlodipine besylate, metoprolol tartrate. The salt part "
                "is separated here so the main molecule can be "
                "prepared. Review that choice.\n\n"
                "2. CHARGE STANDARDISED. Charges are neutralised where "
                "possible. The form present at a given pH can differ.\n\n"
                "Both are workflow choices. Every change is listed "
                "for you to check - because a docking score for the "
                "wrong molecule looks exactly like a score for the right "
                "one."
            ),
            not_mean=[
                "The 3D shape built here is one plausible conformation, not "
                "the shape the molecule adopts when bound. Docking will "
                "search other shapes; the starting one matters less than "
                "people expect.",
                "Neutralising the charge is a convention, not a claim about "
                "what exists at pH 7.4. A carboxylic acid is ionised in the "
                "body, and if the charge is important for binding at your "
                "site, you need to think about it rather than accept the "
                "default.",
                "A molecule that prepares cleanly is not thereby a good "
                "molecule. This step says nothing at all about activity.",
            ],
            learn_more=[
                ("What a SMILES string is",
                 "https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html"),
                ("RDKit molecule standardisation",
                 "https://www.rdkit.org/docs/RDKit_Book.html"),
            ],
        )

    def run(self, molecules: list[dict[str, str]], **kwargs: Any) -> StageResult:
        """
        `molecules` is a list of {"smiles": ..., "name": ...}.
        """
        if not molecules:
            return StageResult(
                success=False,
                message="No molecules were given. Add at least one.",
            )

        folder = self.project.step_folder(self.stage_id)
        results: list[Prepared] = []

        for item in molecules:
            prep = prepare_one(
                item.get("smiles", ""),
                item.get("name", ""),
                seed=self.SEED,
            )
            if prep.ok:
                safe = "".join(
                    c if c.isalnum() or c in "-_" else "_" for c in prep.name
                ) or "molecule"
                path = folder / f"{safe}.sdf"
                writer = Chem.SDWriter(str(path))
                mol3d = getattr(prep, "_mol3d")
                mol3d.SetProp("_Name", prep.name)
                writer.write(mol3d)
                writer.close()
                prep.sdf_path = path
            results.append(prep)

        good = [r for r in results if r.ok]
        bad = [r for r in results if not r.ok]

        self.record(
            ProvenanceEntry(
                stage=self.stage_id,
                tool="RDKit",
                version=self._rdkit_version(),
                parameters={
                    "conformer_method": "ETKDGv3",
                    "minimisation": "MMFF94, 500 steps",
                    "desalted": True,
                    "neutralised": True,
                    "tautomer": "input retained (canonicalisation off)",
                },
                seed=self.SEED,
                citation="RDKit: Open-source cheminformatics. https://www.rdkit.org",
            )
        )
        for entry in self.provenance:
            self.project.add_provenance(entry)

        control = self._change_control(results)

        return StageResult(
            success=bool(good),
            summary=(
                f"Prepared {len(good)} molecule(s)"
                + (f"; {len(bad)} could not be read" if bad else "")
            ),
            values={
                "prepared": [
                    {
                        "name": r.name,
                        "smiles": r.smiles,
                        "input_smiles": r.input_smiles,
                        "changes": r.changes,
                        "properties": r.properties,
                        "file": str(r.sdf_path) if r.sdf_path else "",
                    }
                    for r in good
                ],
                "failed": [
                    {"name": r.name, "input_smiles": r.input_smiles,
                     "error": r.error}
                    for r in bad
                ],
            },
            artefacts=[r.sdf_path for r in good if r.sdf_path],
            control=control,
            message="" if good else "None of the molecules could be read.",
        )

    # -- the control ------------------------------------------------------

    @staticmethod
    def _change_control(results: list[Prepared]) -> ControlOutcome:
        """
        Specification Stage 5: every modification listed, and the user
        must confirm the molecule is still the one they intended.

        Changed molecules produce a WARN that blocks progress -- not
        because anything went wrong, but because the user has to look.
        """
        changed = [r for r in results if r.ok and r.changed]
        failed = [r for r in results if not r.ok]

        if failed and not any(r.ok for r in results):
            return ControlOutcome(
                name="Is this still your molecule?",
                verdict=Verdict.FAIL,
                plain_verdict=(
                    "No molecule could be read from what you typed."
                ),
                what_to_do=(
                    "Check the SMILES text. The commonest mistakes are a "
                    "missing bracket, a capital letter where a small one "
                    "belongs (Cl is chlorine, CL is wrong), or an extra "
                    "space pasted from Word. The safest route is to search "
                    "your drug on PubChem and copy the 'Canonical SMILES' "
                    "line from there."
                ),
                detail={"failed": [r.name for r in failed]},
                blocks_progress=True,
            )

        if changed:
            lines = {r.name: r.changes for r in changed}
            names = ", ".join(lines)
            return ControlOutcome(
                name="Is this still your molecule?",
                verdict=Verdict.WARN,
                plain_verdict=(
                    f"{len(changed)} molecule(s) were altered during "
                    f"preparation ({names}). This is almost always "
                    "correct and expected."
                ),
                what_to_do=(
                    "Read the list of changes below once. Removing a salt "
                    "(besylate, tartrate, hydrochloride) or neutralising a "
                    "charge is fine - press Next. If instead you see that a "
                    "whole ring or group has gone, something was wrong with "
                    "the SMILES you typed; go back and copy it again from "
                    "PubChem."
                ),
                detail={"changes": lines},
                blocks_progress=True,
            )

        return ControlOutcome(
            name="Is this still your molecule?",
            verdict=Verdict.PASS,
            plain_verdict=(
                "Nothing had to be changed - the prepared molecules are "
                "exactly what you entered."
            ),
            what_to_do="Nothing to check. Press Next and dock them.",
            detail={},
        )

    # -- helpers ----------------------------------------------------------

    @staticmethod
    def _rdkit_version() -> str:
        try:
            from rdkit import rdBase
            return rdBase.rdkitVersion
        except Exception:  # noqa: BLE001
            return "unknown"

    def draw_comparison(self, prep: Prepared, path: Path) -> Path | None:
        """
        Side-by-side picture of what the user entered and what will
        actually be docked. Specification Stage 5, control.
        """
        a = Chem.MolFromSmiles(prep.input_smiles)
        b = Chem.MolFromSmiles(prep.smiles)
        if a is None or b is None:
            return None
        img = Draw.MolsToGridImage(
            [a, b],
            legends=["What you entered", "What will be docked"],
            molsPerRow=2,
            subImgSize=(320, 260),
        )
        data = getattr(img, "data", None)
        if data is None:                    # PIL image rather than SVG
            img.save(str(path))
        else:
            path.write_bytes(data if isinstance(data, bytes) else data.encode())
        return path

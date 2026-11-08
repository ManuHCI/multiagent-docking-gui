"""
Stage 8 -- Pharmacokinetics and toxicity.

The stage this audience understands best, and therefore the one where
over-reading the output is easiest and most costly.

What it does, in order of how much the numbers can be trusted:

  MEASURED PROPERTIES  size, lipophilicity, polar surface area, counts
                       of donors and acceptors. These are calculated
                       from the structure and are not predictions at
                       all. They are as reliable as the structure.

  RULES               Lipinski, Veber. Descriptions of what orally
                      absorbed drugs have tended to look like. Not
                      standards, and a great many useful drugs fail
                      them.

  PREDICTIONS         absorption, brain penetration, CYP inhibition,
                      toxicity. Models fitted to measured data. Useful
                      for triage, weak for anything finer.

The control is the applicability domain: where this molecule sits
relative to approved drugs on each property. A prediction about a
molecule inside that space is an interpolation. One about a molecule
outside it is an extrapolation dressed in the same decimal places, and
nothing in the output distinguishes the two unless something says so.

Runs offline on RDKit alone. ADMETlab 3.0 and ProTox add predictions
when there is a connection; SwissADME and pkCSM are stations the user
visits, because their terms do not permit automated access.
"""

from __future__ import annotations

import json
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

REFERENCE = Path(__file__).parent.parent / "assets" / "reference_drugs.json"
ADMETLAB = "https://admetlab3.scbdd.com/api/predict"

# What the properties mean, in words a clinician already uses.
PROPERTY_NOTES = {
    "mw": ("Molecular weight",
           "How big the molecule is. Above about 500 it gets harder to "
           "absorb from the gut, though plenty of useful drugs are "
           "bigger."),
    "logp": ("Lipophilicity (logP)",
             "How greasy it is. Too low and it will not cross a "
             "membrane; too high and it sticks to everything, including "
             "places you did not intend."),
    "tpsa": ("Polar surface area",
             "A measure of how polar the surface is. Above about 140 "
             "square angstroms gut absorption usually falls off; above "
             "about 90 it rarely enters the brain."),
    "hbd": ("Hydrogen bond donors",
            "Each one makes the molecule harder to pull through a "
            "membrane."),
    "hba": ("Hydrogen bond acceptors", "As above."),
    "rot": ("Rotatable bonds",
            "How floppy it is. More than about 10 and oral "
            "bioavailability tends to drop."),
}


def load_reference() -> list[dict[str, Any]]:
    try:
        return json.loads(REFERENCE.read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001
        return []


def describe(smiles: str) -> dict[str, Any] | None:
    """Properties calculated from the structure. Not predictions."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Crippen, Descriptors

    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "mw": round(Descriptors.MolWt(mol), 2),
        "logp": round(Crippen.MolLogP(mol), 2),
        "tpsa": round(Descriptors.TPSA(mol), 2),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "rot": Descriptors.NumRotatableBonds(mol),
        "heavy": mol.GetNumHeavyAtoms(),
        "fsp3": round(Descriptors.FractionCSP3(mol), 3),
        "rings": Descriptors.RingCount(mol),
    }


def rules(props: dict[str, Any]) -> dict[str, Any]:
    """
    Lipinski and Veber, reported as what they are: descriptions.
    """
    lip = []
    if props["mw"] > 500:
        lip.append("molecular weight over 500")
    if props["logp"] > 5:
        lip.append("logP over 5")
    if props["hbd"] > 5:
        lip.append("more than 5 hydrogen bond donors")
    if props["hba"] > 10:
        lip.append("more than 10 hydrogen bond acceptors")

    veber = []
    if props["rot"] > 10:
        veber.append("more than 10 rotatable bonds")
    if props["tpsa"] > 140:
        veber.append("polar surface area over 140")

    return {
        "lipinski_violations": len(lip),
        "lipinski_detail": lip,
        "veber_violations": len(veber),
        "veber_detail": veber,
    }


def percentile(value: float, values: list[float]) -> float:
    if not values:
        return float("nan")
    below = sum(1 for v in values if v < value)
    return round(100.0 * below / len(values), 1)


def applicability(props: dict[str, Any],
                  reference: list[dict[str, Any]]) -> dict[str, Any]:
    """
    Where does this molecule sit among approved drugs?

    Not a verdict on the molecule. A statement about how much the
    predictions below can be trusted: a model fitted to drug-like
    molecules is interpolating inside this space and extrapolating
    outside it, and the output looks identical either way.
    """
    out: dict[str, Any] = {}
    outside: list[str] = []
    for key in ("mw", "logp", "tpsa", "hbd", "hba", "rot"):
        values = [r[key] for r in reference if key in r]
        if not values:
            continue
        lo, hi = min(values), max(values)
        value = props[key]
        out[key] = {
            "value": value,
            "percentile_among_approved_drugs": percentile(value, values),
            "approved_range": [lo, hi],
            "outside": value < lo or value > hi,
        }
        if value < lo or value > hi:
            label = PROPERTY_NOTES.get(key, (key, ""))[0]
            out[key]["how_far"] = (
                f"{label} is {value}, outside the {lo}-{hi} seen in the "
                "approved drugs used for comparison"
            )
            outside.append(label.lower())
    out["_outside"] = outside
    return out


# ----------------------------------------------------------------------

class AdmetStage(Stage):
    stage_id = "08_admet"
    title = "Pharmacokinetics"
    tier = Tier.BUNDLED
    requires = ["05_ligands"]

    def check(self) -> CheckResult:
        if not self.project.is_complete("05_ligands"):
            return CheckResult(
                readiness=Readiness.NEEDS_INPUT,
                reason="No prepared molecules yet.",
                fix="Complete the molecule preparation step first.",
            )
        return CheckResult(readiness=Readiness.READY)

    def explain(self) -> Teaching:
        return Teaching(
            rule="Calculated descriptors help frame ADME questions; reference ranges do not validate prediction models.",
            rule_detail=[
                "This offline step calculates descriptors and compares them with a small reference-drug set.",
                "Optional external predictions require their own validation and applicability information.",
                "Neither a favourable descriptor profile nor a docking score establishes clinical benefit or safety.",
            ],
            button="Calculate descriptors and compare reference ranges.",
            simple=(
                "Target activity is only part of a drug's behaviour. Absorption, distribution, metabolism and excretion affect exposure.\n\n"
                "This step calculates molecular weight, calculated logP, polar surface area and related descriptors from the input structure. These are not measurements of behaviour in the body.\n\n"
                "Calculated logP estimates lipophilicity using an empirical method. Optional external services can estimate specific endpoints such as absorption or enzyme interactions; read each model's validation and applicability information.\n\n"
                "The check compares descriptors with a small set of approved drugs. It flags unusual properties and does not validate an absorption or toxicity model."
            ),
            panel=(
                "Molecular weight and atom counts follow the entered structure. Hydrogen-bond counts and polar surface area depend on descriptor definitions. Calculated logP is an estimate rather than a measured partition coefficient. Ionisation and tautomer choices can change the inputs and results.\n\n"
                "Lipinski's rules describe tendencies useful for screening. Passing does not establish oral bioavailability, and a violation is not automatically disqualifying. Interpret the rules in the context of the molecule and intended use.\n\n"
                "External prediction models cover particular endpoints and have their own training data, validation, uncertainty and limits. Check which endpoints a service actually reports. A comparison with reference-drug descriptor ranges cannot establish a prediction model's applicability."
            ),
            not_mean=[
                "Prediction reliability depends on the specific model and input. The descriptor-range check does not establish that reliability.",
                "A prediction of gut absorption is not oral bioavailability; exposure also depends on metabolism and other processes.",
                "Passing Lipinski's rules does not prove oral bioavailability or efficacy.",
                "A favourable toxicity prediction covers only its modelled endpoints and does not establish safety.",
                "Favourable descriptors do not establish target engagement or clinical benefit.",
            ],
            learn_more=[
                ("SwissADME - the original paper", "https://www.nature.com/articles/srep42717"),
                ("Lipinski's rule of five, as originally published", "https://pubmed.ncbi.nlm.nih.gov/11259830/"),
            ],
        )

    # -- run ---------------------------------------------------------------

    def run(self, use_online: bool = True, **kwargs: Any) -> StageResult:
        step = self.project.get_step("05_ligands")
        prepared = json.loads(
            step["values_json"] or "{}").get("prepared", []) if step else []
        if not prepared:
            return StageResult(success=False,
                               message="No prepared molecules were found.")

        reference = load_reference()
        molecules: list[dict[str, Any]] = []
        for item in prepared:
            props = describe(item["smiles"])
            if props is None:
                continue
            entry = {
                "name": item.get("name") or "molecule",
                "smiles": item["smiles"],
                "properties": props,
                "rules": rules(props),
                "applicability": applicability(props, reference),
            }
            molecules.append(entry)

        if not molecules:
            return StageResult(
                success=False,
                message="None of the molecules could be read.",
            )

        self.record(
            ProvenanceEntry(
                stage=self.stage_id,
                tool="RDKit descriptors",
                version=self._rdkit_version(),
                parameters={"molecules": len(molecules),
                            "reference_drugs": len(reference)},
                citation="RDKit: Open-source cheminformatics. "
                         "https://www.rdkit.org",
            )
        )
        for e in self.provenance:
            self.project.add_provenance(e)

        return StageResult(
            success=True,
            summary=(
                f"Profiled {len(molecules)} molecule(s). "
                + ", ".join(
                    f"{m['name']}: {m['rules']['lipinski_violations']} "
                    "Lipinski violation(s)" for m in molecules[:3])
            ),
            values={"molecules": molecules,
                    "reference_drugs": len(reference)},
            control=self._domain_control(molecules, len(reference)),
        )

    # -- the control -------------------------------------------------------

    @staticmethod
    def _domain_control(molecules: list[dict[str, Any]],
                        n_reference: int) -> ControlOutcome:
        """
        Specification Stage 8: say when a prediction is an
        extrapolation rather than an interpolation.

        Every molecule is reported, not just one. An earlier version
        picked a single compound to describe and silently said nothing
        about the rest, which is the wrong shape for a check: a user
        with five molecules needs to know which of them the numbers
        apply to.
        """
        detail: dict[str, Any] = {
            "approved drugs compared against": n_reference,
        }
        offenders: list[str] = []
        for mol in molecules:
            app = mol["applicability"]
            out = app.get("_outside", [])
            if out:
                offenders.append(f"{mol['name']} ({', '.join(out)})")
            parts = []
            for key, info in app.items():
                if key.startswith("_"):
                    continue
                label = PROPERTY_NOTES.get(key, (key, ""))[0]
                parts.append(
                    f"{label} {info['value']} "
                    f"({info['percentile_among_approved_drugs']:.0f}th pct"
                    + (", OUTSIDE range)" if info["outside"] else ")")
                )
            detail[mol["name"]] = "; ".join(parts)

        if not n_reference:
            return ControlOutcome(name="Descriptor reference ranges", verdict=Verdict.WARN,
                plain_verdict="The reference-drug set is unavailable. No range comparison was performed.",
                what_to_do="Use the calculated descriptors, and restore the reference file before comparing ranges.", detail=detail)
        if not offenders:
            return ControlOutcome(
                name="Are these descriptors within the reference-drug ranges?",
                verdict=Verdict.PASS,
                plain_verdict=(
                    f"All {len(molecules)} molecule(s) sit inside the "
                    "range of approved drugs on every property compared. "
                    "This is a descriptor-range comparison only. It "
                    "does not establish the applicability or accuracy "
                    "of an absorption, toxicity or other prediction model."
                ),
                what_to_do=(
                    "Report these as calculated descriptors, with the "
                    "method used. Treat calculated logP as an estimate. "
                    "Evaluate any external prediction using that model's "
                    "own validation and applicability information."
                ),
                detail=detail,
            )

        return ControlOutcome(
            name="Are these descriptors within the reference-drug ranges?",
            verdict=Verdict.WARN,
            plain_verdict=(
                f"{'; '.join(offenders)} - outside the range of approved "
                "drugs on those properties. This small reference set "
                "is not a prediction model's training set. Being outside "
                "its ranges flags a question to investigate; it does "
                "not by itself establish model validity or invalidity."
            ),
            what_to_do=(
                "Inspect the unusual properties in light of the intended "
                "route and target. For external predictions, read the "
                "model's validation and domain information. Do not infer "
                "oral bioavailability or safety from this range check."
            ),
            detail=detail,
            blocks_progress=True,
        )

    @staticmethod
    def _rdkit_version() -> str:
        try:
            from rdkit import rdBase
            return rdBase.rdkitVersion
        except Exception:  # noqa: BLE001
            return "unknown"

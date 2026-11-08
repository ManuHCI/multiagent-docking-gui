"""
SwissADME station.

SwissADME has no programming interface and its terms of use prohibit
automated access, so this software never calls it. It prepares the
input, hands the user to the site, and reads back the CSV the user
downloads.

Reference: Daina A, Michielin O, Zoete V. SwissADME: a free web tool to
evaluate pharmacokinetics, drug-likeness and medicinal chemistry
friendliness of small molecules. Sci Rep. 2017;7:42717.
"""

from __future__ import annotations

import csv
from pathlib import Path

from core.contract import (
    CheckResult,
    ControlOutcome,
    ProvenanceEntry,
    Readiness,
    StageResult,
    Teaching,
    Verdict,
)
from core.station import StationSpec, StationStage

# Columns worth putting in front of a clinician, with the plain reading.
HEADLINE = {
    "MW": "Molecular weight",
    "TPSA": "Polar surface area",
    "iLOGP": "Lipophilicity (iLOGP)",
    "GI absorption": "Gut absorption",
    "BBB permeant": "Crosses blood-brain barrier",
    "Lipinski #violations": "Lipinski violations",
    "Bioavailability Score": "Bioavailability score",
}


class SwissADMEStation(StationStage):
    stage_id = "08_admet"
    title = "ADME profile (SwissADME)"
    requires = ["05_ligands"]

    # -- the service ------------------------------------------------------

    def spec(self) -> StationSpec:
        return StationSpec(
            key="swissadme",
            name="SwissADME",
            url="https://www.swissadme.ch/index.php",
            method="post",
            field_name="smiles",
            why=(
                "It estimates how your molecule is likely to behave in the "
                "body - whether it is likely to be absorbed from the gut, "
                "whether it is likely to reach the brain, how lipophilic it "
                "is, and whether it looks like an orally active drug."
            ),
            cannot_tell_you=(
                "These are predictions from the properties of the molecule, "
                "not measurements. They say nothing about whether the "
                "molecule is active at your target, and nothing about dose. "
                "A molecule can pass every one of these filters and still "
                "be useless, and several marketed drugs fail them."
            ),
            instructions=[
                "Your molecules are already in the box. Press Run.",
                "Wait a few seconds - a panel appears for each molecule.",
                "Scroll to the top of the results and click CSV.",
                "Save the file somewhere you can find it, such as Downloads.",
            ],
            screenshot="swissadme_run.png",
            expect_back="You will have a file called swissadme.csv.",
            citation="Daina A, Michielin O, Zoete V. Sci Rep. 2017;7:42717.",
            terms_note=(
                "SwissADME is free for private study and internal research. "
                "This software does not contact it automatically - you do, "
                "by pressing the button."
            ),
        )

    def payload(self) -> str:
        """One SMILES per line, optionally followed by a name."""
        ligands = self.project.get_step("05_ligands")
        if not ligands:
            return ""
        import json

        values = json.loads(ligands["values_json"] or "{}")
        lines = []
        for item in values.get("prepared", []):
            smiles = item.get("smiles", "").strip()
            name = item.get("name", "").strip()
            if smiles:
                lines.append(f"{smiles} {name}".strip())
        return "\n".join(lines)

    # -- contract ---------------------------------------------------------

    def check(self) -> CheckResult:
        if not self.project.is_complete("05_ligands"):
            return CheckResult(
                readiness=Readiness.NEEDS_INPUT,
                reason="No prepared molecules yet.",
                fix="Complete the ligand preparation step first.",
            )
        if not self.payload():
            return CheckResult(
                readiness=Readiness.BLOCKED,
                reason="The ligand step finished but produced no molecules.",
                fix="Go back to ligand preparation and add at least one molecule.",
            )
        return CheckResult(readiness=Readiness.READY)

    def explain(self) -> Teaching:
        return Teaching(
            button="Send your molecules to SwissADME and bring the results back.",
            panel=(
                "Before a molecule can act on a target it has to get there and "
                "stay there long enough to matter. That is what ADME describes: "
                "absorption, distribution, metabolism and excretion. You already "
                "think in these terms clinically - why an aminoglycoside cannot "
                "be given by mouth, why a drug that crosses into the brain causes "
                "sedation when it also acts on relevant brain targets.\n\n"
                "SwissADME estimates these behaviours from the structure alone. "
                "It calculates size, lipophilicity and polar surface area, then "
                "uses models built from known drugs to predict gut absorption, "
                "brain penetration, whether the molecule is handled by P-"
                "glycoprotein, and which cytochrome P450 enzymes it is likely to "
                "inhibit, which can contribute to drug interactions.\n\n"
                "Treat the output as a triage tool. It is good at telling you a "
                "molecule has an obvious problem. It is much weaker at telling "
                "you a molecule is good."
            ),
            not_mean=[
                "A predicted value is an interpolation from the molecules the "
                "model was trained on. An unusual molecule may fall outside "
                "the model's reliable domain; review its applicability.",
                "'High GI absorption' is not oral bioavailability. It ignores "
                "first-pass metabolism, which is often the thing that decides.",
                "Passing Lipinski's rule of five is not a mark of quality. It "
                "was a description of orally absorbed drugs, not a standard to "
                "meet, and antibiotics and many other useful classes fail it.",
                "These predictions say nothing about whether the molecule binds "
                "your target. A perfect ADME profile on an inactive compound is "
                "a perfect profile of nothing.",
            ],
            learn_more=[
                ("SwissADME - the original paper",
                 "https://www.nature.com/articles/srep42717"),
                ("SwissADME help page", "https://www.swissadme.ch/help.php"),
            ],
        )

    # -- bringing the result back ----------------------------------------

    def ingest(self, path: Path) -> StageResult:
        path = Path(path)
        if not path.exists():
            return StageResult(
                success=False,
                message=f"Could not find {path.name}. Try dropping the file again.",
            )

        try:
            text = path.read_text(encoding="utf-8-sig")
            rows = list(csv.DictReader(text.splitlines()))
        except Exception as exc:  # noqa: BLE001
            return StageResult(
                success=False,
                message=(
                    "That file could not be read as a SwissADME CSV. "
                    "Make sure you clicked CSV on the results page rather "
                    "than saving the web page itself."
                ),
            )

        if not rows:
            return StageResult(
                success=False,
                message="That CSV is empty. Re-run SwissADME and download again.",
            )

        if "Canonical SMILES" not in rows[0] and "SMILES" not in rows[0]:
            return StageResult(
                success=False,
                message=(
                    "That looks like a CSV, but not one from SwissADME - "
                    "the expected columns are missing."
                ),
            )

        molecules = []
        for row in rows:
            entry = {"name": row.get("Molecule", "").strip()}
            for key, label in HEADLINE.items():
                if key in row:
                    entry[label] = row[key].strip()
            entry["_all"] = {k: v for k, v in row.items() if k}
            molecules.append(entry)

        expected = len([ln for ln in self.payload().splitlines() if ln.strip()])
        control = self._completeness(len(molecules), expected)

        self.record(
            ProvenanceEntry(
                stage=self.stage_id,
                tool="SwissADME",
                version="web service, accessed by the user",
                parameters={"molecules_submitted": expected},
                citation=self.spec().citation,
            )
        )
        for entry in self.provenance:
            self.project.add_provenance(entry)

        saved = self.project.step_folder(self.stage_id) / "swissadme.csv"
        saved.write_text(text, encoding="utf-8")

        return StageResult(
            success=True,
            summary=f"Read ADME predictions for {len(molecules)} molecule(s).",
            values={"molecules": molecules, "source": "SwissADME"},
            artefacts=[saved],
            control=control,
        )

    @staticmethod
    def _completeness(got: int, expected: int) -> ControlOutcome:
        if expected and got < expected:
            return ControlOutcome(
                name="All molecules accounted for",
                verdict=Verdict.WARN,
                plain_verdict=(
                    f"You sent {expected} molecules but the file has {got}. "
                    "Some may have failed on the website - check before you "
                    "rely on this."
                ),
                detail={"submitted": expected, "returned": got},
            )
        return ControlOutcome(
            name="All molecules accounted for",
            verdict=Verdict.PASS,
            plain_verdict=f"All {got} molecule(s) came back with predictions.",
            detail={"submitted": expected, "returned": got},
        )

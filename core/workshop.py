"""Workshop curriculum and participant files, independent of the interface."""
from __future__ import annotations

import json
import shutil
from copy import deepcopy
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

from .contract import ControlOutcome, ProvenanceEntry, StageResult, Verdict
from .project import Project
from .runtime import workshop_examples

EXAMPLES = workshop_examples()
MOLECULES = [
    {"name": "Erlotinib", "smiles": "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC"},
    {"name": "Gefitinib", "smiles": "COc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OCCCN1CCOCC1"},
    {"name": "Risperidone", "smiles": "Cc1nc2n(c(=O)c1CCN1CCC(CC1)c1noc3cc(F)ccc13)CCCC2"},
]

LESSONS = [
    {"key": "welcome", "title": "Ask a question", "time": "0–20 min", "stage": None,
     "why": "We will investigate possible ways molecules sit in a protein pocket, then judge the evidence.",
     "task": "Work in pairs. Predict which example molecule will rank highest in EGFR, and give one reason.",
     "look": "A prediction with a reason. You can change your conclusion when you see the results.",
     "question": "What would you need before saying a molecule works as a medicine?",
     "takeaway": "A computer model helps us choose questions for experiments. It does not measure clinical benefit."},
    {"key": "structure", "title": "Meet the protein and pocket", "time": "20–45 min", "stage": "02_structure",
     "why": "A protein has a 3D shape. A pocket is a region where a molecule may bind.",
     "task": "Rotate EGFR. Find the green reference molecule, then use Focus on pocket. Identify the protein and ligand to your partner.",
     "look": "Blue protein and green erlotinib from structure 1M17. This measured position gives us a test for our setup.",
     "question": "Why is a structure with a known bound molecule useful?",
     "takeaway": "The reference position helps locate the pocket and lets us test whether docking can recover a known pose."},
    {"key": "prepare", "title": "Prepare the protein", "time": "45–55 min", "stage": "03_receptor",
     "why": "A downloaded structure needs checking before a docking program can use it.",
     "task": "Compare the original and prepared views. Inspect what was kept, removed or repaired in the saved preparation record.",
     "look": "The reference ligand is kept for checking the setup; it will be excluded from the docking receptor.",
     "question": "Why should removing every non-protein component automatically make you cautious?",
     "takeaway": "Some metals, cofactors and waters matter. Preparation choices can change the question we are asking."},
    {"key": "molecules", "title": "Prepare three molecules", "time": "55–70 min", "stage": "05_ligands",
     "why": "SMILES is text describing a chemical structure. Preparation builds a possible 3D shape.",
     "task": "Select the example molecules and press Prepare selected molecules. Rotate each result and inspect any reported changes.",
     "look": "One 3D molecule per selection. The names help you read the results; the chemical structures are the actual inputs.",
     "question": "Is a prepared free molecule already in its bound position?",
     "takeaway": "No. A prepared shape is an input. Docking searches for possible positions and orientations in the pocket."},
    {"key": "break", "title": "Take a break", "time": "70–80 min", "stage": None,
     "why": "Pause and compare notes with your partner.", "task": "Explain target, ligand and pose without looking at the glossary.",
     "look": "Each partner can explain one term in their own words.", "question": "Which step is still unclear?",
     "takeaway": "Your notes are saved in your own workshop folder."},
    {"key": "dock", "title": "Choose the region and dock", "time": "80–115 min", "stage": "06_dock",
     "why": "The search box tells the computer where to try positions. Docking ranks the positions it finds using an approximate formula.",
     "task": "Inspect the amber search box. Predict the ranking, then run docking. Compare green reference erlotinib with its orange predicted pose.",
     "look": "The re-docking RMSD describes the positional difference. This app uses 2 Å as a teaching check, not a universal guarantee.",
     "question": "What would you do if the setup could not recover a known pose?",
     "takeaway": "Investigate the setup before interpreting a ranking. A favourable score alone does not establish binding or potency."},
    {"key": "interpret", "title": "Judge the evidence", "time": "115–145 min", "stage": "06_dock",
     "why": "A score needs context: the reference pose, comparison molecules and limitations of the method.",
     "task": "Read all checks. Explain the ranking to your partner, then inspect the saved SV2A case with a failed check.",
     "look": "A finished run can have failed checks. The comparison molecules are a small teaching baseline, not proof of inactivity or binding.",
     "question": "Does the most negative docking score establish the most potent drug? Why?",
     "takeaway": "Report what the model produced, which checks passed, and what remains uncertain."},
    {"key": "adme", "title": "Connect to pharmacology", "time": "145–165 min", "stage": "08_admet",
     "why": "Binding is only one part of a drug's behaviour. Size, polarity and lipophilicity affect the questions we ask about ADME.",
     "task": "Calculate properties. Compare molecular weight, calculated logP and polar surface area. Suggest one ADME question for an experiment.",
     "look": "These are calculated descriptors and rules. Online absorption or toxicity predictions are separate and optional.",
     "question": "Does satisfying Lipinski's rules prove oral bioavailability or safety?",
     "takeaway": "No. Rules describe tendencies. Calculated logP is an estimate, and a property-range check does not validate a prediction model."},
    {"key": "conclude", "title": "Write a defensible conclusion", "time": "165–180 min", "stage": None,
     "why": "A useful result includes evidence, limitations and a sensible next experiment.",
     "task": "Complete the worksheet and download it. Each pair shares one conclusion and one limitation. Open the optional MD preview if time allows.",
     "look": "Your prediction, observations, checks and proposed experiment are all recorded.",
     "question": "What is the strongest statement your results actually support?",
     "takeaway": "A pose is a hypothesis. Molecular dynamics explores motion under a model; staying in a pocket does not prove potency or safety."},
]

GLOSSARY = {
    "Target": "The protein we are investigating.",
    "Ligand": "The molecule we try placing in its pocket.",
    "Structure": "A 3D model of the protein or molecule.",
    "SMILES": "Text describing a chemical structure; you can copy it rather than memorise it.",
    "Search box": "The region the docking program explores.",
    "Pose": "One predicted position and orientation of a molecule.",
    "Re-docking": "Testing whether our setup can recover a known position.",
    "RMSD": "A measure of positional difference, here in angstroms (Å). Lower means closer positions.",
    "ADME": "Absorption, distribution, metabolism and excretion.",
}


def rebase(value, old: Path, new: Path):
    if isinstance(value, str):
        prefix = str(old.resolve())
        if value.lower().startswith(prefix.lower() + "\\") or value.startswith(prefix + "/"):
            return str(new.resolve()) + value[len(prefix):]
        return value
    if isinstance(value, dict):
        return {k: rebase(v, old, new) for k, v in value.items()}
    if isinstance(value, list):
        return [rebase(v, old, new) for v in value]
    return value


def copy_steps(source: Project, destination: Project, stage_ids: list[str]) -> None:
    """Copy actual results and files, with portable references; never invent scores."""
    for sid in stage_ids:
        row = source.get_step(sid)
        if not row or row["status"] != "complete":
            raise ValueError(f"The saved example needs a current {sid} result.")
        old_folder = source.folder / "steps" / sid
        if old_folder.exists():
            shutil.copytree(old_folder, destination.step_folder(sid), dirs_exist_ok=True)
        values = rebase(json.loads(row["values_json"] or "{}"), source.folder, destination.folder)
        c = json.loads(row["control_json"] or "null")
        if c:
            c["verdict"] = Verdict(c["verdict"])
        destination.save_step(sid, row["title"], StageResult(
            success=True, summary=row["summary"], values=values,
            control=ControlOutcome(**c) if c else None))
        destination.add_provenance(ProvenanceEntry(
            stage=sid, tool="Saved workshop example", version="1",
            parameters={"source_project": source.name, "source_updated": row["updated"]},
            citation="Copied saved calculation; not a new live run."))


def new_workshop(root: Path) -> Project:
    template = Project.open(EXAMPLES / "egfr")
    folder = root / ("Workshop EGFR " + datetime.now().strftime("%Y%m%d-%H%M%S") + "-" + uuid4().hex[:6])
    p = Project.create(folder, name=folder.name)
    try:
        copy_steps(template, p, ["01_target", "02_structure", "03_receptor", "04_site"])
        if (template.folder / "cache").exists():
            shutil.copytree(template.folder / "cache", p.folder / "cache", dirs_exist_ok=True)
        state = {"version": 1, "lesson": 0, "answers": {}, "done": [], "sources": {
            sid: "Saved preparation" for sid in ["01_target", "02_structure", "03_receptor", "04_site"]}}
        save_state(p, state)
        return p
    except Exception:
        p.close()
        raise
    finally:
        template.close()


def load_state(project: Project) -> dict:
    return json.loads((project.folder / "workshop.json").read_text(encoding="utf-8"))


def save_state(project: Project, state: dict) -> None:
    path = project.folder / "workshop.json"
    temp = path.with_suffix(".tmp")
    temp.write_text(json.dumps(state, indent=2, ensure_ascii=False), encoding="utf-8")
    temp.replace(path)


def use_saved_results(project: Project) -> None:
    template = Project.open(EXAMPLES / "egfr")
    try:
        copy_steps(template, project, ["01_target", "02_structure", "03_receptor", "04_site",
                                      "05_ligands", "06_dock", "08_admet"])
    finally:
        template.close()
    state = load_state(project)
    state["sources"].update({sid: "Saved example · not a live run"
                             for sid in ["01_target", "02_structure", "03_receptor", "04_site",
                                         "05_ligands", "06_dock", "08_admet"]})
    state["selected"] = [m["name"] for m in MOLECULES]
    state["done"] = [key for key in state["done"] if key in {"welcome", "break"}]
    save_state(project, state)


def evidence_source(project: Project, sid: str, state: dict | None = None) -> str:
    """Describe both where a result came from and whether it is current."""
    step = project.get_step(sid)
    if not step:
        return "No calculation recorded yet"
    if state is None:
        state = load_state(project)
    source = state.get("sources", {}).get(sid, "Source not recorded")
    if step["status"] == "stale":
        return source + " · Needs updating: inputs changed"
    if step["status"] != "complete":
        return source + " · Run failed: no current result"
    return source


def teaching_check(sid: str, outcome: dict, check_name: str = "", recorded: dict | None = None) -> dict:
    """Current teaching prose over an untouched recorded verdict and evidence.

    Historical files retain their original wording. This presentation copy
    avoids treating old explanatory claims as new scientific evidence.
    """
    display = deepcopy(outcome)
    if not outcome:
        return display
    if sid in {"01_target", "02_structure", "03_receptor", "04_site"}:
        status = {"pass": "passed its recorded criteria", "warn": "raised a warning",
                  "fail": "failed its recorded criteria"}.get(
                      outcome.get("verdict"), "did not provide an applicable conclusion")
        detail = outcome.get("detail", {})
        if sid == "01_target":
            text = (f"The target-evidence check {status}. Existing structures and drug associations can help choose a study target. They do not predict success for a new molecule or establish that the selected structure and pocket suit this question.")
            if "pdb_structures" in detail:
                text += f" The saved record lists {detail['pdb_structures']} structures; this is a recorded count, not a live database query."
            action = "Review the recorded evidence and biological question, then choose a suitable structure. Keep target-level evidence separate from evidence about your own molecules."
        elif sid == "02_structure":
            text = (f"The structure check {status}. Method, resolution, missing residues and local pocket quality need review together. Missing coordinates can have several causes; a global quality measure does not establish the accuracy of every part of a structure.")
            if "amino acids missing" in detail:
                text += f" Missing amino acids recorded: {detail['amino acids missing']}."
            action = "Inspect the reference ligand and pocket, including nearby missing residues and validation information. If the local evidence is insufficient, consider another structure or state the limitation."
        elif sid == "03_receptor":
            text = (f"The preparation check {status}. The retained components, removals and unresolved gaps are recorded below. These automated checks do not establish that every biologically necessary component was retained, or that protonation and other chemical choices are correct.")
            action = "Review retained and removed components, chains, repairs and gap distances in the context of this protein. Check relevant waters, metals, cofactors and protonation before docking."
        else:
            text = (f"The search-box geometry check {status}. Volume, enclosure and available-space estimates describe this chosen region under the method used. They do not prove a biological binding site, a snug molecular fit or druggability; a shallow surface region can still have biological binding partners.")
            action = "Review the box source and dimensions against the reference ligand or other biological evidence. Confirm that the intended search region is covered, and use re-docking when a suitable reference pose is available."
        display["plain_verdict"] = text
        display["what_to_do"] = action
        return display
    if sid != "06_dock":
        return display
    name = (check_name or outcome.get("name", "")).lower()
    verdict = outcome.get("verdict", "not_applicable")
    detail = display.get("detail", {})
    if "re-dock" in name or "known pose" in name:
        measured = detail.get("RMSD to the crystal pose")
        if measured is None and recorded and recorded.get("redock_rmsd") is not None:
            measured = str(recorded["redock_rmsd"]) + " Å"
        evidence = f" Recorded RMSD: {measured}." if measured is not None else ""
        meaning = {
            "pass": "The reference-pose check passed its recorded threshold. This supports recovery of this ligand's pose under this setup; it does not validate other poses or score rankings.",
            "warn": "The reference-pose check needs review. It does not provide a clear pass for recovery of this known pose.",
            "fail": "The reference-pose check failed. Review this setup before interpreting other docked poses or score rankings.",
        }.get(verdict, "Reference-pose reproduction was not established by this check.")
        action = "Report the RMSD, reference ligand, settings and check outcome. Review preparation, the search box and pose plausibility before drawing conclusions."
    elif "baseline" in name or "arbitrary molecule" in name:
        evidence = ""
        meaning = {
            "pass": "The baseline comparison passed the exercise's score-separation check. This small comparison set does not establish binding or rule out scoring bias.",
            "warn": "The baseline comparison raised a warning. It did not show clear separation under this exercise's check; a small score difference is not a reliable potency difference.",
            "fail": "The leading score did not outperform the recorded comparison baseline. This weakens a score-based claim and does not establish why the result occurred.",
        }.get(verdict, "No usable baseline conclusion was recorded. Review comparison availability and matching before interpreting the score.")
        action = "Report the scores, comparison set and recorded outcome. Review the poses and preparation; crossing a score threshold alone does not establish improved binding."
    elif "size" in name or "about fit" in name:
        evidence = ""
        meaning = {
            "pass": "This set passed the recorded size-association check. That does not establish that score ordering measures pocket fit or binding affinity.",
            "warn": "The size-association check raised a warning. Size may influence these scores; an association does not establish its cause or a potency ranking.",
            "fail": "The size-association check failed its exercise criterion. Interpret the ordering cautiously: correlation with size does not establish the cause of the ranking.",
        }.get(verdict, "The size-association check was not applicable to this recorded set.")
        action = "Report the association and score spread, inspect ligand efficiency and poses, and compare suitably similar compounds. These diagnostics do not validate potency rankings."
    else:
        evidence = ""
        meaning = "Review this recorded overall docking outcome alongside the individual checks and numerical evidence. A favourable score does not establish binding, potency or safety."
        action = "Report the outcome and its limitations. Investigate warnings or failed checks before interpreting rankings."
    display["plain_verdict"] = meaning + evidence
    display["what_to_do"] = action
    # Old aggregate rows embedded other checks as narrative detail. Keep all
    # numerical evidence, but apply the same teaching prose to those summaries.
    for key in list(detail):
        if key.lower().endswith(" check") and isinstance(detail[key], str):
            child_name = key[:-6]
            child = (recorded or {}).get("checks", {}).get(child_name)
            if child:
                detail[key] = teaching_check(sid, child, child_name)["plain_verdict"]
            else:
                detail[key] = "A historical narrative summary was recorded. See the individual check's current interpretation of the recorded numbers; its original wording remains in the project."
    return display


def docking_checks(recorded: dict, control: dict) -> dict:
    """Expose recorded checks and clearly label interpretations of legacy values."""
    checks = {name: teaching_check("06_dock", check, name, recorded)
              for name, check in recorded.get("checks", {}).items()}
    if not checks and control and "arbitrary molecule" in control.get("name", "").lower():
        baseline = teaching_check("06_dock", control, "baseline", recorded)
        baseline["detail"] = {key: value for key, value in baseline.get("detail", {}).items()
                              if not key.lower().endswith(" check")}
        checks["baseline"] = baseline
    if "re-docking" not in checks and "redock_rmsd" in recorded:
        from stages.dock import redock_verdict
        check = teaching_check("06_dock", asdict(redock_verdict(
            recorded["redock_rmsd"], recorded.get("redock_molecule", ""))), "re-docking", recorded)
        check["plain_verdict"] = ("Current teaching interpretation of recorded RMSD; no separate historical verdict was stored. "
                                  + check["plain_verdict"])
        checks["re-docking"] = check
    if "size" not in checks and recorded.get("size_confound"):
        from stages.dock import size_verdict
        check = teaching_check("06_dock", asdict(size_verdict(recorded["size_confound"])), "size", recorded)
        check["plain_verdict"] = ("Current teaching interpretation of recorded size data; no separate historical verdict was stored. "
                                  + check["plain_verdict"])
        if recorded["size_confound"].get("why"):
            check["plain_verdict"] += " " + recorded["size_confound"]["why"] + "."
        checks["size"] = check
    return checks


def worksheet(project: Project) -> str:
    state = load_state(project)
    out = ["EGFR workshop worksheet", project.name, "",
           "Activity ticks record practice; scientific checks are separate.", "",
           "Evidence sources:"]
    for sid in ["01_target", "02_structure", "03_receptor", "04_site",
                "05_ligands", "06_dock", "08_admet"]:
        out.append(sid + ": " + evidence_source(project, sid, state))
    out.append("")
    for lesson in LESSONS:
        answers = state["answers"].get(lesson["key"], {})
        completed = lesson["key"] in state.get("done", [])
        out.extend([lesson["title"] + " (" + lesson["time"] + ")",
                    "Activity: " + ("tried and discussed" if completed else "not ticked"),
                    lesson["question"]])
        for key, label in [("predict", "Prediction"), ("observe", "Observation"), ("explain", "Explanation")]:
            out.append(label + ": " + answers.get(key, "[not recorded]"))
        out.append("")

    def describe_check(label: str, check: dict) -> None:
        out.append(label + ": " + check.get("verdict", "not_applicable").upper())
        out.append(check.get("plain_verdict", "No scientific check was recorded."))
        if check.get("what_to_do"):
            out.append("Next action: " + check["what_to_do"])
        for name, value in check.get("detail", {}).items():
            if isinstance(value, (dict, list)):
                value = json.dumps(value, ensure_ascii=False)
            out.append("  " + name + ": " + str(value))

    out.append("CALCULATIONS AND SCIENTIFIC CHECKS")
    for sid in ["01_target", "02_structure", "03_receptor", "04_site", "05_ligands", "06_dock", "08_admet"]:
        step = project.get_step(sid)
        out.append(sid + ": " + project.display_status(sid)[2])
        if not step:
            out.extend(["No current result.", ""])
            continue
        if step["status"] != "complete":
            out.extend(["No current result. Rerun this step before interpreting it.",
                        "Earlier output files remain in the project for reference.", ""])
            continue
        out.append(step["summary"] or "Calculation finished")
        values = json.loads(step["values_json"] or "{}")
        original_control = json.loads(step["control_json"] or "{}")
        control = teaching_check(sid, original_control, recorded=values)
        if control:
            describe_check("Overall check", control)
        else:
            out.append("No scientific check was recorded.")
        if sid == "06_dock":
            rmsd = values.get("redock_rmsd")
            out.append("Re-docking RMSD: " + (str(rmsd) + " Å" if rmsd is not None else "not available"))
            out.append("Molecule | Docking score (kcal/mol)")
            for molecule in values.get("results", []):
                out.append(molecule["name"] + " | " + str(molecule.get("score", "Failed")))
            checks = docking_checks(values, original_control)
            for name, check in checks.items():
                describe_check(name.upper() + " check", check)
            size = values.get("size_confound", {})
            if "size" not in checks and size.get("usable") is False:
                out.append("SIZE check: NOT_APPLICABLE — " + size.get("why", "Too few usable molecules."))
            out.append("Re-docking tests recovery of this reference pose under this setup. A pass does not validate every molecule or establish binding, potency or safety.")
            out.append("Docking scores are model estimates, not measured potencies or a clinical ranking.")
        elif sid == "08_admet":
            out.append("Molecule | MW (g/mol) | Calculated logP | Polar surface (Å²) | Lipinski violations")
            for molecule in values.get("molecules", []):
                properties = molecule.get("properties", {})
                rules = molecule.get("rules", {})
                out.append(" | ".join(str(value) for value in [molecule["name"],
                    properties.get("mw", "—"), properties.get("logp", "—"),
                    properties.get("tpsa", "—"), rules.get("lipinski_violations", "—")]))
            out.append("These descriptors and rules do not establish oral bioavailability or safety. A reference-range check does not validate prediction models.")
        out.append("")
    out.append("Next experiment: " + state["answers"].get("next_experiment", "[not recorded]"))
    return "\n".join(out)

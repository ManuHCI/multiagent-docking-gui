"""Plain language study reports and portable, selected-pair research packages."""
from __future__ import annotations

import hashlib
import html
import io
import json
import math
from pathlib import Path
import zipfile

from . import viewer
from .runtime import resource_root

SWISSADME_URL = "https://www.swissadme.ch/index.php"


def preparation_review(project) -> list[dict]:
    r = viewer.values(project, "03_receptor")
    if not r:
        return []
    rows = [
        ("Hydrogen atoms added", r.get("hydrogens_added"), "These small atoms help describe chemical bonding and hydrogen bonds."),
        ("Missing atoms added", r.get("missing_atoms_added"), "Missing atoms in existing amino acids were modelled. Their positions are estimates."),
        ("Water molecules removed", r.get("waters_removed"), "These water molecules are absent from the prepared receptor. Important pocket waters may need a separate review."),
        ("Missing amino acids modelled", r.get("gaps_repaired"), "These amino acids were modelled rather than experimentally observed."),
        ("Protein segments left missing", len(r["gaps_left_open"]) if "gaps_left_open" in r else None, "These missing protein segments remain missing. Check whether they are near the drug pocket."),
    ]
    return [{"change": name, "value": "Not recorded" if value is None else value, "meaning": meaning}
            for name, value, meaning in rows]


def drug_review(drug: dict) -> str:
    changes = drug.get("changes", [])
    lines = list(changes) if changes else ["No salt or charge changes were recorded."]
    details = drug.get("preparation", {})
    if details:
        lines += [f"Hydrogen atoms added: {details['hydrogens_added']}. A 3D shape was built using {details['shape_method']}.",
                  "Shape relaxation using " + details["relaxation_method"] + (" finished within the set limit." if details["relaxation_converged"] else " reached its iteration limit. Review the starting geometry.")]
    else:
        lines.append("The saved file contains a prepared 3D structure. Detailed shape relaxation settings were not recorded by this earlier version.")
    lines += ["Docking can search other shapes.", "The prepared charge is a workflow choice. It may differ from the form present in the body."]
    return "\n".join(lines)


def _display(study, text) -> str:
    text = str(text)
    for drug in sorted(study.data["drugs"], key=lambda d: len(d["internal_name"]), reverse=True):
        text = text.replace(drug["internal_name"], drug["name"])
    return text


def report_text(study) -> str:
    from .study import results_rows, target_project
    from .workshop import docking_checks
    out = ["In Silico Docking Studio", "Docking report", study.data["name"], "Study status: " + study.data.get("status", "Not recorded"), "",
           "What this workflow does", "The software prepares protein and drug structures, then searches for ways each drug could fit in the selected protein pocket. Completed results are listed below.",
           "A docking score is a model estimate. A more negative number is favoured by this scoring method. It is not a measured binding strength or a recommended dose.", ""]
    for target in study.data["targets"]:
        out += ["Protein: " + target["name"], "PDB entry: " + (target.get("pdb_id") or "User supplied structure")]
        project = target_project(study, target["id"])
        try:
            out += ["", "Preparation changes"]
            for row in preparation_review(project):
                out.append(f"{row['change']}: {row['value']}. {row['meaning']}")
            receptor = viewer.values(project, "03_receptor")
            for label, key in (("Kept components", "kept"), ("Removed components", "removed")):
                if key in receptor:
                    out.append(label + ": " + (", ".join(str(x.get("label", "")) + " " + str(x.get("what", "")) for x in receptor[key]) or "None recorded"))
            for drug in viewer.values(project, "05_ligands").get("prepared", []):
                out += ["", "Prepared drug: " + _display(study, drug["name"]), _display(study, drug_review(drug)),
                        "Prepared SMILES: " + drug.get("smiles", "")]
            out += ["", "Docking results"]
            for row in results_rows(study):
                if row["target_id"] == target["id"]:
                    out.append(f"{row['drug']}: {row['status']}." + (f" Score: {row['score_kcal_mol']} kcal/mol." if row.get("score_kcal_mol") != "" else ""))
                    if row.get("error"):
                        out.append(row["error"])
            record = project.get_step("06_dock")
            dock = viewer.values(project, "06_dock")
            if dock:
                checks = docking_checks(dock, json.loads(record["control_json"] or "{}"))
                for name, check in checks.items():
                    label = {"baseline": "Score comparison", "re-docking": "Known pose recovery", "size": "Effect of molecular size"}.get(name, name)
                    state = {"pass": "Recorded criterion met", "fail": "Needs review", "warn": "Needs review", "not_applicable": "Not checked"}.get(check.get("verdict"), "Not checked")
                    out += ["", f"{label}: {state}", _display(study, check.get("plain_verdict", ""))]
                    if name == "baseline":
                        out.append("Comparison molecules are not proven nonbinders. A failed score comparison does not overturn experimental evidence for a known drug and target.")
                    if check.get("detail"):
                        out.append("Recorded numbers: " + _display(study, json.dumps(check["detail"], ensure_ascii=False)))
                out += ["", "Run settings: " + json.dumps(dock.get("settings", study.data.get("settings", {})), ensure_ascii=False),
                        "Search box: " + json.dumps(dock.get("box", {}), ensure_ascii=False)]
            else:
                out.append("No current completed docking result. Earlier results, if present, need updating.")
        finally:
            project.close()
        out.append("")
    out += ["Drug behaviour in the body", "The local drug descriptors are calculated chemical properties. SwissADME can provide additional ADME predictions.",
            "No SwissADME prediction was obtained automatically by this report. Paste the exported SMILES into " + SWISSADME_URL,
            "SwissADME predictions do not establish clinical pharmacokinetics, dose or safety.", "",
            "What to report", "Report the protein entry, prepared drug identity, search box, settings, scores, pose recovery and all comparison outcomes.",
            "Treat new docking results as hypotheses for further testing. A good score or a stable MD trajectory alone does not establish binding or activity."]
    return "\n".join(out) + "\n"


def report_html(study) -> bytes:
    return ("<!doctype html><html lang='en'><meta charset='utf-8'><title>Docking report</title>"
            "<style>body{max-width:900px;margin:40px auto;padding:20px;font:16px/1.65 system-ui;color:#21354b}pre{white-space:pre-wrap;overflow-wrap:anywhere}"
            "@media print{body{margin:0}}</style><pre>" + html.escape(report_text(study)) + "</pre></html>").encode("utf-8")


def swissadme_input(study) -> bytes:
    from .study import target_project
    prepared = {}
    # Use the prepared identities consistently; before preparation use explicit inputs.
    for target in study.data["targets"]:
        project = target_project(study, target["id"])
        try:
            for drug in viewer.values(project, "05_ligands").get("prepared", []):
                prepared.setdefault(drug["name"], drug["smiles"])
        finally:
            project.close()
    return ("\n".join(prepared.get(d["internal_name"], d["smiles"]) + " " + d["name"].replace("\n", " ").replace("\r", " ")
                       for d in study.data["drugs"]) + "\n").encode("utf-8")


def _file(project, filename: str) -> Path:
    path = Path(filename).resolve()
    if not filename or not path.is_relative_to(project.folder.resolve()) or not path.is_file():
        raise ValueError("A current structure file is missing or outside this project. Prepare and dock this pair again.")
    return path


def notebook(title: str, script: str, md: bool = False) -> bytes:
    def cell(kind, source):
        return {"cell_type": kind, "id": hashlib.sha256(source.encode("utf-8")).hexdigest()[:12], "metadata": {}, "source": source.splitlines(keepends=True),
                **({"execution_count": None, "outputs": []} if kind == "code" else {})}
    install = ("%pip -q install condacolab\nimport condacolab\ncondacolab.install()\n"
               if md else "%pip -q install vina==1.2.7\n")
    cells = [cell("markdown", "# " + title + "\nUpload this notebook at https://colab.research.google.com/. Run the cells in order.\n"
                  + ("For MD, choose a GPU in Runtime > Change runtime type when available. Installation restarts the runtime. After that restart, continue with the next cell.\n"
                     "This notebook is a starting protocol. Review job.json and README.txt. Membrane proteins such as SV2A need a prepared lipid membrane system.\n" if md else
                     "This uses the recorded receptor, drug input, box and settings. Scores can differ between platforms and processor counts. Comparison checks are recorded in the report and are not rerun by this selected-pair notebook.\n")),
             cell("code", install)]
    if md:
        cells.append(cell("code", "import subprocess\nsubprocess.run(['conda', 'install', '-y', '-c', 'conda-forge', 'openmm>=8.1,<9', 'openmmforcefields', 'openff-toolkit', 'ambertools', 'rdkit', 'mdtraj', 'matplotlib'], check=True)\n"))
    cells += [cell("markdown", "Upload the matching ZIP downloaded from Docking Studio. The ZIP contains the selected protein, selected drug, code and recorded settings.\n"),
              cell("code", "from google.colab import files\nfrom pathlib import Path\nimport zipfile, os\nuploaded = files.upload()\narchives = [name for name in uploaded if name.endswith('.zip')]\nassert len(archives) == 1, 'Upload exactly one Docking Studio ZIP.'\nwork = Path('/content/docking_studio_job').resolve()\nwork.mkdir(exist_ok=True)\nwith zipfile.ZipFile(archives[0]) as archive:\n    for member in archive.infolist():\n        assert (work / member.filename).resolve().is_relative_to(work), 'Unsafe archive path'\n    archive.extractall(work)\nos.chdir(work)\nprint(Path('README.txt').read_text())\n"),
              cell("code", "import subprocess, sys\nsubprocess.run([sys.executable, '" + script + "'], check=True)\n"),
              cell("code", "import shutil\nfrom google.colab import files\nshutil.make_archive('/content/studio_results', 'zip', 'output')\nfiles.download('/content/studio_results.zip')\n")]
    return json.dumps({"nbformat": 4, "nbformat_minor": 5, "metadata": {"kernelspec": {"name": "python3", "display_name": "Python 3"}, "colab": {"provenance": []}}, "cells": cells}, indent=2).encode("utf-8")


def _zip(files: dict[str, bytes]) -> bytes:
    stream = io.BytesIO()
    with zipfile.ZipFile(stream, "w", zipfile.ZIP_DEFLATED) as archive:
        for name, data in files.items():
            archive.writestr(name, data)
    return stream.getvalue()


def export_pair(study, target_id: str, drug_id: str, kind: str = "docking", *, system_type: str = "water", nanoseconds: float = 1.0) -> bytes:
    """Export exactly the selected current pair; never fall back to a different drug."""
    from .study import target_project, results_rows
    target = next((t for t in study.data["targets"] if t["id"] == target_id), None)
    drug = next((d for d in study.data["drugs"] if d["id"] == drug_id), None)
    if target is None or drug is None:
        raise ValueError("Select a protein and a drug from this study.")
    pair = next(r for r in results_rows(study) if r["target_id"] == target_id and r["drug_id"] == drug_id)
    if pair["status"] != "completed":
        raise ValueError("This pair needs a current completed docking run before export.")
    if kind not in {"docking", "md"} or system_type not in {"water", "membrane"}:
        raise ValueError("Choose docking or MD and a supported system type.")
    if kind == "md":
        try:
            nanoseconds = float(nanoseconds)
        except (TypeError, ValueError):
            raise ValueError("Enter an MD duration greater than zero and no more than 1000 ns.") from None
        if not math.isfinite(nanoseconds) or not 0 < nanoseconds <= 1000:
            raise ValueError("Enter an MD duration greater than zero and no more than 1000 ns.")
    project = target_project(study, target_id)
    try:
        dock = viewer.values(project, "06_dock")
        result = next(r for r in dock["results"] if r["name"] == drug["internal_name"] and "score" in r)
        pose = _file(project, result["file"])
        ligand_input = _file(project, result.get("input_file") or str(pose.with_name(pose.name.removesuffix("_docked.pdbqt") + ".pdbqt")))
        receptor = _file(project, dock["receptor"])
        from stages.dock import DEFAULT_SEED, DEFAULT_POSES
        settings = dock.get("settings")
        if not settings:
            # Historical settings must come from this target's own run provenance.
            runs = [p for p in project.provenance() if p.get("stage_id") == "06_dock" and p.get("tool") == "AutoDock Vina"]
            if not runs:
                raise ValueError("This historical run has no recorded search settings. Run docking again before exporting code.")
            run = runs[-1]
            parameters = run.get("parameters", {})
            if isinstance(parameters, str):
                parameters = json.loads(parameters)
            settings = {"exhaustiveness": parameters.get("exhaustiveness"), "seed": run.get("seed", DEFAULT_SEED), "n_poses": DEFAULT_POSES, "vina_version": run.get("version", "unknown")}
        if not settings.get("exhaustiveness"):
            raise ValueError("Search effort was not recorded. Run docking again before exporting code.")
        job = {"target": target["name"], "pdb_id": target.get("pdb_id", ""), "drug": drug["name"], "smiles": result["smiles"],
               "recorded_score_kcal_mol": result["score"], "box": dock["box"], "settings": settings,
               "checks": dock.get("checks", {}), "redock_rmsd_angstrom": dock.get("redock_rmsd"),
               "purpose": "Selected pair rerun; comparison checks are preserved, not rerun."}
        files = {"receptor.pdbqt": receptor.read_bytes(), "ligand.pdbqt": ligand_input.read_bytes(), "recorded_pose.pdbqt": pose.read_bytes(),
                 "Docking_report.txt": report_text(study).encode("utf-8")}
        if kind == "docking":
            files["run_docking.py"] = (resource_root() / "assets/templates/run_docking.py").read_bytes()
            files["Run_Docking.ipynb"] = notebook("Dock your selected protein and drug", "run_docking.py")
            files["README.txt"] = ("Selected protein: " + target["name"] + "\nSelected drug: " + drug["name"] + "\n\n"
                "Upload Run_Docking.ipynb to Google Colab. Run its cells and upload this complete ZIP when asked.\n"
                "The code uses the exact prepared PDBQT inputs, recorded box and run settings. No new protein or drug is downloaded.\n"
                "It runs AutoDock Vina through its Python interface. It does not repeat preparation or the comparison controls.\n"
                "The desktop uses the Vina command line interface. Differences in processors, platforms and interfaces can affect results.\n"
                "The original checks remain in job.json and Docking_report.txt. Interpret scores with those checks.\n"
                "Local use: install vina==1.2.7, then run python run_docking.py in this folder.\n").encode("utf-8")
        else:
            files.update(_md_files(project, pose, drug, target, system_type, nanoseconds))
            job.update(system_type=system_type, nanoseconds=nanoseconds, timestep_fs=2.0, temperature_K=310,
                       salt_M=0.15, padding_nm=1.0, seed=2026, protein_forcefield="AMBER ff14SB", ligand_forcefield="OpenFF 2.2.1 with AM1-BCC",
                       purpose="MD starting protocol, not an executed simulation")
            job["water_model"] = "TIP3P"
            if system_type == "membrane":
                job.update(protein_forcefield="Defined in the supplied membrane System; record builder details",
                           ligand_forcefield="Defined in the supplied membrane System; record builder details",
                           water_model="Defined in the supplied membrane System", salt_M=None, padding_nm=None)
        job["input_sha256"] = {name: hashlib.sha256(data).hexdigest() for name, data in files.items() if name.endswith((".pdbqt", ".pdb", ".sdf"))}
        files["job.json"] = json.dumps(job, indent=2, ensure_ascii=False).encode("utf-8")
        return _zip(files)
    finally:
        project.close()


def _md_files(project, pose, drug, target, system_type, nanoseconds):
    from meeko import PDBQTMolecule, RDKitMolCreate
    from rdkit import Chem
    molecules = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(pose.read_text(encoding="utf-8"), poses_to_read=1, skip_typing=True))
    if len(molecules) != 1 or molecules[0] is None or molecules[0].GetNumConformers() != 1:
        raise ValueError("The docked drug cannot be reconstructed with its bond orders. Run docking again with this version.")
    mol = Chem.AddHs(molecules[0], addCoords=True)
    mol.SetProp("_Name", drug["name"])
    ligand_sdf = (Chem.MolToMolBlock(mol) + "\n$$$$\n").encode("utf-8")
    receptor = viewer.values(project, "03_receptor")
    prepared = _file(project, receptor["prepared_file"]).read_text(encoding="utf-8")
    raw = viewer.local_text(project, viewer.values(project, "02_structure").get("chosen", {}).get("local_file", ""))
    import re
    if system_type == "water" and re.search(r"\bSV2[ABC]\b|MEMBRANE PROTEIN", target["name"] + "\n" + raw[:1500], re.IGNORECASE):
        raise ValueError("This target is a membrane protein. Choose membrane mode and prepare a lipid system before MD.")
    reference = receptor.get("ligand", "").split()
    def keep(line):
        if not line.startswith(("ATOM  ", "HETATM", "TER", "SSBOND", "LINK")):
            return False
        if len(reference) >= 2 and line.startswith("HETATM"):
            label = line[17:20].strip() + " " + line[21:22] + line[22:27].strip()
            return label != " ".join(reference[:2])
        return True
    protein = "\n".join(line for line in prepared.splitlines() if keep(line)) + "\nEND\n"
    instruction = ("Membrane mode: build a lipid membrane system for this selected complex using a suitable membrane builder, such as CHARMM-GUI (https://www.charmm-gui.org/).\n"
                   "Export an OpenMM System as membrane_system.xml and its matching topology and coordinates as membrane.pdb into this folder.\n"
                   "The system must contain the selected drug as residue LIG, with correct parameters, solvent, ions and a periodic box.\n"
                   "Review protein identity, drug chemistry, membrane orientation and equilibration. The runner checks particle counts and the selected drug's element counts, but these checks cannot establish correct system construction.\n"
                   "The script stops until those files are provided. It does not simulate SV2A in water alone.\n" if system_type == "membrane" else
                   "Water mode is for a soluble protein with a standard, parameterisable drug. It is not suitable for membrane proteins such as SV2A.\n"
                   "Review the retained components and missing protein regions first. Unknown cofactors need appropriate parameters; the script stops on unsupported chemistry.\n")
    return {"protein.pdb": protein.encode("utf-8"), "ligand.sdf": ligand_sdf,
            "run_md.py": (resource_root() / "assets/templates/run_md.py").read_bytes(),
            "Run_MD.ipynb": notebook("MD for your selected protein and drug", "run_md.py", md=True),
            "README.txt": (f"Selected protein: {target['name']}\nSelected drug: {drug['name']}\nRequested production duration: {nanoseconds:g} ns\n\n" + instruction +
                "\nThis is downloadable code. No MD simulation has been run by the desktop app.\n"
                "ligand.sdf contains the docked coordinates and chemical bond orders recovered from Meeko mapping.\n"
                "Upload Run_MD.ipynb to Colab and follow its cells. Installation restarts the runtime before you upload this ZIP.\n"
                "The default soluble protocol uses AMBER ff14SB, TIP3P water, OpenFF 2.2.1, AM1-BCC charges, 310 K and a 2 fs timestep.\n"
                "It minimises energy, runs 100 ps of NVT equilibration and 100 ps of NPT equilibration, then production. These durations require review for your system.\n"
                "Membrane mode preserves the supplied System's force fields and barostat. Its settings must be documented separately.\n"
                "Results and checkpoints are saved in output. Download them before a temporary Colab session ends.\n"
                "Resume requires the saved system, topology and checkpoint on a compatible OpenMM platform.\n"
                "A stable trajectory alone does not prove binding. Review trajectories, replicates, contacts and force field assumptions.\n"
                "References: https://docs.openmm.org/ ; https://github.com/openmm/openmmforcefields ; https://docs.openforcefield.org/\n").encode("utf-8")}

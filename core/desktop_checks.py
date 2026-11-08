"""Offline functional acceptance checks for a distributed desktop build.

All calculations and project edits happen in a temporary writable directory.
Bundled reference projects are copied before opening because opening a project
may update its portable paths. Numerical results here are runtime checks,
not a new scientific validation of the teaching example.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import tempfile
import time
import traceback

from .runtime import resource_root


def _small_peptide(source: Path, destination: Path) -> None:
    """Take three contiguous heavy-atom residues for a fast preparation check."""
    residues: list[str] = []
    lines: list[str] = []
    for line in source.read_text(encoding="utf-8").splitlines():
        if not line.startswith("ATOM  ") or line[21] != "A":
            continue
        if line[76:78].strip() in {"H", "D"}:
            continue
        key = line[21:27]
        if key not in residues:
            if len(residues) == 3:
                break
            residues.append(key)
        lines.append(line)
    if len(residues) != 3:
        raise ValueError("The bundled EGFR structure lacks the peptide fixture.")
    destination.write_text("\n".join(lines) + "\nTER\nEND\n", encoding="utf-8")


def _stored_values(case: Path, stage_id: str) -> dict:
    # Read-only: never relocate or modify the distributed reference database.
    with sqlite3.connect((case / "project.db").resolve().as_uri() + "?mode=ro",
                         uri=True) as connection:
        row = connection.execute("SELECT values_json FROM steps WHERE stage_id=?",
                                 (stage_id,)).fetchone()
    if row is None:
        raise ValueError(f"The reference case is missing {stage_id}.")
    return json.loads(row[0])


def run_scientific_checks(output_path: str) -> int:
    """Write a JSON report and return 0 only when every functional check passes."""
    started = time.monotonic()
    checks: list[dict] = []
    resources = resource_root()
    case = resources / "assets" / "workshop" / "egfr"

    def record(name, function):
        check_started = time.monotonic()
        try:
            detail = function()
            checks.append({"name": name, "passed": True, "detail": detail,
                           "seconds": round(time.monotonic() - check_started, 3)})
        except Exception as exc:  # Each missing resource remains a separate finding.
            checks.append({"name": name, "passed": False,
                           "error": f"{type(exc).__name__}: {exc}",
                           "traceback": traceback.format_exc(),
                           "seconds": round(time.monotonic() - check_started, 3)})

    temporary_parent = os.environ.get("INSILICO_DATA_DIR") or None
    if temporary_parent:
        Path(temporary_parent).mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="desktop_acceptance_", dir=temporary_parent) as temp:
        work = Path(temp)

        def chemistry():
            from rdkit import Chem, RDConfig
            from rdkit.Chem import ChemicalFeatures
            from rdkit.Chem.Draw import rdMolDraw2D
            from .project import Project
            from stages.admet import AdmetStage
            from stages.ligands import LigandStage

            feature_file = Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"
            factory = ChemicalFeatures.BuildFeatureFactory(str(feature_file))
            project = Project.create(work / "chemistry", name="Desktop acceptance")
            try:
                result = LigandStage(project).run(molecules=[{
                    "smiles": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]", "name": "Aspirin sodium"}])
                if not result.success or len(result.values["prepared"]) != 1:
                    raise RuntimeError(result.message or "3D molecule preparation failed.")
                prepared = result.values["prepared"][0]
                if not prepared["changes"]:
                    raise RuntimeError("The preparation change record is empty.")
                molecule = next(Chem.SDMolSupplier(prepared["file"], removeHs=False))
                if molecule is None or not molecule.GetConformer().Is3D():
                    raise RuntimeError("The prepared SDF lacks a 3D conformer.")
                features = factory.GetFeaturesForMol(molecule)
                # Match the UI's SVG route; this RDKit build has no Cairo renderer.
                drawer = rdMolDraw2D.MolDraw2DSVG(300, 220)
                drawer.DrawMolecule(molecule)
                drawer.FinishDrawing()
                image = work / "aspirin.svg"
                image.write_text(drawer.GetDrawingText(), encoding="utf-8")
                if image.stat().st_size < 500 or "<svg" not in image.read_text():
                    raise RuntimeError("The 2D molecule image is empty.")
                project.save_step("05_ligands", "Prepared molecule", result)
                adme = AdmetStage(project).run(use_online=False)
                if not adme.success or adme.values["reference_drugs"] < 1:
                    raise RuntimeError(adme.message or "ADME reference data are missing.")
                mw = adme.values["molecules"][0]["properties"]["mw"]
                if abs(mw - 180.16) > 0.1:
                    raise RuntimeError(f"Unexpected prepared aspirin molecular weight: {mw}")
                return {"molecular_weight": mw, "features": len(features),
                        "reference_drugs": adme.values["reference_drugs"],
                        "preparation_changes": prepared["changes"]}
            finally:
                project.close()

        record("3D chemistry, RDKit data, molecular drawing and ADME", chemistry)

        def protein():
            import openmm
            from openmm.app import ForceField, PDBFile
            from .project import Project
            from stages.receptor import ReceptorStage

            fixture = work / "peptide.pdb"
            _small_peptide(case / "steps" / "02_structure" / "1M17.pdb", fixture)
            project = Project.create(work / "protein", name="Protein acceptance")
            try:
                result = ReceptorStage(project).run(pdb_path=fixture, repair_gaps=False)
                if not result.success or result.values.get("hydrogens_added", 0) < 1:
                    raise RuntimeError(result.message or "Protein hydrogens were not added.")
                prepared = PDBFile(result.values["prepared_file"])
                forcefield = ForceField("amber14-all.xml")
                system = forcefield.createSystem(prepared.topology)
                integrator = openmm.VerletIntegrator(0.001)
                context = openmm.Context(system, integrator,
                                         openmm.Platform.getPlatformByName("Reference"))
                context.setPositions(prepared.positions)
                energy = context.getState(getEnergy=True).getPotentialEnergy()
                value = energy.value_in_unit(openmm.unit.kilojoules_per_mole)
                if not math.isfinite(value):
                    raise RuntimeError("The prepared peptide has a non-finite energy.")
                del context, integrator
                return {"hydrogens_added": result.values["hydrogens_added"],
                        "atoms": prepared.topology.getNumAtoms(),
                        "platforms": [openmm.Platform.getPlatform(i).getName()
                                      for i in range(openmm.Platform.getNumPlatforms())],
                        "energy_kj_mol": round(value, 4)}
            finally:
                project.close()

        record("Protein preparation and OpenMM force-field runtime", protein)

        def docking_input():
            from stages.dock import DockStage
            path, error = DockStage._ligand_pdbqt("CCO", work, "ethanol")
            if path is None or not path.is_file() or "ATOM" not in path.read_text():
                raise RuntimeError(error or "Meeko did not write docking input.")
            return {"pdbqt_bytes": path.stat().st_size}

        record("Meeko docking-input generation", docking_input)

        def vina_score():
            from stages.dock import find_vina
            executable = find_vina()
            if not executable:
                raise RuntimeError("The bundled Vina executable was not found.")
            docking = case / "steps" / "06_dock"
            source = (docking / "Erlotinib_docked.pdbqt").read_text(encoding="utf-8")
            first_pose = "\n".join(line for line in source.split("ENDMDL", 1)[0].splitlines()
                                   if not line.startswith("MODEL")) + "\n"
            pose = work / "erlotinib_first_pose.pdbqt"
            pose.write_text(first_pose, encoding="utf-8")
            box = _stored_values(case, "04_site")["box"]
            command = [executable, "--receptor", str(docking / "receptor.pdbqt"),
                       "--ligand", str(pose), "--score_only", "--cpu", "1"]
            for key in ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z"):
                command.extend(["--" + key, str(box[key])])
            result = subprocess.run(command, capture_output=True, text=True, timeout=45,
                                    creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            text = result.stdout + result.stderr
            match = re.search(r"Estimated Free Energy of Binding\s*:\s*([-+\d.]+)", text)
            if result.returncode != 0 or match is None:
                raise RuntimeError(text[-1800:] or "Vina produced no score.")
            score = float(match.group(1))
            if not math.isfinite(score):
                raise RuntimeError("Vina produced a non-finite score.")
            return {"score_kcal_mol": score, "mode": "score_only; no docking search repeated"}

        record("Bundled Vina real score calculation", vina_score)

        def workshop():
            from . import workshop as curriculum
            from .viewer import document, payload

            original_examples = curriculum.EXAMPLES
            copied_examples = work / "examples"
            shutil.copytree(resources / "assets" / "workshop", copied_examples)
            curriculum.EXAMPLES = copied_examples
            project = None
            try:
                project = curriculum.new_workshop(work / "projects")
                state = curriculum.load_state(project)
                state["answers"]["welcome"] = {"predict": "Acceptance-test note"}
                curriculum.save_state(project, state)
                curriculum.use_saved_results(project)
                if curriculum.load_state(project)["answers"]["welcome"]["predict"] != "Acceptance-test note":
                    raise RuntimeError("Loading saved results lost the student's notes.")
                data = payload(project, "06_dock", "Erlotinib")
                if not data["protein"] or not data["ligand"] or not data["overlay"]:
                    raise RuntimeError("The saved reference overlay is incomplete.")
                if "/assets/vendor/3Dmol-min.js" not in document(data):
                    raise RuntimeError("The molecular view lost its local 3D library.")
                worksheet = curriculum.worksheet(project)
                if "Acceptance-test note" not in worksheet or "Re-docking RMSD: 1.072" not in worksheet:
                    raise RuntimeError("The worksheet omitted notes or numerical evidence.")
                if "not a live run" not in worksheet:
                    raise RuntimeError("Saved results are not labelled in the worksheet.")
                return {"worksheet_characters": len(worksheet),
                        "protein_atoms": sum(line.startswith("ATOM") for line in data["protein"].splitlines()),
                        "saved_result_labels": True}
            finally:
                if project is not None:
                    project.close()
                curriculum.EXAMPLES = original_examples

        record("Offline workshop, saved results, notes, overlay and worksheet", workshop)

    report = {"passed": bool(checks) and all(check["passed"] for check in checks),
              "seconds": round(time.monotonic() - started, 3), "checks": checks}
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return 0 if report["passed"] else 1

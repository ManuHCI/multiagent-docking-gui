"""Batch studies retain pair identity, target geometry and failed results."""
from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch
import zipfile

from core.contract import ControlOutcome, StageResult, Verdict
from core.study import (DrugSpec, TargetSpec, create_example_study, create_study,
                        define_binding_sites, export_csv, export_package, list_studies,
                        open_study, results_rows, run_docking, run_preparation, target_project)


def good(values=None, verdict=Verdict.PASS, blocks=False):
    return StageResult(True, summary="Calculated", values=values or {},
                       control=ControlOutcome("Scientific check", verdict, "Recorded scientific outcome", blocks_progress=blocks))


class StudyTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="study_test_")
        self.root = Path(self.temporary.name)
        self.protein = self.root / "input.pdb"
        self.protein.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n", encoding="utf-8")
        self.study = create_study(self.root, "Two targets", [TargetSpec("Target A", file=str(self.protein)),
                                  TargetSpec("Target B", file=str(self.protein))],
                                 [DrugSpec("Drug A", "CCO"), DrugSpec("Drug B", "CCN")])

    def tearDown(self):
        self.temporary.cleanup()

    def seed(self):
        for index, target in enumerate(self.study.data["targets"]):
            project = target_project(self.study, target["id"])
            try:
                box = {"center_x": index * 20, "center_y": 0, "center_z": 0,
                       "size_x": 20, "size_y": 20, "size_z": 20}
                project.save_step("01_target", "Target", good())
                project.save_step("02_structure", "Structure", good())
                project.save_step("03_receptor", "Protein", good({"prepared_file": str(self.protein), "ligand": "LIG A2"}))
                project.save_step("04_site", "Site", good({"box": box}))
                project.save_step("05_ligands", "Drugs", good({"prepared": [
                    {"name": drug["internal_name"], "smiles": drug["smiles"]} for drug in self.study.data["drugs"]]}))
                target.update(status="prepared", box=box)
            finally:
                project.close()
        self.study.save()

    def mock_dock(self, failed_target=""):
        class Dock:
            stage_id = "06_dock"
            title = "Dock"
            def __init__(inner, project):
                inner.project = project
            def run(inner, **kwargs):
                if inner.project.name == failed_target:
                    return StageResult(False, message="This target could not be docked")
                molecules = json.loads(inner.project.get_step("05_ligands")["values_json"])["prepared"]
                box = json.loads(inner.project.get_step("04_site")["values_json"])["box"]
                results = []
                for number, molecule in enumerate(molecules):
                    pose = inner.project.step_folder("06_dock") / (molecule["name"] + ".pdbqt")
                    pose.write_text("MODEL 1\nENDMDL\n", encoding="utf-8")
                    results.append({"name": molecule["name"], "score": -5.0-number-box["center_x"]/20,
                                    "file": str(pose), "smiles": molecule["smiles"]})
                return good({"results": results, "box": box, "checks": {"baseline": {"verdict": "warn"}}}, Verdict.WARN)
        return Dock

    def test_matrix_has_every_target_drug_pair_before_execution(self):
        rows = results_rows(self.study)
        self.assertEqual(len(rows), 4)
        self.assertEqual(len({(r["target_id"], r["drug_id"]) for r in rows}), 4)
        self.assertTrue(all(row["status"] == "pending" for row in rows))
        self.assertEqual(list_studies(self.root)[0]["name"], "Two targets")

    def test_two_by_two_run_uses_each_targets_own_box_and_keeps_warning(self):
        self.seed()
        with patch("stages.dock.DockStage", self.mock_dock()):
            run_docking(self.study, exhaustiveness=1, n_decoys=1)
        rows = results_rows(self.study)
        self.assertEqual([row["score_kcal_mol"] for row in rows], [-5.0, -6.0, -6.0, -7.0])
        self.assertTrue(all(row["scientific_status"] == "warn" for row in rows))
        self.assertEqual(self.study.data["status"], "completed")
        reopened = open_study(self.study.folder)
        self.assertEqual(results_rows(reopened), rows)

    def test_failed_target_does_not_stop_another_target_or_drop_pairs(self):
        self.seed()
        with patch("stages.dock.DockStage", self.mock_dock(failed_target="Target A")):
            run_docking(self.study, exhaustiveness=1, n_decoys=1)
        rows = results_rows(self.study)
        self.assertEqual(len(rows), 4)
        self.assertEqual([r["status"] for r in rows], ["failed", "failed", "completed", "completed"])
        self.assertTrue(all("could not be docked" in row["error"] for row in rows[:2]))
        self.assertEqual(self.study.data["status"], "partial")

    def test_sanitized_or_duplicate_names_cannot_overwrite_drug_or_target_files(self):
        study = create_study(self.root, "Duplicates", [TargetSpec("A B", file=str(self.protein)), TargetSpec("A_B", file=str(self.protein))],
                             [DrugSpec("Same name", "CCO"), DrugSpec("Same name", "CCN")])
        self.assertEqual(len({t["project"] for t in study.data["targets"]}), 2)
        self.assertEqual(len({d["internal_name"] for d in study.data["drugs"]}), 2)
        self.assertEqual(len(results_rows(study)), 4)

    def test_changed_site_suppresses_previous_scores_and_poses(self):
        self.seed()
        with patch("stages.dock.DockStage", self.mock_dock()):
            run_docking(self.study, exhaustiveness=1, n_decoys=1)
        from stages.site import SiteStage
        with patch.object(SiteStage, "run", return_value=good({"box": {"center_x": 99}})):
            define_binding_sites(self.study)
        rows = results_rows(self.study)
        self.assertTrue(all(row["status"] == "needs_updating" for row in rows))
        self.assertTrue(all(row["score_kcal_mol"] == "" and row["pose_file"] == "" for row in rows))
        self.assertTrue(all(row["redock_rmsd_angstrom"] == "" and row["baseline_check"] == "not_current" for row in rows))

    def test_reopened_interrupted_run_does_not_present_previous_scores_or_checks_as_current(self):
        self.seed()
        with patch("stages.dock.DockStage", self.mock_dock()):
            run_docking(self.study, exhaustiveness=1, n_decoys=1)
        self.study.data["status"] = "docking"
        self.study.data["targets"][0]["status"] = "docking"
        self.study.save()
        reopened = open_study(self.study.folder)
        self.assertEqual(reopened.data["status"], "interrupted")
        rows = results_rows(reopened)
        self.assertEqual(rows[0]["status"], "interrupted")
        self.assertEqual(rows[0]["score_kcal_mol"], "")
        self.assertEqual(rows[0]["baseline_check"], "not_current")
        self.assertEqual(rows[2]["status"], "completed")

    def test_input_source_metadata_and_exact_download_links_survive_export(self):
        study = create_study(self.root, "Sources", [TargetSpec("EGFR", pdb_id="1M17")],
                             [DrugSpec("Ethanol", "CCO", cid=702, source="PubChem", source_url="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/702/property/ConnectivitySMILES/JSON", retrieved="2026-11-07")])
        reopened = open_study(study.folder)
        row = results_rows(reopened)[0]
        self.assertEqual(row["drug_cid"], 702)
        self.assertEqual(row["drug_source"], "PubChem")
        self.assertIn("/cid/702/", row["drug_download_url"])
        self.assertEqual(row["protein_download_url"], "https://files.rcsb.org/download/1M17.cif")
        self.assertIn("2026-11-07", json.dumps(reopened.data["drugs"]))

    def test_nonfatal_preparation_warning_proceeds_and_fatal_input_fails(self):
        self.seed()
        project = target_project(self.study, "target-001")
        prepared = json.loads(project.get_step("05_ligands")["values_json"])
        project.save_step("05_ligands", "Drugs", good(prepared, Verdict.WARN, blocks=True))
        project.close()
        project = target_project(self.study, "target-002")
        project.save_step("04_site", "Site", good({"box": {}}, Verdict.FAIL, blocks=True))
        project.close()
        with patch("stages.dock.DockStage", self.mock_dock()):
            run_docking(self.study, exhaustiveness=1, n_decoys=1)
        self.assertEqual([r["status"] for r in results_rows(self.study)], ["completed", "completed", "failed", "failed"])

    def test_export_preserves_pair_identity_failures_and_moved_pose_paths(self):
        self.seed()
        with patch("stages.dock.DockStage", self.mock_dock(failed_target="Target B")):
            run_docking(self.study, exhaustiveness=1, n_decoys=1)
        moved_folder = self.root / "Moved study"
        shutil.copytree(self.study.folder, moved_folder)
        moved = open_study(moved_folder)
        rows = list(csv.DictReader(io.StringIO(export_csv(moved).decode("utf-8-sig"))))
        self.assertEqual(len(rows), 4)
        self.assertEqual(sum(row["status"] == "failed" for row in rows), 2)
        self.assertTrue((moved_folder / rows[0]["pose_file"]).exists())
        with patch("core.references.render_references", return_value="Actual tool references"):
            archive = export_package(moved)
        with zipfile.ZipFile(archive) as package:
            self.assertIn("references.txt", package.namelist())
            self.assertIn("methods_and_checks.txt", package.namelist())
            self.assertIn(rows[0]["pose_file"], package.namelist())
            self.assertIn("Target run", package.read("results.csv").decode("utf-8-sig"))

    def test_offline_example_preparation_uses_writable_snapshot_and_fresh_drugs(self):
        from core.workshop import EXAMPLES
        source = EXAMPLES / "egfr" / "project.db"
        before = hashlib.sha256(source.read_bytes()).hexdigest()
        study = create_example_study(self.root)
        run_preparation(study)
        self.assertEqual(study.data["targets"][0]["status"], "prepared", study.data["targets"][0].get("error"))
        self.assertEqual(len(results_rows(study)), 3)
        self.assertTrue(all(row["status"] == "pending" for row in results_rows(study)))
        self.assertEqual(before, hashlib.sha256(source.read_bytes()).hexdigest())
        project = target_project(study, "target-001")
        try:
            self.assertIsNone(project.get_step("06_dock"))
            self.assertEqual(len(json.loads(project.get_step("05_ligands")["values_json"])["prepared"]), 3)
        finally:
            project.close()


if __name__ == "__main__":
    unittest.main()

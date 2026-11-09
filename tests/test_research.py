"""Selected-pair exports must preserve chemistry, provenance and limitations."""
import ast
import io
import json
from pathlib import Path
import tempfile
import unittest
import zipfile
from unittest.mock import patch

from core import research, study, viewer
from core.contract import StageResult
from stages.dock import DockStage, decoy_verdict


class ResearchExportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.study = study.create_example_study(Path(self.temp.name))
        study.run_preparation(self.study)

    def record_dock(self):
        project = study.target_project(self.study, self.study.data["targets"][0]["id"])
        self.addCleanup(project.close)
        drug = self.study.data["drugs"][1]
        prepared = next(x for x in viewer.values(project, "05_ligands")["prepared"] if x["name"] == drug["internal_name"])
        folder = project.step_folder("06_dock")
        ligand, problem = DockStage._ligand_pdbqt(prepared["smiles"], folder, prepared["name"], prepared_file=prepared["file"])
        self.assertFalse(problem)
        # A fixture pose retains Meeko's original connectivity and coordinate mapping.
        pose = folder / (ligand.stem + "_docked.pdbqt")
        pose.write_text("MODEL 1\nREMARK VINA RESULT: -6.86 0 0\n" + ligand.read_text() + "\nENDMDL\n")
        receptor = folder / "receptor.pdbqt"
        receptor.write_text("REMARK test receptor input\n")
        box = {"center_x": 1.0, "center_y": 2.0, "center_z": 3.0, "size_x": 15.0, "size_y": 16.0, "size_z": 17.0}
        project.save_step("06_dock", "Docking", StageResult(success=True,
            values={"receptor": str(receptor), "box": box, "settings": {"seed": 77, "exhaustiveness": 4, "n_poses": 9, "vina_version": "AutoDock Vina v1.2.7"},
                    "redock_rmsd": 0.503, "redock_molecule": drug["internal_name"],
                    "results": [{"name": drug["internal_name"], "smiles": prepared["smiles"], "score": -6.86, "file": str(pose), "input_file": str(ligand)}]},
            control=decoy_verdict(-6.86, [-7.50], drug["internal_name"])))
        self.project, self.drug, self.pose = project, drug, pose
        return project, drug, pose

    def archive(self, kind="docking", **kwargs):
        self.project, self.drug, self.pose = self.record_dock()
        blob = research.export_pair(self.study, self.study.data["targets"][0]["id"], self.drug["id"], kind, **kwargs)
        return zipfile.ZipFile(io.BytesIO(blob))

    def test_code_export_uses_exact_selected_input_and_settings(self):
        with self.archive() as archive:
            job = json.loads(archive.read("job.json"))
            self.assertEqual(job["drug"], self.drug["name"])
            self.assertEqual(job["settings"]["seed"], 77)
            self.assertEqual(job["settings"]["exhaustiveness"], 4)
            self.assertEqual(archive.read("recorded_pose.pdbqt"), self.pose.read_bytes())
            ast.parse(archive.read("run_docking.py").decode())
            notebook = json.loads(archive.read("Run_Docking.ipynb"))
            try:
                import nbformat
                nbformat.validate(nbformat.from_dict(notebook))
            except ImportError:
                self.assertEqual(notebook["nbformat"], 4)
            self.assertIn("not rerun", job["purpose"])

    def test_md_preserves_docked_chemistry_and_removes_reference(self):
        from rdkit import Chem
        with self.archive("md", system_type="water", nanoseconds=0.1) as archive:
            sdf = archive.read("ligand.sdf").decode()
            actual = Chem.MolFromMolBlock(sdf.split("$$$$")[0], removeHs=True)
            self.assertIsNotNone(actual)
            self.assertEqual(Chem.MolToSmiles(actual), Chem.MolToSmiles(Chem.MolFromSmiles(self.drug["smiles"])))
            protein = archive.read("protein.pdb").decode()
            label = viewer.values(self.project, "03_receptor")["ligand"]
            self.assertFalse(viewer.reference_atoms(protein, label))
            ast.parse(archive.read("run_md.py").decode())
            notebook = json.loads(archive.read("Run_MD.ipynb"))
            install = next(i for i,c in enumerate(notebook["cells"]) if "condacolab.install" in "".join(c["source"]))
            upload = next(i for i,c in enumerate(notebook["cells"]) if "files.upload" in "".join(c["source"]))
            self.assertLess(install, upload)

    def test_export_refuses_stale_pair_and_different_drug(self):
        project, drug, pose = self.record_dock()
        target_id = self.study.data["targets"][0]["id"]
        with self.assertRaisesRegex(ValueError, "current completed"):
            research.export_pair(self.study, target_id, self.study.data["drugs"][0]["id"])
        project.conn.execute("UPDATE steps SET status='stale' WHERE stage_id='06_dock'")
        project.conn.commit()
        with self.assertRaisesRegex(ValueError, "current completed"):
            research.export_pair(self.study, target_id, drug["id"])

    def test_report_keeps_negative_comparison_and_positive_pose_evidence(self):
        self.record_dock()
        text = research.report_text(self.study)
        self.assertIn("0.503", text)
        self.assertIn("-7.5", text)
        self.assertIn("not proven nonbinders", text)
        self.assertIn("does not overturn experimental evidence", text)
        self.assertNotIn("Do not report this molecule as a hit", text)
        self.assertIn("No SwissADME prediction was obtained automatically", text)

    def test_membrane_target_cannot_export_water_only_md(self):
        project, drug, pose = self.record_dock()
        target = self.study.data["targets"][0]
        target["name"] = "Human SV2A"
        with self.assertRaisesRegex(ValueError, "membrane protein"):
            research.export_pair(self.study, target["id"], drug["id"], "md", system_type="water")
        with zipfile.ZipFile(io.BytesIO(research.export_pair(self.study, target["id"], drug["id"], "md", system_type="membrane"))) as archive:
            self.assertIn("membrane_system.xml", archive.read("README.txt").decode())

    def test_prepared_drug_coordinates_are_used_without_reembedding(self):
        self.record_dock()
        drug = next(x for x in viewer.values(self.project, "05_ligands")["prepared"] if x["name"] == self.drug["internal_name"])
        with patch("rdkit.Chem.AllChem.EmbedMolecule", side_effect=AssertionError("Must use prepared coordinates")):
            path, error = DockStage._ligand_pdbqt(drug["smiles"], self.project.step_folder("06_dock"), "reused", prepared_file=drug["file"])
        self.assertIsNotNone(path)
        self.assertFalse(error)


if __name__ == "__main__":
    unittest.main()

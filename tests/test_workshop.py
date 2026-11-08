"""Behaviour checks for result safety, offline practice and project copies."""
import json
from copy import deepcopy
from pathlib import Path
import shutil
import tempfile
import unittest

from core.contract import ControlOutcome, StageResult, Verdict
from core.project import Project
from core.viewer import document, local_text, payload
from core.workshop import (EXAMPLES, MOLECULES, docking_checks, evidence_source, load_state, new_workshop,
                           save_state, teaching_check, use_saved_results, worksheet)


def result(value=1, verdict=Verdict.PASS):
    return StageResult(success=True, values={"value":value}, summary="Result",
                       control=ControlOutcome("Check",verdict,"Scientific check"))


class DependencyTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="insilico_workshop_test_")
        self.root = Path(self.tmp.name)
        self.p = Project.create(self.root/"project",name="test")

    def tearDown(self):
        self.p.close(); self.tmp.cleanup()

    def seed(self):
        for sid in ["01_target","02_structure","03_receptor","04_site","05_ligands","06_dock","08_admet","08_admet_swissadme","09_md"]:
            self.p.save_step(sid,sid,result())

    def test_ligand_change_invalidates_docking_admet_and_md_but_not_protein(self):
        self.seed()
        self.p.save_step("05_ligands","Molecules",result(2))
        for sid in ["06_dock","08_admet","08_admet_swissadme","09_md"]:
            self.assertEqual(self.p.get_step(sid)["status"],"stale")
            self.assertFalse(self.p.is_complete(sid))
        self.assertTrue(self.p.is_complete("03_receptor"))
        self.assertEqual(json.loads(self.p.get_step("06_dock")["values_json"]),{"value":1})

    def test_site_change_does_not_invalidate_independent_admet(self):
        self.seed(); self.p.save_step("04_site","Site",result(2))
        self.assertFalse(self.p.is_complete("06_dock"))
        self.assertTrue(self.p.is_complete("08_admet"))

    def test_unchanged_rerun_preserves_timestamp_and_dependents_after_reopen(self):
        self.seed()
        previous=self.p.get_step("05_ligands")["updated"]
        self.p.save_step("05_ligands","Molecules",result())
        self.assertEqual(self.p.get_step("05_ligands")["updated"],previous)
        folder=self.p.folder; self.p.close(); self.p=Project.open(folder)
        self.assertTrue(self.p.is_complete("06_dock"))

    def test_old_project_with_newer_input_is_audited_on_open(self):
        self.seed()
        self.p.conn.execute("UPDATE steps SET updated='2000-01-01T00:00:00+00:00' WHERE stage_id IN ('08_admet','09_md')")
        self.p.conn.commit(); folder=self.p.folder; self.p.close();self.p=Project.open(folder)
        self.assertFalse(self.p.is_complete("08_admet")); self.assertFalse(self.p.is_complete("09_md"))

    def test_finished_failure_is_not_green_and_unchecked_is_not_passed(self):
        self.p.save_step("05_ligands","Molecules",result(verdict=Verdict.FAIL))
        icon,colour,label=self.p.display_status("05_ligands")
        self.assertNotEqual(icon,"check_circle");self.assertIn("check failed",label)
        self.p.save_step("05_ligands","Molecules",StageResult(success=True))
        self.assertIn("not checked",self.p.display_status("05_ligands")[2])

    def test_whole_project_copy_rebases_file_references(self):
        f=self.p.step_folder("05_ligands")/"ligand.sdf"; f.write_text("actual data")
        self.p.save_step("05_ligands","Molecules",StageResult(success=True,values={"file":str(f)}))
        self.p.close(); copied=self.root/"copied";shutil.copytree(self.root/"project",copied)
        q=Project.open(copied)
        try:
            path=Path(json.loads(q.get_step("05_ligands")["values_json"])["file"])
            self.assertTrue(path.is_relative_to(copied));self.assertEqual(path.read_text(),"actual data")
        finally:q.close()


@unittest.skipUnless((EXAMPLES/"egfr"/"project.db").exists(),"Prepared workshop pack required")
class OfflineWorkshopTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix="insilico_workshop_test_")
        self.p=new_workshop(Path(self.tmp.name))

    def tearDown(self):
        self.p.close();self.tmp.cleanup()

    def test_fresh_exercise_has_protein_but_no_completed_docking(self):
        self.assertTrue(self.p.is_complete("04_site"));self.assertIsNone(self.p.get_step("06_dock"))
        data=payload(self.p,"02_structure")
        self.assertIn("ATOM",data["protein"]);self.assertIn("AQ4",data["reference"])
        self.assertIn("ATOM",payload(self.p,"04_site")["protein"])

    def test_saved_fallback_keeps_notes_and_uses_consistent_molecules(self):
        state=load_state(self.p);state["answers"]["welcome"]={"predict":"My prediction"}
        state["selected"]=["Erlotinib","Gefitinib"]
        state["done"]=["welcome","molecules","dock","conclude"]
        save_state(self.p,state)
        use_saved_results(self.p)
        state=load_state(self.p)
        self.assertEqual(state["answers"]["welcome"]["predict"],"My prediction")
        self.assertIn("not a live run",state["sources"]["06_dock"])
        self.assertEqual(state["selected"],[m["name"] for m in MOLECULES])
        self.assertEqual(state["done"],["welcome"])
        for sid,key in [("05_ligands","prepared"),("06_dock","results"),("08_admet","molecules")]:
            values=json.loads(self.p.get_step(sid)["values_json"])
            self.assertEqual({m["name"] for m in values[key]},{m["name"] for m in MOLECULES})
        self.assertIn("My prediction",worksheet(self.p))
        self.assertIn("Saved example",worksheet(self.p))

    def test_worksheet_exports_checks_numbers_activity_and_notes(self):
        use_saved_results(self.p)
        state=load_state(self.p)
        state["done"]=["interpret"]
        state["answers"]["interpret"]={"explain":"A model result needs experiments."}
        state["answers"]["next_experiment"]="Measure target binding."
        save_state(self.p,state)
        text=worksheet(self.p)
        values=json.loads(self.p.get_step("06_dock")["values_json"])
        for name,check in values["checks"].items():
            self.assertIn(name.upper()+" check: "+check["verdict"].upper(),text)
            self.assertIn(teaching_check("06_dock",check,name,values)["plain_verdict"],text)
        self.assertIn("SIZE check: NOT_APPLICABLE",text)
        self.assertIn(values["size_confound"]["why"],text)
        self.assertIn("Re-docking RMSD: "+str(values["redock_rmsd"])+" Å",text)
        for molecule in values["results"]:
            self.assertIn(molecule["name"]+" | "+str(molecule["score"]),text)
        self.assertIn("Calculated logP",text)
        self.assertIn("Activity: tried and discussed",text)
        self.assertIn("A model result needs experiments.",text)
        self.assertIn("Next experiment: Measure target binding.",text)

    def test_teaching_prose_preserves_recorded_results_without_exporting_old_claims(self):
        use_saved_results(self.p)
        step=self.p.get_step("06_dock")
        original_control=json.loads(step["control_json"])
        original_values=json.loads(step["values_json"])
        before=deepcopy((original_control,original_values))
        shown=teaching_check("06_dock",original_control,recorded=original_values)
        self.assertEqual(shown["verdict"],original_control["verdict"])
        for key,value in original_control["detail"].items():
            if not key.endswith(" check"):
                self.assertEqual(shown["detail"][key],value)
        checks=docking_checks(original_values,original_control)
        for name,check in original_values["checks"].items():
            self.assertEqual(checks[name]["verdict"],check["verdict"])
            self.assertEqual(checks[name]["detail"],check["detail"])
        self.assertEqual((original_control,original_values),before)
        self.assertEqual(self.p.get_step("06_dock"),step)
        text=worksheet(self.p)
        self.assertNotIn("settings are all working",text)
        self.assertNotIn("docking was done properly",text)
        self.assertNotIn("a margin above 1.0 is what you are aiming for",text)

    def test_preparation_teaching_preserves_records_and_avoids_historical_overclaims(self):
        for sid in ["01_target","02_structure","03_receptor","04_site"]:
            step=self.p.get_step(sid)
            original=json.loads(step["control_json"])
            before=deepcopy(original)
            shown=teaching_check(sid,original)
            self.assertEqual(original,before)
            self.assertEqual(self.p.get_step(sid),step)
            for key in ["name","verdict","detail","blocks_progress"]:
                self.assertEqual(shown[key],original[key])
            self.assertNotEqual(shown["plain_verdict"],original["plain_verdict"])
        text=worksheet(self.p)
        self.assertNotIn("single best predictor",text)
        self.assertNotIn("moving too much to be photographed",text)
        self.assertNotIn("nothing the protein needs was deleted",text)
        self.assertNotIn("would not hold a drug",text)
        self.assertNotIn("That is a snug fit",text)
        self.assertIn("do not predict success for a new molecule",text)
        self.assertIn("protonation",text)

    def test_legacy_check_interpretations_are_labelled_and_keep_numerical_evidence(self):
        p=Project.open(EXAMPLES/"sv2a")
        try:
            step=p.get_step("06_dock")
            control=json.loads(step["control_json"])
            values=json.loads(step["values_json"])
            before=deepcopy((control,values))
            checks=docking_checks(values,control)
            self.assertEqual(checks["baseline"]["verdict"],control["verdict"])
            self.assertEqual(checks["re-docking"]["verdict"],"warn")
            self.assertIn(str(values["redock_rmsd"]),checks["re-docking"]["plain_verdict"])
            self.assertIn("no separate historical verdict",checks["re-docking"]["plain_verdict"])
            self.assertIn("no separate historical verdict",checks["size"]["plain_verdict"])
            self.assertEqual(checks["size"]["detail"]["molecules compared"],values["size_confound"]["molecules"])
            self.assertEqual((control,values),before)
            self.assertEqual(p.get_step("06_dock"),step)
        finally:
            p.close()

    def test_stale_exports_and_source_labels_do_not_present_old_results_as_current(self):
        use_saved_results(self.p)
        old=json.loads(self.p.get_step("06_dock")["values_json"])
        self.p.save_step("05_ligands","Changed molecules",result(2))
        self.assertIn("Needs updating: inputs changed",evidence_source(self.p,"06_dock"))
        text=worksheet(self.p)
        self.assertIn("06_dock: Needs updating",text)
        self.assertIn("No current result. Rerun this step",text)
        self.assertNotIn("Re-docking RMSD: "+str(old["redock_rmsd"]),text)
        self.assertNotIn("Molecule | Docking score",text)
        self.assertNotIn("Molecule | MW",text)

    def test_saved_overlay_contains_one_predicted_erlotinib_pose_in_reference_frame(self):
        use_saved_results(self.p);data=payload(self.p,"06_dock","Erlotinib")
        self.assertTrue(data["overlay"])
        self.assertEqual(sum(line.startswith("HETATM") for line in data["reference"].splitlines()),29)
        self.assertEqual(sum(line.startswith("HETATM") for line in data["ligand"].splitlines()),29)
        self.assertTrue(data["box"]["size_x"]>0)
        self.assertFalse(payload(self.p,"06_dock","Risperidone")["overlay"])

    def test_copied_exercises_do_not_share_files_or_notes(self):
        q=new_workshop(Path(self.tmp.name))
        try:
            a=load_state(self.p);a["answers"]["welcome"]={"explain":"Only this pair"};save_state(self.p,a)
            self.assertNotEqual(self.p.folder,q.folder)
            self.assertNotIn("welcome",load_state(q)["answers"])
            qfile=Path(json.loads(q.get_step("03_receptor")["values_json"])["prepared_file"])
            self.assertTrue(qfile.is_relative_to(q.folder));self.assertTrue(qfile.exists())
        finally:q.close()

    def test_no_current_view_for_stale_results_and_no_files_outside_project(self):
        use_saved_results(self.p)
        self.p.save_step("05_ligands","Changed molecules",result(2))
        self.assertFalse(payload(self.p,"06_dock")["ligand"])
        self.assertEqual(local_text(self.p,str(EXAMPLES/"egfr"/"project.json")),"")


class ViewerSecurityTests(unittest.TestCase):
    def test_molecule_text_cannot_break_out_of_script(self):
        page=document({"name":"</script><img src=x onerror=alert(1)>"})
        self.assertNotIn("<img src=x",page)
        self.assertIn("/assets/vendor/3Dmol-min.js",page)
        self.assertNotIn("https://",page)


if __name__ == "__main__":
    unittest.main()

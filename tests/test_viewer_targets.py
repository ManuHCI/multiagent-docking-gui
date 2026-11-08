"""A familiar drug name cannot imply reference identity on another target."""
import json
from pathlib import Path
import tempfile
import unittest

from core.contract import StageResult
from core.viewer import payload
from core.workshop import new_workshop, use_saved_results


class TargetOverlayTests(unittest.TestCase):
    def test_overlay_requires_this_targets_recorded_redocking_identity(self):
        with tempfile.TemporaryDirectory(prefix="viewer_target_") as temporary:
            project = new_workshop(Path(temporary))
            try:
                use_saved_results(project)
                self.assertTrue(payload(project, "06_dock", "Erlotinib")["overlay"])
                row = project.get_step("06_dock")
                values = json.loads(row["values_json"])
                values["redock_molecule"] = "Another reference drug"
                project.save_step("06_dock", "Dock", StageResult(True, values=values))
                self.assertFalse(payload(project, "06_dock", "Erlotinib")["overlay"])
                values["redock_molecule"] = "Erlotinib"
                values["redock_rmsd"] = None
                project.save_step("06_dock", "Dock", StageResult(True, values=values))
                self.assertFalse(payload(project, "06_dock", "Erlotinib")["overlay"])
            finally:
                project.close()

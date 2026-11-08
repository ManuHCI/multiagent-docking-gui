"""Desktop storage isolation and resource-location checks."""
from __future__ import annotations

import json
import os
from pathlib import Path
import socket
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, call, patch

from core.runtime import configure, data_root, resource_root
from desktop import (desktop_origin, free_port, native_page_ready, native_page_state, server_process_options,
                     splash, stop_server, webview_backend)


class DesktopStorageTests(unittest.TestCase):
    def test_mac_projects_use_application_support_not_a_windows_environment_path(self):
        with patch("core.runtime.sys.platform", "darwin"), patch.dict(os.environ, {
                "INSILICO_DATA_DIR": "", "LOCALAPPDATA": "unrelated Windows path"}):
            self.assertEqual(data_root(), Path.home() / "Library" / "Application Support" / "InSilicoTeachingSuite")

    def test_mac_frozen_runtime_finds_private_plugins_and_rdkit_data(self):
        with tempfile.TemporaryDirectory(prefix="insilico_mac_libraries_") as temporary:
            root = Path(temporary)
            plugins = root / "lib" / "plugins"
            rdkit = root / "share" / "RDKit"
            plugins.mkdir(parents=True)
            (rdkit / "Data").mkdir(parents=True)
            with patch("core.runtime.sys.platform", "darwin"), patch.object(sys, "frozen", True, create=True), \
                    patch.object(sys, "_MEIPASS", str(root), create=True), patch.dict(os.environ):
                configure()
                self.assertEqual(os.environ["OPENMM_PLUGIN_DIR"], str(plugins))
                self.assertEqual(os.environ["RDBASE"], str(rdkit))
                self.assertEqual(os.environ["PATH"].split(os.pathsep)[0], str(root / "bin"))

    def test_explicit_user_data_folder_controls_settings_and_new_projects(self):
        with tempfile.TemporaryDirectory(prefix="insilico_desktop_settings_") as temporary:
            root = Path(temporary).resolve()
            env = os.environ.copy()
            env["INSILICO_DATA_DIR"] = str(root / "student data")
            env.pop("INSILICO_PROJECTS", None)
            script = '''import json
from core.settings import SETTINGS_FILE, projects_root, set_projects_root
initial = projects_root()
ok, message = set_projects_root(initial.parent / "chosen projects")
print(json.dumps({"settings":str(SETTINGS_FILE), "initial":str(initial),
                  "changed":str(projects_root()), "ok":ok}))'''
            result = subprocess.run([sys.executable, "-c", script], env=env,
                                    cwd=Path(__file__).resolve().parent.parent,
                                    capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)
            values = json.loads(result.stdout)
            self.assertEqual(Path(values["settings"]), root / "student data" / "settings.json")
            self.assertEqual(Path(values["initial"]), root / "student data" / "Projects")
            self.assertTrue(values["ok"])
            self.assertEqual(Path(values["changed"]), root / "student data" / "chosen projects")
            saved = json.loads(Path(values["settings"]).read_text(encoding="utf-8"))
            self.assertEqual(Path(saved["projects_root"]), Path(values["changed"]))

    def test_packaged_resources_are_independent_of_working_directory(self):
        with tempfile.TemporaryDirectory(prefix="insilico_resources_") as temporary:
            bundled = Path(temporary).resolve() / "_internal"
            with patch.object(sys, "_MEIPASS", str(bundled), create=True):
                self.assertEqual(resource_root(), bundled)
            self.assertEqual(resource_root(), Path(__file__).resolve().parent.parent)

    def test_data_folder_override_is_resolved_without_reading_developer_settings(self):
        with tempfile.TemporaryDirectory(prefix="insilico_data_") as temporary:
            with patch.dict(os.environ, {"INSILICO_DATA_DIR": temporary}):
                self.assertEqual(data_root(), Path(temporary).resolve())


class DesktopStartupTests(unittest.TestCase):
    def test_native_origin_stays_local_without_a_webkit_network_exception(self):
        with patch("desktop.sys.platform", "darwin"):
            self.assertEqual(desktop_origin(43210), "http://localhost:43210")
        with patch("desktop.sys.platform", "win32"):
            self.assertEqual(desktop_origin(43210), "http://127.0.0.1:43210")

    def test_mac_probe_reads_json_string_and_requires_the_working_home_page(self):
        window = MagicMock()
        state = {"url": "http://127.0.0.1:8111/study", "ready_state": "complete",
                 "navigation": ["Home", "Projects", "Reference", "About"], "start_project": True,
                 "product_title": False}
        window.run_js.return_value = json.dumps(state)
        with patch("desktop.sys.platform", "darwin"):
            self.assertEqual(native_page_state(window), state)
        window.evaluate_js.assert_not_called()
        self.assertTrue(native_page_ready(state))
        self.assertFalse(native_page_ready({**state, "url": "about:blank"}))
        self.assertFalse(native_page_ready({**state, "navigation": []}))
        self.assertFalse(native_page_ready({**state, "start_project": False}))

    def test_native_renderer_matches_operating_system(self):
        with patch("desktop.sys.platform", "darwin"):
            self.assertEqual(webview_backend(), ("cocoa", "Cocoa WebKit"))
            self.assertEqual(server_process_options(), {"start_new_session": True})
        with patch("desktop.sys.platform", "win32"):
            self.assertEqual(webview_backend(), ("edgechromium", "Edge WebView2"))
            self.assertIn("creationflags", server_process_options())
            self.assertNotIn("start_new_session", server_process_options())

    def test_mac_cleanup_stops_the_owned_group_even_if_server_exited_first(self):
        process = MagicMock(pid=4321)
        process.poll.return_value = 0
        with patch("desktop.sys.platform", "darwin"), patch("desktop.os.killpg", create=True) as kill_group, \
                patch("desktop.signal.SIGKILL", 9, create=True):
            stop_server(process)
        self.assertEqual(kill_group.call_args_list, [call(4321, 15), call(4321, 9)])
        process.terminate.assert_not_called()
        process.kill.assert_not_called()

    def test_port_selection_does_not_use_an_occupied_server_port(self):
        with socket.socket() as existing:
            existing.bind(("127.0.0.1", 0))
            existing.listen()
            selected = free_port()
            self.assertNotEqual(selected, existing.getsockname()[1])
            with socket.socket() as candidate:
                candidate.bind(("127.0.0.1", selected))

    def test_startup_failure_text_is_escaped_and_has_no_busy_animation(self):
        page = splash("Cannot open <student> & notes", error=True)
        self.assertIn("Cannot open &lt;student&gt; &amp; notes", page)
        self.assertNotIn('<div class="bar">', page)
        self.assertIn("MEDICAL EDUCATION", page)


class PortableDockingTests(unittest.TestCase):
    def test_mac_bundle_prefers_its_native_vina_and_never_windows_exe(self):
        from stages.dock import find_vina
        with tempfile.TemporaryDirectory(prefix="insilico_mac_vina_") as temporary:
            root = Path(temporary)
            (root / "bin").mkdir()
            native = root / "bin" / "vina"
            native.write_bytes(b"native build fixture")
            with patch("stages.dock.sys.platform", "darwin"), patch.object(sys, "frozen", True, create=True), \
                    patch("core.runtime.resource_root", return_value=root), patch("stages.dock.os.access", return_value=True), \
                    patch("stages.dock.shutil.which", return_value="/usr/local/bin/vina"):
                self.assertEqual(find_vina(), str(native))
            native.unlink()
            (root / "vina.exe").write_bytes(b"Windows fixture")
            with patch("stages.dock.sys.platform", "darwin"), patch.object(sys, "frozen", True, create=True), \
                    patch("core.runtime.resource_root", return_value=root), patch("stages.dock.shutil.which", return_value=None):
                self.assertIsNone(find_vina())

    def test_windows_recorded_project_paths_reopen_as_native_paths(self):
        from core.contract import StageResult
        from core.project import Project
        with tempfile.TemporaryDirectory(prefix="insilico_portable_paths_") as temporary:
            folder = Path(temporary) / "Copied project"
            project = Project.create(folder, "Copied project")
            artifact = project.step_folder("03_receptor") / "protein.pdbqt"
            artifact.write_text("ATOM\n", encoding="utf-8")
            project.save_step("03_receptor", "Prepared protein", StageResult(True, values={
                "prepared_file": r"E:\Original\Study\steps\03_receptor\protein.pdbqt",
                "external_file": r"E:\Original\Study sibling\keep.pdb"}))
            summary = project.summary
            summary["location"] = r"E:\Original\Study"
            project._write_summary(summary)
            project.close()
            reopened = Project.open(folder)
            try:
                values = json.loads(reopened.get_step("03_receptor")["values_json"])
                self.assertEqual(Path(values["prepared_file"]), artifact)
                self.assertTrue(Path(values["prepared_file"]).is_file())
                self.assertEqual(values["external_file"], r"E:\Original\Study sibling\keep.pdb")
            finally:
                reopened.close()


if __name__ == "__main__":
    unittest.main()

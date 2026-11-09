"""Install, exercise and uninstall the setup in an isolated workspace folder.

This checks the installer on the build computer; it is not a clean-PC test.
No desktop shortcuts or uninstall-registry entries are created in this mode.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / ".verification" / "studio"


def main() -> int:
    output = OUTPUT.resolve()
    output.mkdir(parents=True, exist_ok=True)
    parent = ROOT / ".qa-temp"
    parent.mkdir(exist_ok=True)
    case = Path(tempfile.mkdtemp(prefix="setup-acceptance-", dir=parent)).resolve()
    installed = case / "Teaching Suite"
    if not installed.resolve().is_relative_to(parent.resolve()):
        raise RuntimeError("Installation test must stay in the workspace QA directory.")
    from desktop import VERSION
    setup = ROOT / "dist" / f"InSilicoDockingStudio-{VERSION}-Windows-Setup.exe"
    # NSIS requires /D last, without quotes around its value, even with spaces.
    command = subprocess.list2cmdline([str(setup), "/S", "/NOINTEGRATION"]) + " /D=" + str(installed)
    subprocess.run(command, cwd=ROOT, timeout=180, check=True)
    manifest = json.loads((installed / "build-manifest.json").read_text(encoding="utf-8"))
    for relative, record in manifest["files"].items():
        path = installed / relative
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != record["sha256"]:
            raise AssertionError(f"Installed file differs: {relative}")
    environment = os.environ.copy()
    for name in ("PYTHONHOME", "PYTHONPATH", "CONDA_PREFIX", "CONDA_DEFAULT_ENV", "RDBASE", "OPENMM_PLUGIN_DIR", "INSILICO_PROJECTS"):
        environment.pop(name, None)
    environment["PATH"] = r"C:\Windows\System32;C:\Windows"
    user_data = case / "student data"
    environment["INSILICO_DATA_DIR"] = str(user_data)
    notes = user_data / "Projects" / "keep-my-notes.txt"
    notes.parent.mkdir(parents=True)
    notes.write_text("Student notes must survive uninstall.\n", encoding="utf-8")
    executable = installed / "InSilicoTeachingSuite.exe"
    for flag, filename in [("--verify-output", "installed-diagnostics.json"),
                           ("--verify-science", "installed-science.json"),
                           ("--verify-study", "installed-batch.json"),
                           ("--window-test-output", "installed-native-window.json")]:
        completed = subprocess.run([str(executable), flag, str(output / filename)], cwd=case,
                                   env=environment, check=False, timeout=180)
        report_path = output / filename
        if not report_path.is_file():
            raise AssertionError(f"Installed check wrote no report: {flag}; exit code {completed.returncode}")
        result = json.loads(report_path.read_text(encoding="utf-8"))
        print(f"Installed check {flag}: {json.dumps(result, indent=2)}", flush=True)
        completed.check_returncode()
        if not result.get("can_start", result.get("passed", False)):
            raise AssertionError(f"Installed acceptance failed: {filename}")
        if flag == "--window-test-output" and not (result.get("downloads_enabled") and result.get("server_stopped")):
            raise AssertionError("Native downloads or server shutdown is unavailable")
    subprocess.run([sys.executable, str(ROOT / "verify_live_docking.py"),
                    "--resources", str(installed / "_internal"),
                    "--output", str(output / "installed-live-docking.json")],
                   cwd=ROOT, env=environment, check=True, timeout=120)
    subprocess.run([sys.executable, str(ROOT / "verify_desktop_ui.py"),
                    "--executable", str(executable)], cwd=ROOT, check=True, timeout=300)
    ui_report = json.loads((output / "ui-acceptance.json").read_text(encoding="utf-8"))
    study_data = (Path(ui_report["user_data_dir"]) / "Projects").resolve()
    if not study_data.is_relative_to((ROOT / ".qa-temp").resolve()):
        raise AssertionError("The UI test's saved studies are outside the isolated QA directory")
    study_manifests = list(study_data.glob("*/study.json"))
    if not study_manifests:
        raise AssertionError("The installed UI did not persist a docking study")
    saved_hashes = {path: hashlib.sha256(path.read_bytes()).hexdigest()
                    for path in study_data.rglob("*")
                    if path.is_file() and path.suffix.lower() in {".json", ".db", ".pdbqt", ".sdf"}}
    for relative, record in manifest["files"].items():
        if relative.startswith("_internal/assets/workshop/"):
            if hashlib.sha256((installed / relative).read_bytes()).hexdigest() != record["sha256"]:
                raise AssertionError(f"The app modified a bundled reference: {relative}")
    subprocess.run([str(installed / "Uninstall.exe"), "/S"], cwd=case,
                   env=environment, timeout=120, check=True)
    deadline = time.monotonic() + 20
    while executable.exists() and time.monotonic() < deadline:
        time.sleep(0.2)
    if executable.exists():
        raise AssertionError("Uninstall did not remove the application executable.")
    if not notes.exists() or notes.read_text(encoding="utf-8") != "Student notes must survive uninstall.\n":
        raise AssertionError("Uninstall changed student notes.")
    if any(not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest
           for path, digest in saved_hashes.items()):
        raise AssertionError("Uninstall changed a saved study or its molecular results")
    (output / "setup-acceptance.json").write_text(json.dumps({
        "passed": True, "installed_files_verified": len(manifest["files"]),
        "installation_with_spaces": True, "system_python_path_removed": True,
        "native_window": True, "uninstall_kept_notes": True,
        "saved_studies_preserved": len(study_manifests), "saved_files_unchanged": len(saved_hashes),
        "bundled_references_unchanged": True,
        "scope": "isolated setup test on the build computer; separate clean-PC test still required",
    }, indent=2), encoding="utf-8")
    print(f"Setup install, application and uninstall checks passed: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

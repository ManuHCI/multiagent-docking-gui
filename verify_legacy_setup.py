"""Isolated legacy installer/upgrade/uninstall acceptance on the current host.

All application and data paths belong to a new workspace .qa-temp case.
/NOINTEGRATION avoids real shortcuts and uninstall registry entries. This
does not claim execution on original Windows 7, 8 or 8.1 machines.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import time
import traceback

ROOT = Path(__file__).resolve().parent
DEFAULT_SETUP = ROOT / "dist/InSilicoDockingStudio-1.3.0-Legacy-Windows-x64-Setup.exe"
DEFAULT_OUTPUT = ROOT / ".verification/studio/legacy-setup-acceptance.json"


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def host_description():
    record = {"platform": platform.platform(), "machine": platform.machine(), "python": sys.version}
    if os.name == "nt":
        version = sys.getwindowsversion()
        record.update(major=version.major, minor=version.minor, build=version.build)
        record["windows_family"] = "Windows 11" if version.major == 10 and version.build >= 22000 else "Windows " + str(version.major) + "." + str(version.minor)
    return record


def case_path(folder):
    """Keep every recursive inspection and installer destination in this case."""
    folder = Path(folder).resolve()
    folder.relative_to((ROOT / ".qa-temp").resolve())
    return folder


def installed_hashes(installed):
    installed = case_path(installed)
    return {str(path.relative_to(installed)): sha256(path)
            for path in installed.rglob("*") if path.is_file()}


def verify_installed_files(installed, expected):
    installed = case_path(installed)
    for relative, record in expected["files"].items():
        path = (installed / relative).resolve()
        path.relative_to(installed)
        if not path.is_file() or sha256(path) != record["sha256"]:
            raise AssertionError("Installed payload mismatch: " + relative)
    # Prove the legacy app contains the pinned older runtimes rather than
    # resolving to whichever VC/UCRT version is installed on this host.
    files_by_name = {Path(relative).name.lower(): relative for relative in expected["files"]}
    runtime_count = 0
    for metadata_file in ("runtime-files.json", "ucrt-files.json"):
        records = json.loads((ROOT / ".build-tools/legacy-vc-runtime" / metadata_file).read_text(encoding="utf-8-sig"))
        for record in records:
            relative = files_by_name.get(record["name"].lower())
            if not relative or sha256(installed / relative) != record["sha256"]:
                raise AssertionError("Pinned Microsoft runtime missing or changed: " + record["name"])
            runtime_count += 1
    return {"payload_files": len(expected["files"]), "pinned_runtime_files": runtime_count}


def stop_owned_idle(process):
    """Stop only the idle app launched by this test, never an active Vina job."""
    if process is None:
        return
    try:
        import psutil
        children = psutil.Process(process.pid).children(recursive=True)
        if any(child.name().lower().startswith("vina") for child in children):
            raise RuntimeError("The test app has an active Vina calculation; refusing to stop it.")
    except ImportError:
        pass
    except Exception as error:
        # An already exited root is harmless. Other failures must not turn
        # into a blind process-tree kill.
        if process.poll() is None:
            raise error
    if process.poll() is None:
        process.terminate()
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait(timeout=10)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--installer", type=Path, default=DEFAULT_SETUP)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--skip-science", action="store_true", help="Installer-only developer check; omits repeated installed science/batch checks.")
    args = parser.parse_args()
    started = time.monotonic()
    report = {"passed": False, "scope": "Isolated current-host legacy setup acceptance",
              "host": host_description(), "old_windows_execution_tested": False,
              "all_target_operating_systems_validated": False,
              "validation_limit": "Current-host execution and static checks do not certify Windows 7/8/8.1.",
              "modern_installation_touched": False, "integration_disabled": True}
    idle = None
    idle_log = None
    installed = None
    environment = None
    try:
        if os.name != "nt":
            raise RuntimeError("Run this installer test on Windows.")
        setup = args.installer.resolve(strict=True)
        expected = json.loads((ROOT / "dist/InSilicoDockingStudioLegacy/build-manifest.json").read_text(encoding="utf-8"))
        expected_exe_hash = expected["files"]["InSilicoDockingStudioLegacy.exe"]["sha256"]
        report.update(installer=str(setup), installer_sha256=sha256(setup), application_version=expected["version"])
        qa = ROOT / ".qa-temp"
        qa.mkdir(exist_ok=True)
        case = case_path(Path(tempfile.mkdtemp(prefix="legacy-setup-acceptance-", dir=str(qa))))
        installed = case / "Legacy Teaching Suite"
        executable = installed / "InSilicoDockingStudioLegacy.exe"
        data = case / "student data"
        data.mkdir()
        project_root = data / "Projects"
        project_root.mkdir()
        marker = data / "saved-notes-marker.txt"
        marker.write_text("Saved student work must survive upgrade and uninstall.\n", encoding="utf-8")
        # A real project folder, not only a loose marker, tests preservation.
        from legacy.compat import install_compatibility
        install_compatibility()
        from core.project import Project
        project = Project.create(project_root / "Preserved teaching project", "Preserved teaching project")
        project.close()
        preserved = {str(path): sha256(path) for path in data.rglob("*") if path.is_file()}
        report.update(test_directory=str(case), project_files_before=len(preserved))
        environment = os.environ.copy()
        for name in ("PYTHONHOME", "PYTHONPATH", "CONDA_PREFIX", "CONDA_DEFAULT_ENV", "RDBASE", "OPENMM_PLUGIN_DIR"):
            environment.pop(name, None)
        windows = Path(environment.get("SystemRoot", "C:/Windows"))
        environment["PATH"] = os.pathsep.join(str(path) for path in (windows / "System32", windows, windows / "System32/Wbem"))
        environment["INSILICO_DATA_DIR"] = str(data)
        environment["INSILICO_PROJECTS"] = str(project_root)
        # NSIS /D must be last, and its directory argument must be unquoted.
        command = subprocess.list2cmdline([str(setup), "/S", "/NOINTEGRATION"]) + " /D=" + str(installed)
        subprocess.run(command, cwd=case, env=environment, timeout=240, check=True)
        report["fresh_install"] = verify_installed_files(installed, expected)
        ini = installed / "installation.ini"
        if "NoIntegration=1" not in ini.read_text(encoding="utf-16" if ini.read_bytes().startswith(b"\xff\xfe") else "utf-8-sig"):
            raise AssertionError("The isolated test installation did not disable integration.")

        # Mark the installed manifest so the successful retry proves a real
        # replacement occurred, while a rejected in-use attempt changes none.
        manifest_file = installed / "build-manifest.json"
        marked = json.loads(manifest_file.read_text(encoding="utf-8"))
        marked["upgrade_test_marker"] = "Must survive rejection and disappear after successful retry."
        manifest_file.write_text(json.dumps(marked, indent=2), encoding="utf-8")
        before = installed_hashes(installed)
        idle_log = (case / "idle-native.log").open("w", encoding="utf-8")
        idle = subprocess.Popen([str(executable)], cwd=case, env=environment,
                                stdout=idle_log, stderr=subprocess.STDOUT,
                                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        time.sleep(2)
        if idle.poll() is not None:
            raise RuntimeError("The installed idle native app exited before the in-use guard test.")
        rejected = subprocess.run(command, cwd=case, env=environment, timeout=60, check=False)
        if rejected.returncode != 10:
            raise AssertionError("In-use legacy update should reject with exit 10; got " + str(rejected.returncode))
        if installed_hashes(installed) != before:
            raise AssertionError("A rejected update changed installed application files or metadata.")
        if idle.poll() is not None:
            raise AssertionError("A rejected update stopped the running test app.")
        stop_owned_idle(idle)
        idle = None
        idle_log.close()
        idle_log = None
        deadline = time.monotonic() + 15
        retries = 0
        while True:
            result = subprocess.run(command, cwd=case, env=environment, timeout=240, check=False)
            if result.returncode == 0:
                break
            if result.returncode != 10 or time.monotonic() >= deadline:
                raise AssertionError("Stopped-app upgrade retry failed: " + str(result.returncode))
            retries += 1
            time.sleep(.3)
        report["upgraded_install"] = verify_installed_files(installed, expected)
        if "upgrade_test_marker" in json.loads(manifest_file.read_text(encoding="utf-8")):
            raise AssertionError("The successful upgrade did not replace the old manifest.")
        if sha256(executable) != expected_exe_hash:
            raise AssertionError("The installed executable differs from the verified build.")
        report.update(in_use_upgrade_rejected_exit=10, rejected_update_changed_no_files=True,
                      running_app_remained_alive=True, successful_retry=True, transient_lock_retries=retries)

        checks = [("diagnostics", "--verify-output", "can_start", 180),
                  ("native-window", "--window-test-output", "passed", 180)]
        if not args.skip_science:
            checks[1:1] = [("science", "--verify-science", "passed", 240),
                           ("multi-target-docking", "--verify-study", "passed", 300)]
        report["installed_checks"] = {}
        for name, flag, success_key, timeout in checks:
            destination = case / (name + ".json")
            subprocess.run([str(executable), flag, str(destination)], cwd=case, env=environment,
                           timeout=timeout, check=True)
            result = json.loads(destination.read_text(encoding="utf-8"))
            if not result.get(success_key):
                raise AssertionError("Installed legacy check failed: " + name)
            report["installed_checks"][name] = result

        # NSIS has an explicit generated file list and never deletes project
        # storage. Its temporary uninstaller copy may outlive the first process.
        subprocess.run([str(installed / "Uninstall.exe"), "/S"], cwd=case, env=environment,
                       timeout=120, check=True)
        deadline = time.monotonic() + 30
        while executable.exists() and time.monotonic() < deadline:
            time.sleep(.2)
        if executable.exists():
            raise AssertionError("The isolated legacy application did not uninstall.")
        for path, expected_hash in preserved.items():
            if not Path(path).is_file() or sha256(Path(path)) != expected_hash:
                raise AssertionError("A saved project or notes file changed: " + path)
        report.update(passed=True, uninstall_removed_application=True, user_data_preserved=True,
                      preserved_project_files=len(preserved), executable_sha256=expected_exe_hash)
        print("Legacy isolated setup, in-use upgrade, retry, installed checks and saved-data preservation passed.", flush=True)
        return 0
    except Exception as error:
        report.update(error=type(error).__name__ + ": " + str(error), traceback=traceback.format_exc())
        print(report["error"], file=sys.stderr, flush=True)
        return 1
    finally:
        if idle is not None:
            try:
                stop_owned_idle(idle)
            except Exception as error:
                report["idle_cleanup_error"] = str(error)
        if idle_log is not None:
            idle_log.close()
        report["seconds"] = round(time.monotonic() - started, 3)
        args.output.resolve().parent.mkdir(parents=True, exist_ok=True)
        args.output.resolve().write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())

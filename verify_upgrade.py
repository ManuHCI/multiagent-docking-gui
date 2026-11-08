"""Reject an in-use upgrade before copying, then verify a successful retry.

All installation, server and data paths are isolated under this workspace's
.qa-temp directory. The test never changes the real installed application.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import socket
import subprocess
import sys
import tempfile
import time
import traceback
import urllib.request

ROOT = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stop_owned_server(server, report: dict) -> None:
    """Stop only this test's server tree, including frozen helper processes."""
    try:
        import psutil
        children = psutil.Process(server.pid).children(recursive=True)
        report["owned_server_children_before_stop"] = [child.pid for child in children]
    except (ImportError, OSError):
        children = []
    except Exception:
        children = []
    if server.poll() is None:
        if os.name == "nt":
            taskkill = Path(os.environ.get("SystemRoot", r"C:\Windows")) / "System32/taskkill.exe"
            stopped = subprocess.run([str(taskkill), "/PID", str(server.pid), "/T", "/F"],
                                     capture_output=True, timeout=15,
                                     creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            report["owned_tree_stop_exit"] = stopped.returncode
        if server.poll() is None:
            server.terminate()
        try:
            server.wait(timeout=15)
        except subprocess.TimeoutExpired:
            server.kill()
            server.wait(timeout=10)
    for child in children:
        try:
            child.wait(timeout=10)
        except Exception:
            # Only descendants captured while the owned root was alive.
            if child.is_running():
                child.kill()
                child.wait(timeout=10)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--installer", type=Path,
                        default=ROOT / "dist/InSilicoTeachingSuite-1.3.0-Setup.exe")
    parser.add_argument("--output", type=Path,
                        default=ROOT / ".verification/studio/upgrade-acceptance.json")
    arguments = parser.parse_args()
    report = {"passed": False, "scope": "isolated upgrade test on the build computer",
              "setup_version": "1.3.0", "application_payload_version": "1.3.0"}
    started = time.monotonic()
    server = None
    log = None
    installed = None
    try:
        setup = arguments.installer.resolve()
        if not setup.is_file():
            raise FileNotFoundError(f"The guarded setup is missing: {setup}")
        parent = (ROOT / ".qa-temp").resolve()
        parent.mkdir(exist_ok=True)
        case = Path(tempfile.mkdtemp(prefix="upgrade-acceptance-", dir=parent)).resolve()
        report["test_directory"] = str(case)
        report["installer_sha256"] = sha256(setup)
        installed = case / "Teaching Suite"
        if not installed.resolve().is_relative_to(parent):
            raise RuntimeError("The upgrade test must remain inside the workspace QA folder")
        command = subprocess.list2cmdline([str(setup), "/S", "/NOINTEGRATION"]) + " /D=" + str(installed)
        subprocess.run(command, cwd=case, timeout=180, check=True)
        executable = installed / "InSilicoTeachingSuite.exe"
        expected = sha256(ROOT / "dist/InSilicoTeachingSuite/InSilicoTeachingSuite.exe")
        if sha256(executable) != expected:
            raise AssertionError("Fresh installation executable differs from the verified payload")
        manifest_path = installed / "build-manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["upgrade_test_marker"] = "This must remain unchanged while the app is in use."
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        watched = [manifest_path, executable, installed / "installation.ini", installed / "Getting Started.md"]
        before = {str(path.relative_to(installed)): sha256(path) for path in watched}
        environment = os.environ.copy()
        for name in ("PYTHONPATH", "PYTHONHOME", "CONDA_PREFIX", "CONDA_DEFAULT_ENV", "RDBASE", "OPENMM_PLUGIN_DIR", "INSILICO_PROJECTS"):
            environment.pop(name, None)
        environment["PATH"] = r"C:\Windows\System32;C:\Windows"
        data = case / "student data"
        environment["INSILICO_DATA_DIR"] = str(data)
        data.mkdir()
        saved = data / "saved-study-marker.txt"
        saved.write_text("Existing study data stays outside the application folder.\n", encoding="utf-8")
        with socket.socket() as selected:
            selected.bind(("127.0.0.1", 0))
            port = selected.getsockname()[1]
        base = f"http://127.0.0.1:{port}"
        log = (case / "server.log").open("w", encoding="utf-8")
        server = subprocess.Popen([str(executable), "--serve", "--port", str(port)],
                                  cwd=case, env=environment, stdout=log, stderr=log,
                                  creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        deadline = time.monotonic() + 90
        while True:
            if server.poll() is not None:
                raise RuntimeError(f"The installed test server exited with {server.returncode}")
            try:
                with urllib.request.urlopen(base + "/assets/logo.png", timeout=1) as response:
                    if response.status == 200:
                        break
            except OSError:
                if time.monotonic() >= deadline:
                    raise TimeoutError("The installed test server did not become ready")
                time.sleep(0.2)
        rejected = subprocess.run(command, cwd=case, timeout=60, check=False)
        if rejected.returncode != 10:
            raise AssertionError(f"In-use upgrade should exit before copying with code 10, got {rejected.returncode}")
        after = {str(path.relative_to(installed)): sha256(path) for path in watched}
        if after != before:
            raise AssertionError("A blocked upgrade changed application bytes or installation metadata")
        with urllib.request.urlopen(base + "/assets/logo.png", timeout=2) as response:
            if response.status != 200:
                raise AssertionError("The blocked upgrade interrupted the running application")
        stop_owned_server(server, report)
        server = None
        log.close()
        log = None
        try:
            with urllib.request.urlopen(base + "/assets/logo.png", timeout=1):
                raise AssertionError("The owned server remained available after tree shutdown")
        except OSError:
            pass
        # A just-exited Windows process can retain its image mapping briefly.
        # A conservative guard may reject that state; wait only for lock release.
        retry_deadline = time.monotonic() + 10
        lock_retries = 0
        while True:
            retried = subprocess.run(command, cwd=case, timeout=180, check=False)
            if retried.returncode == 0:
                break
            if retried.returncode != 10 or time.monotonic() >= retry_deadline:
                raise AssertionError(f"Upgrade retry failed after owned-server shutdown: {retried.returncode}")
            lock_retries += 1
            time.sleep(0.25)
        report["transient_lock_retries"] = lock_retries
        if sha256(executable) != expected:
            raise AssertionError("The retried upgrade did not install the verified executable")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if "upgrade_test_marker" in manifest:
            raise AssertionError("The successful retry left the old marked manifest behind")
        for relative, record in manifest["files"].items():
            path = (installed / relative).resolve()
            if not path.is_relative_to(installed) or not path.is_file() or sha256(path) != record["sha256"]:
                raise AssertionError(f"Retried upgrade payload mismatch: {relative}")
        native_report = case / "native-window.json"
        subprocess.run([str(executable), "--window-test-output", str(native_report)],
                       cwd=case, env=environment, timeout=180, check=True)
        native = json.loads(native_report.read_text(encoding="utf-8"))
        if not (native.get("passed") and native.get("route") == "/study" and native.get("server_stopped")):
            raise AssertionError("The upgraded native app did not load the current workflow cleanly")
        subprocess.run([str(installed / "Uninstall.exe"), "/S"], cwd=case,
                       env=environment, timeout=120, check=True)
        deadline = time.monotonic() + 20
        while executable.exists() and time.monotonic() < deadline:
            time.sleep(0.2)
        if executable.exists():
            raise AssertionError("The isolated test application did not uninstall")
        if not saved.is_file():
            raise AssertionError("The upgrade/uninstall changed saved user data")
        report.update(passed=True, rejected_exit_code=10,
                      no_copy_or_metadata_changes_while_running=True,
                      running_app_remained_available=True, successful_retry=True,
                      installed_files_verified=len(manifest["files"]),
                      executable_sha256=expected, native_window=native,
                      user_data_preserved=True)
        print("Guarded in-use upgrade and successful retry checks passed.", flush=True)
        return 0
    except Exception as error:
        report.update(error=f"{type(error).__name__}: {error}", traceback=traceback.format_exc())
        print(report["error"], file=sys.stderr, flush=True)
        return 1
    finally:
        if server is not None and server.poll() is None:
            try:
                stop_owned_server(server, report)
            except Exception:
                if server.poll() is None:
                    server.kill()
        if log is not None:
            log.close()
        report["seconds"] = round(time.monotonic() - started, 3)
        arguments.output.resolve().parent.mkdir(parents=True, exist_ok=True)
        arguments.output.resolve().write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())

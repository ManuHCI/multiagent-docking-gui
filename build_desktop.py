"""Build and check the self-contained Windows desktop distribution.

Use the Python environment containing the working scientific dependencies:

    python build_desktop.py
    python build_desktop.py --zip

PyInstaller and pywebview may be installed in the project-local .build-tools
directory. End users do not need that directory or a Python installation.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime, timezone


ROOT = Path(__file__).resolve().parent
TOOLS = ROOT / ".build-tools"
WORK = ROOT / "build" / "desktop"
DIST = ROOT / "dist"
PAYLOAD = DIST / "InSilicoTeachingSuite"
EXECUTABLE = PAYLOAD / "InSilicoTeachingSuite.exe"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_environment() -> dict[str, str]:
    environment = os.environ.copy()
    paths = [str(TOOLS), str(ROOT)]
    if environment.get("PYTHONPATH"):
        paths.append(environment["PYTHONPATH"])
    environment["PYTHONPATH"] = os.pathsep.join(paths)
    environment["PYTHONNOUSERSITE"] = "1"
    environment["MATPLOTLIB"] = "false"
    environment["PYINSTALLER_CONFIG_DIR"] = str(WORK / "cache")
    windows = Path(environment.get("SystemRoot", r"C:\Windows"))
    environment["PATH"] = os.pathsep.join(str(p) for p in
        (Path(sys.prefix), Path(sys.prefix) / "Scripts", Path(sys.prefix) / "Library/bin",
         windows / "System32", windows, windows / "System32/Wbem"))
    return environment


def check_prerequisites() -> None:
    if sys.platform != "win32":
        raise RuntimeError("Build this Windows executable on Windows.")
    if TOOLS.is_dir():
        sys.path.insert(0, str(TOOLS))
    missing = [name for name in ("PyInstaller", "webview", "nicegui", "rdkit", "openmm", "pdbfixer", "meeko", "numpy", "scipy")
               if importlib.util.find_spec(name) is None]
    if missing:
        raise RuntimeError("Missing build/runtime packages: " + ", ".join(missing))
    for path in (ROOT / "desktop.py", ROOT / "desktop.spec", ROOT / "vina.exe", ROOT / "assets/workshop/egfr/project.db", ROOT / "assets/vendor/3Dmol-min.js", ROOT / "installer/assets/icon.ico"):
        if not path.is_file():
            raise RuntimeError(f"Required build input is missing: {path}")


def verify(environment: dict[str, str]) -> dict:
    report = WORK / "verification.json"
    if report.exists():
        report.unlink()
    # Exercise the distribution without borrowing Python or DLL search paths
    # from the development environment.
    environment = environment.copy()
    for name in ("PYTHONPATH", "PYTHONHOME", "CONDA_PREFIX", "CONDA_DEFAULT_ENV", "RDBASE", "OPENMM_PLUGIN_DIR"):
        environment.pop(name, None)
    windows = Path(environment.get("SystemRoot", r"C:\Windows"))
    environment["PATH"] = os.pathsep.join(str(path) for path in (windows / "System32", windows, windows / "System32/Wbem"))
    environment["INSILICO_PROJECTS"] = str(WORK / "verification-projects")
    environment["INSILICO_DATA_DIR"] = str(WORK / "verification-user")
    result = subprocess.run(
        [str(EXECUTABLE), "--verify-output", str(report)],
        cwd=str(PAYLOAD), env=environment, timeout=180, check=False,
    )
    if result.returncode:
        raise RuntimeError(f"Packaged self-test exited with code {result.returncode}; inspect {report}")
    if not report.is_file():
        raise RuntimeError(f"Packaged self-test did not write its report: {report}")
    data = json.loads(report.read_text(encoding="utf-8"))
    if data.get("can_start") is False or data.get("success") is False:
        raise RuntimeError(f"Packaged self-test reported a missing dependency: {report}")
    scientific_report = WORK / "scientific-verification.json"
    if scientific_report.exists():
        scientific_report.unlink()
    scientific_result = subprocess.run(
        [str(EXECUTABLE), "--verify-science", str(scientific_report)],
        cwd=str(PAYLOAD), env=environment, timeout=180, check=False)
    if scientific_result.returncode or not scientific_report.is_file():
        raise RuntimeError(f"Packaged scientific checks failed; inspect {scientific_report}")
    science = json.loads(scientific_report.read_text(encoding="utf-8"))
    if not science.get("passed"):
        raise RuntimeError(f"Packaged scientific checks did not pass: {scientific_report}")
    data["scientific_checks"] = science
    batch_report = WORK / "batch-verification.json"
    if batch_report.exists():
        batch_report.unlink()
    batch_result = subprocess.run([str(EXECUTABLE), "--verify-study", str(batch_report)],
        cwd=str(PAYLOAD), env=environment, timeout=240, check=False)
    if batch_result.returncode or not batch_report.is_file():
        raise RuntimeError(f"Packaged multi-target checks failed; inspect {batch_report}")
    batch = json.loads(batch_report.read_text(encoding="utf-8"))
    if not batch.get("passed"):
        raise RuntimeError(f"Packaged multi-target checks did not pass: {batch_report}")
    data["batch_checks"] = batch
    report.write_text(json.dumps(data, indent=2), encoding="utf-8")
    print(f"Packaged self-test passed: {report}", flush=True)
    return data


def write_manifest(verification: dict | None) -> Path:
    versions = {}
    for name in ("pyinstaller", "pywebview", "nicegui", "rdkit", "openmm", "pdbfixer", "meeko", "numpy", "scipy"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = "unknown"
    source_paths = [ROOT / name for name in ("desktop.py", "desktop.spec", "build_desktop.py", "app.py", "forms.py", "workshop_ui.py", "study_ui.py", "vina.exe")]
    for folder in ("core", "stages", "stations"):
        source_paths.extend(sorted((ROOT / folder).glob("*.py")))
    manifest = {
        "product": "In-Silico Docking Studio",
        "built_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "distribution": "Windows x64 directory bundle",
        "versions": versions,
        "verification": verification,
        "source_sha256": {str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path) for path in source_paths},
        "files": {},
    }
    destination = PAYLOAD / "build-manifest.json"
    for path in sorted(PAYLOAD.rglob("*")):
        if path.is_file() and path != destination:
            name = str(path.relative_to(PAYLOAD)).replace("\\", "/")
            manifest["files"][name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    destination.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return destination


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--zip", action="store_true", help="Also produce a ZIP containing the complete application folder.")
    parser.add_argument("--skip-verify", action="store_true", help="Skip execution of packaged diagnostics for an intermediate build.")
    parser.add_argument("--reuse", action="store_true", help="Verify and package an existing dist folder without rebuilding it.")
    arguments = parser.parse_args()
    try:
        check_prerequisites()
        WORK.mkdir(parents=True, exist_ok=True)
        environment = build_environment()
        if not arguments.reuse:
            command = [sys.executable, "-m", "PyInstaller", "--noconfirm", "--distpath", str(DIST), "--workpath", str(WORK), str(ROOT / "desktop.spec")]
            print("Building standalone Windows application...", flush=True)
            subprocess.run(command, cwd=str(ROOT), env=environment, check=True)
        if not EXECUTABLE.is_file():
            raise RuntimeError(f"The build did not produce {EXECUTABLE}")
        shutil.copy2(ROOT / "DESKTOP.md", PAYLOAD / "Getting Started.md")
        python_license = Path(sys.prefix) / "LICENSE.txt"
        if python_license.exists():
            licenses = PAYLOAD / "licenses"
            licenses.mkdir(exist_ok=True)
            shutil.copy2(python_license, licenses / "Python-LICENSE.txt")
        verification = None if arguments.skip_verify else verify(environment)
        manifest = write_manifest(verification)
        print(f"Executable: {EXECUTABLE}", flush=True)
        print(f"Build manifest: {manifest}", flush=True)
        if arguments.zip:
            archive = shutil.make_archive(str(DIST / "InSilicoTeachingSuite-Windows-x64"), "zip", root_dir=str(DIST), base_dir=PAYLOAD.name)
            print(f"Portable ZIP: {archive}", flush=True)
        print("Distribute the complete application folder or its ZIP; the EXE needs its _internal folder.", flush=True)
        return 0
    except (RuntimeError, subprocess.SubprocessError, OSError, ValueError) as error:
        print(f"Build failed: {error}", file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

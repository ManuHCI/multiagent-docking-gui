"""Build the separate Python 3.8 legacy edition and its per-user setup EXE."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
WORK = ROOT / "build/legacy"
PAYLOAD = ROOT / "dist/InSilicoDockingStudioLegacy"
EXE = PAYLOAD / "InSilicoDockingStudioLegacy.exe"
SETUP = ROOT / "dist/InSilicoDockingStudio-1.3.0-Legacy-Windows-x64-Setup.exe"
VERSION = "1.3.0-legacy.1"


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def check_prerequisites():
    if sys.platform != "win32" or sys.version_info[:3] != (3, 8, 10):
        raise RuntimeError("Build with the separate official Windows CPython 3.8.10 x64 prefix.")
    import struct
    if struct.calcsize("P") != 8:
        raise RuntimeError("The legacy scientific environment must be 64-bit.")
    for module in ("PyInstaller", "pefile", "tkinter", "rdkit", "openmm", "pdbfixer", "meeko", "numpy", "scipy", "pandas", "PIL"):
        if importlib.util.find_spec(module) is None:
            raise RuntimeError("Missing legacy build module: " + module)
    if importlib.metadata.version("pyinstaller") != "5.13.2":
        raise RuntimeError("Use the pinned PyInstaller 5.13.2 legacy bootloader.")
    for index, folder, expected in (("runtime-files.json", "files-Minimum", "14.29.30157.0"),
                                     ("ucrt-files.json", "ucrt-x64", "10.0.10240.16384")):
        base = ROOT / ".build-tools/legacy-vc-runtime"
        records = json.loads((base / index).read_text(encoding="utf-8-sig"))
        for record in records:
            if record["signature"] != "Valid" or record["version"] != expected:
                raise RuntimeError("The legacy Microsoft runtime provenance does not match the verified version.")
            if sha256(base / folder / record["name"]) != record["sha256"]:
                raise RuntimeError("A verified Microsoft runtime file changed: " + record["name"])
    for name in ("legacy_desktop.py", "legacy.spec", "legacy/ui.py", "legacy/compat.py", "vina.exe"):
        if not (ROOT / name).is_file():
            raise RuntimeError("Missing legacy build input: " + name)


def isolated_environment():
    environment = os.environ.copy()
    for name in ("PYTHONPATH", "PYTHONHOME", "CONDA_PREFIX", "CONDA_DEFAULT_ENV", "RDBASE", "OPENMM_PLUGIN_DIR"):
        environment.pop(name, None)
    windows = Path(environment.get("SystemRoot", "C:/Windows"))
    environment["PATH"] = os.pathsep.join(str(path) for path in (windows / "System32", windows, windows / "System32/Wbem"))
    environment["INSILICO_DATA_DIR"] = str(WORK / "verification-user")
    environment["INSILICO_PROJECTS"] = str(WORK / "verification-user/Projects")
    return environment


def verify():
    from legacy_pe_audit import audit_payload
    pe = audit_payload(PAYLOAD)
    (WORK / "pe-compatibility.json").write_text(json.dumps(pe, indent=2), encoding="utf-8")
    if not pe["passed"]:
        raise RuntimeError("Static legacy compatibility screening failed. Inspect build/legacy/pe-compatibility.json")
    reports = {"static_pe_screening": {key: value for key, value in pe.items() if key != "files"}}
    for name, flag, timeout in (("diagnostics", "--verify-output", 180), ("science", "--verify-science", 240),
                                 ("multi-target-docking", "--verify-study", 300), ("native-window", "--window-test-output", 180)):
        destination = WORK / (name + ".json")
        if destination.exists():
            destination.unlink()
        result = subprocess.run([str(EXE), flag, str(destination)], cwd=PAYLOAD,
                                env=isolated_environment(), timeout=timeout)
        if result.returncode or not destination.is_file():
            raise RuntimeError("Packaged legacy " + name + " failed. Inspect " + str(destination))
        record = json.loads(destination.read_text(encoding="utf-8"))
        if not record.get("passed", record.get("can_start", False)):
            raise RuntimeError("Packaged legacy " + name + " did not pass. Inspect " + str(destination))
        reports[name] = record
    return reports


def write_manifest(reports):
    source_files = [ROOT / name for name in ("legacy_desktop.py", "legacy.spec", "build_legacy.py", "legacy_pe_audit.py", "requirements-legacy.txt")]
    for folder in ("legacy", "core", "stages", "stations", "packaging/legacy"):
        source_files.extend(sorted((ROOT / folder).glob("*.py")))
    manifest = {"product": "In-Silico Docking Studio Legacy", "version": VERSION,
                "built_utc": datetime.now(timezone.utc).isoformat(), "python": sys.version,
                "intended_systems": ["Windows 7 SP1 x64", "Windows 8 x64", "Windows 8.1 x64"],
                "release_status": "Experimental legacy beta; older-Windows execution is unverified",
                "old_windows_execution_tested": False, "verification": reports,
                "versions": {name: importlib.metadata.version(name) for name in
                             ("pyinstaller", "rdkit-pypi", "OpenMM", "pdbfixer", "meeko", "numpy", "scipy", "pandas", "Pillow", "certifi")},
                "source_sha256": {str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path) for path in source_files if path.is_file()},
                "files": {}}
    for path in sorted(PAYLOAD.rglob("*")):
        if path.is_file() and path.name != "build-manifest.json":
            manifest["files"][str(path.relative_to(PAYLOAD)).replace("\\", "/")] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    destination = PAYLOAD / "build-manifest.json"
    destination.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def build_setup(manifest):
    if not manifest.get("verification"):
        raise RuntimeError("The public legacy installer requires completed packaged checks.")
    for relative, record in manifest["files"].items():
        if sha256(PAYLOAD / relative) != record["sha256"]:
            raise RuntimeError("The verified payload changed: " + relative)
    compiler = ROOT / ".build-tools/nsis-3.11/makensis.exe"
    if not compiler.is_file():
        raise RuntimeError("NSIS 3.11 is required in .build-tools/nsis-3.11.")
    includes = WORK / "uninstall-files.nsh"
    def quote(value):
        return value.replace("$", "$$").replace('"', '$\\"')
    paths = sorted(path for path in PAYLOAD.rglob("*") if path.is_file())
    lines = ['Delete "$INSTDIR\\' + quote(str(path.relative_to(PAYLOAD))) + '"' for path in paths]
    dirs = sorted((path for path in PAYLOAD.rglob("*") if path.is_dir()), key=lambda path: len(path.parts), reverse=True)
    lines.extend('RMDir "$INSTDIR\\' + quote(str(path.relative_to(PAYLOAD))) + '"' for path in dirs)
    includes.write_text("\n".join(lines) + "\n", encoding="utf-8-sig")
    subprocess.run([str(compiler), "/V3", "/INPUTCHARSET", "UTF8", "/DPAYLOAD=" + str(PAYLOAD),
                    "/DEXE_SHA256=" + sha256(EXE), "/DOUTPUT=" + str(SETUP), "/DUNINSTALL_FILES=" + str(includes),
                    str(ROOT / "installer/legacy.nsi")], cwd=ROOT / "installer", check=True)
    SETUP.with_suffix(".exe.sha256").write_text(sha256(SETUP) + "  " + SETUP.name + "\n", encoding="ascii")
    print("Legacy setup: " + str(SETUP), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reuse", action="store_true", help="Check an existing legacy payload instead of rebuilding.")
    parser.add_argument("--skip-verify", action="store_true", help="Intermediate developer build only; cannot create setup.")
    parser.add_argument("--setup", action="store_true", help="Create the separate legacy installer after all checks pass.")
    parser.add_argument("--zip", action="store_true", help="Create a complete portable ZIP.")
    args = parser.parse_args()
    try:
        check_prerequisites()
        WORK.mkdir(parents=True, exist_ok=True)
        if not args.reuse:
            environment = os.environ.copy()
            environment["PYTHONNOUSERSITE"] = "1"
            environment["PYINSTALLER_CONFIG_DIR"] = str(WORK / "cache")
            subprocess.run([sys.executable, "-m", "PyInstaller", "--noconfirm", "--distpath", str(ROOT / "dist"),
                            "--workpath", str(WORK / "pyinstaller"), str(ROOT / "legacy.spec")],
                           env=environment, cwd=ROOT, check=True)
        if not EXE.is_file():
            raise RuntimeError("Legacy build did not produce its executable.")
        shutil.copy2(ROOT / "LEGACY_WINDOWS.md", PAYLOAD / "Getting Started.txt")
        python_license = Path(sys.prefix) / "LICENSE.txt"
        if python_license.is_file():
            (PAYLOAD / "licenses").mkdir(exist_ok=True)
            shutil.copy2(python_license, PAYLOAD / "licenses/Python-3.8.10-LICENSE.txt")
        reports = {} if args.skip_verify else verify()
        manifest = write_manifest(reports)
        if args.zip:
            archive = shutil.make_archive(str(ROOT / "dist/InSilicoDockingStudio-1.3.0-Legacy-Windows-x64"), "zip",
                                          root_dir=str(PAYLOAD.parent), base_dir=PAYLOAD.name)
            Path(archive + ".sha256").write_text(sha256(Path(archive)) + "  " + Path(archive).name + "\n", encoding="ascii")
        if args.setup:
            build_setup(manifest)
        print("Legacy packaged checks completed. This build has not been executed on Windows 7/8/8.1.", flush=True)
        return 0
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print("Legacy build failed: " + str(error), file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

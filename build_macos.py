"""Build, verify and package native arm64 or Intel macOS app/DMG/ZIP releases.

Run on a Mac in environment-macos.yml; cross-compiling from Windows is not
supported. Each architecture is built and tested on its own native runner.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
APP_NAME = "In-Silico Docking Studio.app"
PACKAGES = ("nicegui", "pywebview", "pyinstaller", "rdkit", "openmm", "pdbfixer", "meeko", "gemmi", "numpy", "scipy", "vina")


def sha256(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def run(command: list[str], *, environment=None, timeout=300, cwd=None) -> None:
    subprocess.run(command, check=True, timeout=timeout, env=environment, cwd=cwd or ROOT)


def check_prerequisites(architecture: str) -> Path:
    if sys.platform != "darwin":
        raise RuntimeError("Build the macOS application on a Mac or the GitHub macOS runner, not Windows.")
    if architecture not in {"arm64", "x86_64"} or platform.machine() != architecture:
        raise RuntimeError(f"A native {architecture} Python is required; this Python reports {platform.machine()}.")
    missing = [name for name in ("PyInstaller", "webview", "rdkit", "openmm", "pdbfixer", "meeko", "nicegui", "PIL")
               if importlib.util.find_spec(name) is None]
    if missing:
        raise RuntimeError("Missing build dependencies: " + ", ".join(missing))
    # Meeko's top-level import also requires Gemmi. Check real imports before
    # freezing, rather than merely detecting an installed module's directory.
    from meeko import MoleculePreparation  # noqa: F401
    native_vina = Path(sys.prefix) / "bin/vina"
    if not native_vina.is_file() or not os.access(native_vina, os.X_OK):
        raise RuntimeError("Native conda-forge Vina is missing; create environment-macos.yml first.")
    description = subprocess.check_output(["/usr/bin/file", str(native_vina)], text=True)
    if "Mach-O" not in description or architecture not in description:
        raise RuntimeError(f"Vina must be a native {architecture} Mach-O executable: {description.strip()}")
    run([str(native_vina), "--version"])
    for name in ("desktop.py", "desktop_macos.spec", "assets/logo.png", "assets/vendor/3Dmol-min.js",
                 "assets/workshop/egfr/project.db", "assets/workshop/sv2a/project.db"):
        if not (ROOT / name).is_file():
            raise RuntimeError("Missing build input: " + name)
    return native_vina


def source_manifest(architecture: str, native_vina: Path) -> dict:
    from desktop import VERSION
    sources = [ROOT / name for name in ("desktop.py", "app.py", "forms.py", "study_ui.py", "workshop_ui.py",
                                       "desktop_macos.spec", "build_macos.py", "environment-macos.yml")]
    for folder in ("core", "stages", "stations"):
        sources.extend(sorted((ROOT / folder).glob("*.py")))
    sources.extend(sorted((ROOT / "packaging/macos").glob("*")))
    return {
        "product": "In-Silico Docking Studio", "version": VERSION,
        "architecture": architecture, "macos_minimum": "15.0", "built_macos": platform.mac_ver()[0],
        "python": sys.version, "built_utc": datetime.now(timezone.utc).isoformat(),
        "versions": {name: importlib.metadata.version(name) for name in PACKAGES},
        "source_sha256": {str(path.relative_to(ROOT)): sha256(path) for path in sources if path.is_file()},
        "native_vina_input_sha256": sha256(native_vina),
        "verification": "Packaged verification and final signed payload hashes are in the adjacent release-manifest.json.",
        "signing": "Developer ID" if os.environ.get("INSILICO_CODESIGN_IDENTITY") else "ad-hoc; not Apple notarized",
    }


def collect_licenses(destination: Path) -> None:
    """Retain shipped Python/native package license texts and provenance."""
    destination.mkdir(parents=True, exist_ok=True)
    notices = []
    for name in PACKAGES:
        distribution = importlib.metadata.distribution(name)
        notices.append(f"{name} {distribution.version}\nLicense: {distribution.metadata.get('License', 'See bundled license text')}\n")
        for relative in distribution.files or []:
            if any(token in Path(str(relative)).name.lower() for token in ("license", "copying", "copyright", "notice")):
                source = Path(distribution.locate_file(relative))
                if source.is_file() and source.stat().st_size < 3_000_000:
                    target = destination / name / str(relative).replace("..", "_")
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(source, target)
    # Conda's installed package records retain the cache source for native
    # licenses which are not necessarily copied into Python dist-info.
    for record_file in (Path(sys.prefix) / "conda-meta").glob("*.json"):
        record = json.loads(record_file.read_text())
        package = record.get("name", record_file.stem)
        license_root = Path(record.get("link", {}).get("source", "")) / "info/licenses"
        if license_root.is_dir():
            shutil.copytree(license_root, destination / "conda" / package, dirs_exist_ok=True)
    python_license = Path(sys.prefix) / "LICENSE.txt"
    if python_license.is_file():
        shutil.copy2(python_license, destination / "Python-LICENSE.txt")
    (destination / "THIRD_PARTY_NOTICES.txt").write_text("\n".join(notices), encoding="utf-8")


def make_icon(inputs: Path) -> None:
    from PIL import Image
    iconset = inputs / "icon.iconset"
    iconset.mkdir(parents=True, exist_ok=True)
    with Image.open(ROOT / "installer/assets/icon.ico") as original:
        original = original.convert("RGBA")
        for points in (16, 32, 128, 256, 512):
            for scale in (1, 2):
                size = points * scale
                image = original.resize((size, size), Image.Resampling.LANCZOS)
                image.save(iconset / f"icon_{points}x{points}{'@2x' if scale == 2 else ''}.png")
    run(["/usr/bin/iconutil", "-c", "icns", str(iconset), "-o", str(inputs / "icon.icns")])


def isolated_environment(work: Path) -> dict[str, str]:
    environment = os.environ.copy()
    for name in ("PYTHONPATH", "PYTHONHOME", "CONDA_PREFIX", "CONDA_DEFAULT_ENV", "RDBASE", "OPENMM_PLUGIN_DIR",
                 "DYLD_LIBRARY_PATH", "DYLD_FALLBACK_LIBRARY_PATH", "LD_LIBRARY_PATH", "QT_PLUGIN_PATH"):
        environment.pop(name, None)
    environment["PATH"] = "/usr/bin:/bin:/usr/sbin:/sbin"
    environment["INSILICO_DATA_DIR"] = str(work / "verification-user")
    environment["INSILICO_PROJECTS"] = str(work / "verification-user/Projects")
    return environment


def verify_app(app: Path, work: Path, *, native_window: bool) -> dict:
    executable = app / "Contents/MacOS/InSilicoTeachingSuite"
    environment = isolated_environment(work)
    results = {}
    checks = [("diagnostics", "--verify-output", "can_start", 240),
              ("science", "--verify-science", "passed", 300),
              ("multi-target-docking", "--verify-study", "passed", 360)]
    if native_window:
        checks.append(("native-window", "--window-test-output", "passed", 180))
    for name, argument, success, timeout in checks:
        report = work / f"{name}.json"
        if report.exists():
            report.unlink()
        run([str(executable), argument, str(report)], environment=environment, timeout=timeout, cwd=work)
        result = json.loads(report.read_text())
        if not result.get(success):
            raise RuntimeError(f"Packaged {name} check failed: {report}")
        if name == "native-window" and not all(result.get(key) for key in ("downloads_enabled", "server_stopped")):
            raise RuntimeError("The native app did not enable downloads or stop its server cleanly.")
        results[name] = result
    run(["/usr/bin/codesign", "--verify", "--deep", "--strict", str(app)])
    return results


def payload_manifest(app: Path) -> dict:
    files, links = {}, {}
    for path in sorted(app.rglob("*")):
        relative = str(path.relative_to(app))
        if path.is_symlink():
            links[relative] = os.readlink(path)
        elif path.is_file():
            files[relative] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    return {"files": files, "symlinks": links}


def make_archives(app: Path, distribution: Path, manifest: Path, version: str, architecture: str) -> list[Path]:
    stem = f"InSilicoDockingStudio-{version}-macOS-{architecture}"
    zip_file, dmg = distribution / f"{stem}.zip", distribution / f"{stem}.dmg"
    # ditto preserves macOS executable modes, symbolic links and app metadata.
    # The ZIP contains the app and the release evidence alongside it.
    with tempfile.TemporaryDirectory(prefix="dmg-inputs-", dir=distribution) as temporary:
        stage = Path(temporary)
        run(["/usr/bin/ditto", str(app), str(stage / APP_NAME)], timeout=300)
        shutil.copy2(ROOT / "packaging/macos/READ_ME_FIRST.txt", stage / "READ_ME_FIRST.txt")
        shutil.copy2(manifest, stage / "release-manifest.json")
        run(["/usr/bin/ditto", "-c", "-k", "--sequesterRsrc", str(stage), str(zip_file)], timeout=600)
        (stage / "Applications").symlink_to("/Applications", target_is_directory=True)
        if dmg.exists():
            dmg.unlink()
        run(["/usr/bin/hdiutil", "create", "-volname", "In-Silico Docking Studio", "-srcfolder", str(stage),
             "-ov", "-format", "UDZO", str(dmg)], timeout=600)
    for archive in (zip_file, dmg):
        archive.with_suffix(archive.suffix + ".sha256").write_text(f"{sha256(archive)}  {archive.name}\n", encoding="ascii")
    return [dmg, zip_file]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arch", choices=("arm64", "x86_64"), default=platform.machine())
    parser.add_argument("--skip-native-window", action="store_true", help="Skip only the GUI smoke check on a Mac without a desktop session.")
    arguments = parser.parse_args()
    try:
        vina = check_prerequisites(arguments.arch)
        work = ROOT / "build" / "macos" / arguments.arch
        distribution = ROOT / "dist" / "macos" / arguments.arch
        inputs = work / "build-inputs"
        inputs.mkdir(parents=True, exist_ok=True)
        distribution.mkdir(parents=True, exist_ok=True)
        metadata = source_manifest(arguments.arch, vina)
        (inputs / "build-manifest.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        collect_licenses(inputs / "licenses")
        make_icon(inputs)
        environment = os.environ.copy()
        environment["INSILICO_MAC_BUILD_INPUTS"] = str(inputs)
        environment["PYINSTALLER_CONFIG_DIR"] = str(work / "cache")
        environment["PYTHONNOUSERSITE"] = "1"
        run([sys.executable, "-m", "PyInstaller", "--noconfirm", "--distpath", str(distribution),
             "--workpath", str(work / "pyinstaller"), str(ROOT / "desktop_macos.spec")],
            environment=environment, timeout=1800)
        app = distribution / APP_NAME
        if not app.is_dir():
            raise RuntimeError("PyInstaller did not create the macOS application.")
        results = verify_app(app, work, native_window=not arguments.skip_native_window)
        # Keep the signed .app immutable. Its internal source manifest was
        # bundled before signing; final signed payload hashes live alongside it.
        manifest = distribution / "release-manifest.json"
        manifest.write_text(json.dumps({**metadata, "verification": results,
            "native_window_verified": not arguments.skip_native_window,
            "notarized": False, **payload_manifest(app)}, indent=2), encoding="utf-8")
        outputs = make_archives(app, distribution, manifest, metadata["version"], arguments.arch)
        for output in outputs:
            print(f"Built {output} ({output.stat().st_size / 1e6:.1f} MB)", flush=True)
        print(f"Verified release manifest: {manifest}", flush=True)
        return 0
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print(f"macOS build failed: {error}", file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())


"""Export a public source tree without user projects, settings or build caches.

Only allowlisted source, supporting assets and build inputs are exported.
Bundled teaching-case database paths are made portable in the export copy;
the working examples and saved projects in this checkout are never changed.
"""
from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import shutil
import sqlite3
import zipfile

ROOT_FILES = {
    ".gitignore", "CITATION.cff", "README.md", "LICENSE", "DESKTOP.md", "INSTALL.md", "MACOS.md",
    "STUDIO_VERIFICATION.md", "app.py", "forms.py", "study_ui.py", "workshop_ui.py",
    "desktop.py", "desktop.spec", "desktop_macos.spec", "build_desktop.py",
    "build_setup.py", "build_macos.py", "export_source.py", "environment.yml",
    "environment-macos.yml", "run.bat", "run.command", "vina.exe",
    "verify_desktop_ui.py", "verify_studio_ui.py", "verify_setup.py",
    "verify_upgrade.py", "verify_live_docking.py",
    "legacy_desktop.py", "legacy.spec", "build_legacy.py", "legacy_pe_audit.py",
    "requirements-legacy.txt", "LEGACY_WINDOWS.md", "LEGACY_VERIFICATION.md", "verify_legacy_ui.py", "verify_legacy_setup.py",
}
DIRECTORY_PATTERNS = {
    "legacy": {".py"}, "core": {".py"}, "stages": {".py"}, "stations": {".py"}, "tests": {".py"},
    "assets": {".png", ".js", ".txt", ".md", ".json", ".db", ".pdb", ".sdf", ".pdbqt"},
    "installer": {".nsi", ".ico", ".png"},
    "packaging": {".py", ".sh", ".plist", ".txt", ".md", ".png", ".icns", ".json"},
    ".github/workflows": {".yml", ".yaml"},
}
EXCLUDED_COMPONENTS = {"__pycache__", ".git", ".pytest_cache", ".DS_Store"}


def source_files(root: Path) -> list[Path]:
    """Use an allowlist so a later local file cannot silently become public."""
    files = {root / name for name in ROOT_FILES if (root / name).is_file() and not (root / name).is_symlink()}
    for directory, suffixes in DIRECTORY_PATTERNS.items():
        folder = root / directory
        if not folder.is_dir():
            continue
        for path in folder.rglob("*"):
            if not path.is_file() or path.is_symlink():
                continue
            if any(part in EXCLUDED_COMPONENTS for part in path.relative_to(root).parts):
                continue
            if path.suffix.lower() in suffixes:
                files.add(path)
    return sorted(files, key=lambda path: path.relative_to(root).as_posix())


def _portable(value, previous: str, replacement: str):
    if isinstance(value, str):
        normalized = value.replace("\\", "/")
        origin = previous.replace("\\", "/").rstrip("/")
        if normalized.casefold() == origin.casefold():
            return replacement
        if normalized.casefold().startswith(origin.casefold() + "/"):
            return replacement + normalized[len(origin):]
        return value
    if isinstance(value, dict):
        return {key: _portable(item, previous, replacement) for key, item in value.items()}
    if isinstance(value, list):
        return [_portable(item, previous, replacement) for item in value]
    return value


def sanitize_teaching_case(folder: Path) -> None:
    """Rewrite only exported case metadata, then remove stale SQLite pages."""
    summary_path = folder / "project.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    previous = summary.get("location", "")
    if not previous:
        raise ValueError(f"Teaching case has no recorded location: {folder.name}")
    replacement = "/bundled-examples/" + folder.name
    summary = _portable(summary, previous, replacement)
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    with sqlite3.connect(folder / "project.db") as connection:
        for column in ("values_json", "control_json"):
            for stage_id, text in connection.execute(f"SELECT stage_id,{column} FROM steps").fetchall():
                if not text:
                    continue
                original = json.loads(text)
                portable = _portable(original, previous, replacement)
                if portable != original:
                    connection.execute(f"UPDATE steps SET {column}=? WHERE stage_id=?",
                                       (json.dumps(portable, ensure_ascii=False), stage_id))
        for entry_id, text in connection.execute("SELECT id,parameters FROM provenance").fetchall():
            if not text:
                continue
            try:
                original = json.loads(text)
            except ValueError:
                continue
            portable = _portable(original, previous, replacement)
            if portable != original:
                connection.execute("UPDATE provenance SET parameters=? WHERE id=?",
                                   (json.dumps(portable, ensure_ascii=False), entry_id))
        connection.commit()
        connection.execute("VACUUM")
    # Bytes from the original personal filesystem must not survive in free pages.
    raw = (folder / "project.db").read_bytes()
    if previous.encode() in raw or previous.replace("\\", "\\\\").encode() in raw:
        raise ValueError(f"An original teaching-case path remains in {folder.name}")


def export_source(root: Path, destination: Path, archive: bool = True) -> dict:
    root, destination = root.resolve(), destination.resolve()
    if destination == root or not destination.is_relative_to(root):
        raise ValueError("The export destination must be a new subfolder of the source workspace.")
    if destination.relative_to(root).parts[0] in {name.split("/")[0] for name in DIRECTORY_PATTERNS}:
        raise ValueError("The export destination cannot be inside a source or asset directory.")
    if destination.exists():
        raise FileExistsError("The export destination already exists; select a new folder.")
    selected = source_files(root)
    destination.mkdir(parents=True)
    for source in selected:
        relative = source.relative_to(root)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    for case in (destination / "assets" / "workshop").glob("*/project.db"):
        sanitize_teaching_case(case.parent)
    entries = []
    for source in selected:
        relative = source.relative_to(root)
        target = destination / relative
        entries.append({"path": relative.as_posix(), "bytes": target.stat().st_size,
                        "sha256": hashlib.sha256(target.read_bytes()).hexdigest()})
    manifest = {"product": "In-Silico Docking Studio", "format": "public source export",
                "excluded": ["user projects", "settings.json", "backups", "build outputs",
                             "verification/temporary data", "local environments", "internal handoff notes"],
                "teaching_cases": "Only bundled public-data examples; exported paths are portable.",
                "files": entries}
    (destination / "source-publication-manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if archive:
        archive_path = destination.with_suffix(".zip")
        if archive_path.exists():
            raise FileExistsError("The source ZIP already exists; select a new export folder.")
        with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as bundle:
            for path in sorted(destination.rglob("*")):
                if not path.is_file():
                    continue
                name = path.relative_to(destination).as_posix()
                info = zipfile.ZipInfo(name, date_time=(2026, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                mode = 0o755 if path.suffix in {".sh", ".command"} else 0o644
                info.external_attr = (0o100000 | mode) << 16
                bundle.writestr(info, path.read_bytes())
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None,
                        help="New export folder under this workspace; never overwrites an existing folder.")
    parser.add_argument("--no-zip", action="store_true", help="Export the staging tree without a ZIP.")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    destination = args.output or root / "dist" / ("source-publication-" + datetime.now().strftime("%Y%m%d-%H%M%S"))
    manifest = export_source(root, destination, archive=not args.no_zip)
    print(json.dumps({"folder": str(destination.resolve()), "files": len(manifest["files"]),
                      "bytes": sum(item["bytes"] for item in manifest["files"]),
                      "zip": None if args.no_zip else str(destination.resolve().with_suffix(".zip"))}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

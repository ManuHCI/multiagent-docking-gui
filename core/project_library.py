"""Read-only project summaries and bounded, explicit project deletion.

The library is scanned only when the Projects section is opened. Merely saving
projects does not load their structures or results into application memory.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import sqlite3
import stat
from typing import Any


def human_size(size: int) -> str:
    """Format file bytes for the project list."""
    value = float(max(0, size))
    for suffix in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024 or suffix == "TB":
            return f"{int(value)} B" if suffix == "B" else f"{value:.1f} {suffix}"
        value /= 1024
    return "0 B"


def _linked(path: Path) -> bool:
    info = path.lstat()
    return stat.S_ISLNK(info.st_mode) or bool(
        getattr(info, "st_file_attributes", 0) & getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0x400)
    )


def _kind(folder: Path) -> str | None:
    for marker, kind in (("study.json", "study"), ("project.db", "legacy")):
        path = folder / marker
        if path.is_file() and not _linked(path):
            return kind
    return None


def _tree_size(folder: Path) -> int:
    size = 0
    pending = [folder]
    while pending:
        current = pending.pop()
        try:
            with os.scandir(current) as entries:
                for entry in entries:
                    try:
                        path = Path(entry.path)
                        if _linked(path):
                            continue
                        if entry.is_dir(follow_symlinks=False):
                            pending.append(path)
                        elif entry.is_file(follow_symlinks=False):
                            size += entry.stat(follow_symlinks=False).st_size
                    except OSError:
                        continue
        except OSError:
            continue
    return size


def _read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError):
        return {}


def _legacy_summary(folder: Path) -> dict[str, Any]:
    summary = _read_json(folder / "project.json")
    result = {"name": summary.get("name") or folder.name, "status": "saved",
              "target_count": 1, "drug_count": 0}
    connection = None
    try:
        # Do not call Project.open: relocation and dependency reconciliation
        # deliberately write to a project and are unnecessary for this list.
        connection = sqlite3.connect((folder / "project.db").as_uri() + "?mode=ro", uri=True)
        count = connection.execute("SELECT COUNT(*) FROM steps WHERE status='complete'").fetchone()[0]
        result["status"] = f"{count} steps saved" if count else "saved"
        row = connection.execute("SELECT values_json FROM steps WHERE stage_id='05_ligands'").fetchone()
        if row:
            data = json.loads(row[0] or "{}")
            result["drug_count"] = len(data.get("prepared") or [])
    except (OSError, sqlite3.Error, ValueError, TypeError):
        pass
    finally:
        if connection is not None:
            connection.close()
    return result


def list_projects(root: Path | str) -> list[dict[str, Any]]:
    """Summarize direct-child study and earlier project folders without writes."""
    root = Path(root).resolve()
    if not root.is_dir():
        return []
    projects = []
    for folder in root.iterdir():
        try:
            if folder.name.startswith(".") or _linked(folder) or not folder.is_dir():
                continue
            kind = _kind(folder)
            if kind is None:
                continue
            marker = folder / ("study.json" if kind == "study" else "project.db")
            updated = datetime.fromtimestamp(marker.stat().st_mtime, timezone.utc).isoformat(timespec="seconds")
            if kind == "study":
                data = _read_json(marker)
                details = {"name": data.get("name") or folder.name,
                           "status": data.get("status", "saved") if data else "unavailable",
                           "updated": data.get("updated") or updated,
                           "target_count": len(data.get("targets") or []),
                           "drug_count": len(data.get("drugs") or [])}
            else:
                details = {**_legacy_summary(folder), "updated": updated}
            projects.append({**details, "folder": str(folder), "kind": kind,
                             "size_bytes": _tree_size(folder)})
        except (OSError, ValueError, TypeError):
            continue
    return sorted(projects, key=lambda item: str(item["updated"]), reverse=True)


def _checked_folder(root: Path | str, folder: Path | str) -> tuple[Path, str]:
    supplied_root = Path(root).absolute()
    if _linked(supplied_root):
        raise ValueError("Projects cannot be deleted through a linked storage folder.")
    root_path = supplied_root.resolve(strict=True)
    candidate = Path(folder).absolute()
    if candidate.name.startswith(".") or candidate.parent.resolve() != root_path:
        raise ValueError("Only a project directly inside the project storage folder can be deleted.")
    if _linked(candidate):
        raise ValueError("A linked folder cannot be deleted as a project.")
    resolved = candidate.resolve(strict=True)
    if resolved.parent != root_path or not resolved.is_dir():
        raise ValueError("The selected project folder is outside project storage.")
    kind = _kind(resolved)
    if kind is None:
        raise ValueError("The selected folder does not contain a saved project.")
    return resolved, kind


def _reject_linked_contents(folder: Path) -> None:
    # Refuse a tree with junctions, symbolic links or other reparse points.
    # This also avoids silently deleting shared or externally stored inputs.
    pending = [folder]
    while pending:
        current = pending.pop()
        with os.scandir(current) as entries:
            for entry in entries:
                path = Path(entry.path)
                if _linked(path):
                    raise ValueError("This project contains a linked file or folder. Remove the link before deleting the project.")
                if entry.is_dir(follow_symlinks=False):
                    pending.append(path)


def project_lock(folder: Path | str):
    """The shared calculation/deletion lock for any saved project folder."""
    from .study import Study, _lock
    return _lock(Study(Path(folder), {}))


def delete_project(root: Path | str, folder: Path | str) -> None:
    """Permanently remove one explicitly selected saved project, never siblings.

    The caller presents the project name and deletion confirmation. The same
    in-process lock used by scientific study stages prevents deleting a running
    study, including a calculation started by another UI client.
    """
    candidate, _ = _checked_folder(root, folder)
    lock = project_lock(candidate)
    if not lock.acquire(blocking=False):
        raise ValueError("This project is calculating. Wait for the calculation to finish before deleting it.")
    try:
        # Check again after acquiring the lock, immediately before removal.
        candidate, _ = _checked_folder(root, candidate)
        _reject_linked_contents(candidate)
        shutil.rmtree(candidate)
    finally:
        lock.release()

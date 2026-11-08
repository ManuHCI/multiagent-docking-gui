"""
Project container.

A project is one folder on disk plus a small SQLite index. The whole
folder can be zipped and emailed to a supervisor and will open intact
on another machine -- specification Section 5.3.

Folder layout:

    MyProject/
        project.json          human-readable summary, safe to open in Notepad
        project.db            SQLite index of steps and provenance
        steps/
            01_target/
            02_structure/
            ...               each stage owns its own folder of artefacts
        exports/              MD job packages, reports
        incoming/             files the user brings back from a website
"""

from __future__ import annotations

import json
import sqlite3
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .contract import ProvenanceEntry, StageResult

SCHEMA_VERSION = 1

# Dependencies describe the science, not the order of the sidebar. ADMET
# needs prepared molecules, but does not depend on docking or simulation.
DEPENDENCIES = {
    "01_target": (), "02_structure": ("01_target",),
    "03_receptor": ("02_structure",), "04_site": ("03_receptor",),
    "05_ligands": (),
    "06_dock": ("03_receptor", "04_site", "05_ligands"),
    "07_pose": ("06_dock",), "08_admet": ("05_ligands",),
    "08_admet_swissadme": ("05_ligands",),
    "09_md": ("06_dock",), "09_md_export": ("06_dock",),
    "10_report": ("06_dock", "08_admet", "09_md"),
}


def downstream(stage_id: str) -> set[str]:
    found: set[str] = set()
    pending = [stage_id]
    while pending:
        parent = pending.pop()
        for child, inputs in DEPENDENCIES.items():
            if parent in inputs and child not in found:
                found.add(child)
                pending.append(child)
    return found

STEP_FOLDERS = [
    "01_target",
    "02_structure",
    "03_receptor",
    "04_site",
    "05_ligands",
    "06_dock",
    "07_pose",
    "08_admet",
    "09_md_export",
    "10_report",
]


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


class Project:
    """A single piece of work, saved and resumable."""

    def __init__(self, folder: Path) -> None:
        self.folder = Path(folder).resolve()
        self.db_path = self.folder / "project.db"
        self._conn: sqlite3.Connection | None = None

    # -- creation and opening -------------------------------------------

    @classmethod
    def create(cls, folder: Path, name: str, note: str = "") -> "Project":
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        project = cls(folder)

        for sub in ("steps", "exports", "incoming"):
            (folder / sub).mkdir(exist_ok=True)
        for step in STEP_FOLDERS:
            (folder / "steps" / step).mkdir(exist_ok=True)

        project._init_db()
        project._write_summary({
            "name": name,
            "note": note,
            "created": _now(),
            "schema_version": SCHEMA_VERSION,
            "location": str(folder.resolve()),
        })
        return project

    @classmethod
    def open(cls, folder: Path) -> "Project":
        folder = Path(folder)
        if not (folder / "project.db").exists():
            raise FileNotFoundError(
                f"There is no project in {folder}. "
                "Choose a folder that contains project.db."
            )
        project = cls(folder)
        project.relocate_paths()
        project.reconcile_dependencies()
        return project

    # -- database --------------------------------------------------------

    @property
    def conn(self) -> sqlite3.Connection:
        if self._conn is None:
            # check_same_thread=False because the interface runs each
            # stage in a worker thread to keep the window responsive,
            # while the connection is made on the main thread. Safe
            # here: stages run one at a time, never concurrently.
            self._conn = sqlite3.connect(self.db_path,
                                         check_same_thread=False)
            self._conn.row_factory = sqlite3.Row
        return self._conn

    def _init_db(self) -> None:
        self.conn.executescript(
            """
            CREATE TABLE IF NOT EXISTS steps (
                stage_id   TEXT PRIMARY KEY,
                title      TEXT NOT NULL,
                status     TEXT NOT NULL,
                summary    TEXT,
                values_json TEXT,
                control_json TEXT,
                updated    TEXT NOT NULL
            );

            CREATE TABLE IF NOT EXISTS provenance (
                id          INTEGER PRIMARY KEY AUTOINCREMENT,
                stage_id    TEXT NOT NULL,
                tool        TEXT NOT NULL,
                version     TEXT NOT NULL,
                parameters  TEXT,
                seed        INTEGER,
                source_release TEXT,
                citation    TEXT,
                timestamp   TEXT NOT NULL
            );
            """
        )
        self.conn.commit()

    # -- summary file ----------------------------------------------------

    def _write_summary(self, data: dict[str, Any]) -> None:
        (self.folder / "project.json").write_text(
            json.dumps(data, indent=2), encoding="utf-8"
        )

    def relocate_paths(self) -> None:
        """Rebase file references when a complete project folder is copied."""
        summary = self.summary
        previous = summary.get("location")
        current = str(self.folder.resolve())
        def rewrite(value):
            if isinstance(value, str):
                if previous and (value.lower().startswith(previous.lower() + "\\")
                                 or value.startswith(previous + "/")):
                    # Stored examples can originate on another operating
                    # system. Joining native path components converts both
                    # Windows backslashes and POSIX separators on relocation.
                    parts = value[len(previous):].replace("\\", "/").lstrip("/").split("/")
                    if not any(part == ".." for part in parts):
                        return str(self.folder.joinpath(*parts))
                # Older command-line runs sometimes recorded cwd-relative
                # output paths. Normalise only paths inside this project.
                try:
                    path = Path(value)
                    if not path.is_absolute() and len(path.parts) > 1:
                        resolved = path.resolve()
                        if resolved.is_relative_to(self.folder):
                            return str(resolved)
                except (OSError, ValueError):
                    pass
                return value
            if isinstance(value, dict):
                return {k: rewrite(v) for k, v in value.items()}
            if isinstance(value, list):
                return [rewrite(v) for v in value]
            return value
        for row in self.conn.execute("SELECT stage_id,values_json FROM steps").fetchall():
            original = json.loads(row["values_json"] or "{}")
            value = rewrite(original)
            if value != original:
                self.conn.execute("UPDATE steps SET values_json=? WHERE stage_id=?",
                                  (json.dumps(value), row["stage_id"]))
        self.conn.commit()
        if previous != current:
            summary["location"] = current
            self._write_summary(summary)

    @property
    def summary(self) -> dict[str, Any]:
        path = self.folder / "project.json"
        if not path.exists():
            return {}
        return json.loads(path.read_text(encoding="utf-8"))

    @property
    def name(self) -> str:
        return self.summary.get("name", self.folder.name)

    # -- step state ------------------------------------------------------

    def step_folder(self, stage_id: str) -> Path:
        path = self.folder / "steps" / stage_id
        path.mkdir(parents=True, exist_ok=True)
        return path

    def save_step(self, stage_id: str, title: str, result: StageResult) -> None:
        control = asdict(result.control) if result.control else None
        values_json = json.dumps(result.values, default=str)
        control_json = json.dumps(control, default=str) if control else None
        old = self.get_step(stage_id)
        changed = (old is None or old["values_json"] != values_json
                   or old["control_json"] != control_json
                   or old["status"] != ("complete" if result.success else "failed"))
        # An unchanged rerun does not make its consumers older than its inputs.
        updated = _now() if changed else old["updated"]
        self.conn.execute(
            """
            INSERT INTO steps
                (stage_id, title, status, summary, values_json, control_json, updated)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(stage_id) DO UPDATE SET
                title=excluded.title,
                status=excluded.status,
                summary=excluded.summary,
                values_json=excluded.values_json,
                control_json=excluded.control_json,
                updated=excluded.updated
            """,
            (
                stage_id,
                title,
                "complete" if result.success else "failed",
                result.summary,
                values_json,
                control_json,
                updated,
            ),
        )
        if changed:
            for child in downstream(stage_id):
                self.conn.execute("UPDATE steps SET status='stale' WHERE stage_id=?",
                                  (child,))
        self.conn.commit()

    def reconcile_dependencies(self) -> list[str]:
        """Conservatively flag results saved before their current inputs.

        This also repairs older project indexes that predate invalidation.
        Outputs are retained for provenance, never silently deleted.
        """
        rows = {r["stage_id"]: dict(r) for r in self.conn.execute("SELECT * FROM steps")}
        stale: list[str] = []
        for _ in range(len(DEPENDENCIES)):
            changed = False
            for child, parents in DEPENDENCIES.items():
                row = rows.get(child)
                if not row or row["status"] != "complete":
                    continue
                if any(p not in rows or rows[p]["status"] != "complete"
                       or rows[p]["updated"] > row["updated"] for p in parents):
                    row["status"] = "stale"
                    stale.append(child)
                    self.conn.execute("UPDATE steps SET status='stale' WHERE stage_id=?", (child,))
                    changed = True
            if not changed:
                break
        self.conn.commit()
        return stale

    def display_status(self, stage_id: str) -> tuple[str, str, str]:
        """Execution and scientific verdicts must not share a green tick."""
        step = self.get_step(stage_id)
        if step is None:
            return "radio_button_unchecked", "#9aa0a6", "Not run"
        if step["status"] == "stale":
            return "update", "#b26a00", "Needs updating"
        if step["status"] != "complete":
            return "error", "#b3261e", "Run failed"
        control = json.loads(step["control_json"] or "{}")
        return {
            "pass": ("check_circle", "#2e7d32", "Finished · check passed"),
            "warn": ("warning", "#b26a00", "Finished · review warning"),
            "fail": ("cancel", "#b3261e", "Finished · check failed"),
        }.get(control.get("verdict"), ("help_outline", "#5f6368", "Finished · not checked"))

    def get_step(self, stage_id: str) -> dict[str, Any] | None:
        row = self.conn.execute(
            "SELECT * FROM steps WHERE stage_id = ?", (stage_id,)
        ).fetchone()
        return dict(row) if row else None

    def is_complete(self, stage_id: str) -> bool:
        step = self.get_step(stage_id)
        return bool(step and step["status"] == "complete")

    def completed_steps(self) -> list[str]:
        rows = self.conn.execute(
            "SELECT stage_id FROM steps WHERE status = 'complete' ORDER BY stage_id"
        ).fetchall()
        return [r["stage_id"] for r in rows]

    # -- provenance ------------------------------------------------------

    def add_provenance(self, entry: ProvenanceEntry) -> None:
        self.conn.execute(
            """
            INSERT INTO provenance
                (stage_id, tool, version, parameters, seed,
                 source_release, citation, timestamp)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                entry.stage,
                entry.tool,
                entry.version,
                json.dumps(entry.parameters, default=str),
                entry.seed,
                entry.source_release,
                entry.citation,
                entry.timestamp,
            ),
        )
        self.conn.commit()

    def provenance(self) -> list[dict[str, Any]]:
        rows = self.conn.execute(
            "SELECT * FROM provenance ORDER BY id"
        ).fetchall()
        out = []
        for r in rows:
            d = dict(r)
            d["parameters"] = json.loads(d["parameters"] or "{}")
            out.append(d)
        return out

    def close(self) -> None:
        if self._conn is not None:
            self._conn.close()
            self._conn = None

"""Deleting one project must never delete other work or follow links."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from core.project import Project
from core.project_library import delete_project, human_size, list_projects
from core.study import Study, _lock


class ProjectLibraryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="project_library_")
        self.area = Path(self.temporary.name)
        self.root = self.area / "Projects"
        self.root.mkdir()
        self.study = self.root / "Study one"
        self.study.mkdir()
        (self.study / "study.json").write_text(json.dumps({
            "name": "First study", "version": 1, "updated": "2026-01-01T12:00:00+00:00",
            "status": "completed", "targets": [{}, {}], "drugs": [{}, {}, {}],
        }), encoding="utf-8")
        (self.study / "result.pdbqt").write_bytes(b"saved scientific result")
        self.legacy = self.root / "Older project"
        project = Project.create(self.legacy, "Earlier work")
        project.close()
        self.outside = self.area / "Outside"
        self.outside.mkdir()
        (self.outside / "project.db").write_bytes(b"outside marker")
        (self.outside / "keep.txt").write_bytes(b"never delete this")

    def tearDown(self):
        self.temporary.cleanup()

    def test_list_is_read_only_includes_both_formats_and_correct_counts(self):
        db = self.legacy / "project.db"
        before = (hashlib.sha256(db.read_bytes()).hexdigest(), db.stat().st_mtime_ns)
        rows = list_projects(self.root)
        self.assertEqual({row["kind"] for row in rows}, {"study", "legacy"})
        study = next(row for row in rows if row["kind"] == "study")
        self.assertEqual((study["name"], study["target_count"], study["drug_count"]), ("First study", 2, 3))
        self.assertEqual(study["size_bytes"], sum(p.stat().st_size for p in self.study.iterdir()))
        self.assertEqual(before, (hashlib.sha256(db.read_bytes()).hexdigest(), db.stat().st_mtime_ns))

    def test_delete_study_preserves_siblings_and_outside(self):
        delete_project(self.root, self.study)
        self.assertFalse(self.study.exists())
        self.assertTrue((self.legacy / "project.db").is_file())
        self.assertEqual((self.outside / "keep.txt").read_bytes(), b"never delete this")

    def test_delete_legacy_closes_read_only_list_handle(self):
        list_projects(self.root)
        delete_project(self.root, self.legacy)
        self.assertFalse(self.legacy.exists())
        self.assertTrue((self.study / "study.json").is_file())

    def test_rejects_storage_root_outside_and_unrecognized_folders(self):
        ordinary = self.root / "Unrelated documents"
        ordinary.mkdir()
        (ordinary / "notes.txt").write_text("keep", encoding="utf-8")
        for path in (self.root, self.outside, ordinary, self.root / ".." / "Outside"):
            with self.subTest(path=path), self.assertRaises(ValueError):
                delete_project(self.root, path)
            self.assertTrue(path.exists())

    def test_running_study_cannot_be_deleted_and_lock_remains_held(self):
        lock = _lock(Study(self.study, {}))
        lock.acquire()
        try:
            with self.assertRaisesRegex(ValueError, "calculating"):
                delete_project(self.root, self.study)
            self.assertTrue(lock.locked())
            self.assertTrue(self.study.exists())
        finally:
            lock.release()
        delete_project(self.root, self.study)
        self.assertFalse(self.study.exists())

    def test_refuses_linked_project_or_linked_contents(self):
        linked = self.root / "Linked project"
        def create_link(path):
            try:
                path.symlink_to(self.outside, target_is_directory=True)
            except (OSError, NotImplementedError):
                if os.name != "nt":
                    self.skipTest("Creating symbolic links is unavailable")
                # NTFS junctions do not require Windows developer mode or
                # elevated privileges and exercise the actual reparse guard.
                environment = dict(os.environ, INSILICO_TEST_LINK=str(path),
                                   INSILICO_TEST_TARGET=str(self.outside))
                shell = str(Path(os.environ.get("SystemRoot", "C:/Windows")) / "System32/WindowsPowerShell/v1.0/powershell.exe")
                subprocess.run([shell, "-NoProfile", "-Command",
                                "$ErrorActionPreference='Stop'; New-Item -ItemType Junction -Path $env:INSILICO_TEST_LINK -Target $env:INSILICO_TEST_TARGET | Out-Null"],
                               env=environment, check=True, capture_output=True, timeout=20)
        create_link(linked)
        self.assertNotIn(str(linked), {row["folder"] for row in list_projects(self.root)})
        with self.assertRaisesRegex(ValueError, "linked"):
            delete_project(self.root, linked)
        nested = self.study / "external"
        create_link(nested)
        with self.assertRaisesRegex(ValueError, "linked"):
            delete_project(self.root, self.study)
        self.assertEqual((self.outside / "keep.txt").read_bytes(), b"never delete this")

    def test_corrupt_saved_study_is_listed_and_can_be_explicitly_deleted(self):
        (self.study / "study.json").write_text("{incomplete", encoding="utf-8")
        row = next(row for row in list_projects(self.root) if row["kind"] == "study")
        self.assertEqual(row["status"], "unavailable")
        delete_project(self.root, self.study)
        self.assertTrue(self.legacy.exists())

    def test_size_labels(self):
        self.assertEqual(human_size(0), "0 B")
        self.assertEqual(human_size(1024), "1.0 KB")
        self.assertEqual(human_size(3 * 1024 ** 3), "3.0 GB")


if __name__ == "__main__":
    unittest.main()

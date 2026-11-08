"""Compatibility checks for sharing the scientific backend with Python 3.8."""
from __future__ import annotations

import ast
from pathlib import Path, PurePosixPath, PureWindowsPath
import unittest

from legacy.compat import _is_relative_to, install_compatibility


class LegacyCompatibilityTests(unittest.TestCase):
    def test_containment_backport_handles_siblings_and_multiple_components(self):
        path = PurePosixPath("/projects/Study/poses/drug.pdbqt")
        self.assertTrue(_is_relative_to(path, "/projects", "Study"))
        self.assertTrue(_is_relative_to(path, path))
        self.assertFalse(_is_relative_to(path, "/projects/Study sibling"))
        self.assertFalse(_is_relative_to(path, "/other"))
        self.assertFalse(_is_relative_to(path, "relative/path"))
        windows = PureWindowsPath(r"C:\Projects\Study\poses\drug.pdbqt")
        self.assertTrue(_is_relative_to(windows, r"c:\projects\study"))
        self.assertFalse(_is_relative_to(windows, r"D:\Projects\Study"))

    def test_install_keeps_the_native_method_on_modern_python(self):
        original = getattr(Path, "is_relative_to", None)
        install_compatibility()
        self.assertTrue(Path("a/b").is_relative_to(Path("a")))
        if original is not None:
            self.assertIs(Path.is_relative_to, original)

    def test_shared_backend_has_python38_grammar(self):
        root = Path(__file__).resolve().parent.parent
        for folder in ("core", "stages", "stations"):
            for source in sorted((root / folder).glob("*.py")):
                with self.subTest(source=source.relative_to(root)):
                    ast.parse(source.read_text(encoding="utf-8"), filename=str(source),
                              feature_version=(3, 8))

    def test_progress_alias_can_be_imported_without_pep604_type_union(self):
        from core.study import Progress
        self.assertIn("Callable", str(Progress))


if __name__ == "__main__":
    unittest.main()

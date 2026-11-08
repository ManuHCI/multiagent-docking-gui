# -*- mode: python ; coding: utf-8 -*-
"""Python 3.8 / PyInstaller 5 native-Tk legacy edition, separate from modern."""
import importlib.util
from pathlib import Path
import sys

from PyInstaller.utils.hooks import collect_data_files, collect_dynamic_libs, copy_metadata

ROOT = Path(SPECPATH).resolve()
PREFIX = Path(sys.prefix)
SITE = Path(importlib.util.find_spec("openmm").origin).parent.parent
VC = ROOT / ".build-tools/legacy-vc-runtime/files-Minimum"
UCRT = ROOT / ".build-tools/legacy-vc-runtime/ucrt-x64"
if sys.version_info[:3] != (3, 8, 10):
    raise RuntimeError("The legacy build must use official CPython 3.8.10 x64.")
for name in ("msvcp140.dll", "vcruntime140.dll", "vcruntime140_1.dll"):
    if not (VC / name).is_file():
        raise RuntimeError("Extract the signed Microsoft VC++ 2019 runtime before building: " + name)
if not (UCRT / "ucrtbase.dll").is_file():
    raise RuntimeError("The original Windows SDK app-local UCRT files are required.")

datas = [(str(ROOT / "assets"), "assets"), (str(ROOT / "LICENSE"), "."),
         (str(ROOT / "LEGACY_WINDOWS.md"), "."),
         (str(ROOT / ".build-tools/legacy-vc-runtime/Microsoft-Visual-Cpp-2019-License.rtf"), "licenses")]
binaries = [(str(ROOT / "vina.exe"), ".")]
hiddenimports = [
    "legacy", "legacy.compat", "legacy.ui", "_tkinter", "tkinter", "tkinter.ttk",
    "tkinter.filedialog", "tkinter.messagebox", "PIL.ImageTk", "PIL.PngImagePlugin",
    "core.study", "core.references", "core.project_library", "core.contract", "core.project",
    "core.runtime", "core.settings", "core.viewer", "core.workshop", "core.mmcif", "core.pdbqt",
    "core.receptor_pdbqt", "core.web", "core.install", "core.desktop_checks", "core.study_checks",
    "stages.target", "stages.structure", "stages.receptor", "stages.site", "stages.ligands",
    "stages.dock", "stages.admet", "stages.md", "rdkit.Chem.AllChem", "rdkit.Chem.Descriptors",
    "rdkit.Chem.Draw", "rdkit.Chem.MolStandardize.rdMolStandardize", "rdkit.Chem.rdMolAlign",
    "rdkit.Chem.rdFMCS", "rdkit.Chem.rdMolTransforms", "scipy.spatial", "scipy.stats",
    "scipy.optimize", "openmm.app", "openmm.unit", "pdbfixer", "meeko", "certifi", "unittest",
]
for package in ("rdkit", "openmm", "pdbfixer", "meeko", "certifi", "PIL"):
    datas += collect_data_files(package)
    binaries += collect_dynamic_libs(package)
for package in ("numpy", "scipy", "pandas"):
    binaries += collect_dynamic_libs(package)
for name in ("rdkit-pypi", "OpenMM", "pdbfixer", "meeko", "numpy", "scipy", "pandas", "Pillow", "certifi"):
    datas += copy_metadata(name)

# Preserve wheel-native library folders and mirror OpenMM's CPU runtime into
# the same private paths used by core.runtime.configure(). No GPU plugins.
for folder in SITE.glob("*.libs"):
    if any(name in folder.name.lower() for name in ("rdkit", "pandas")):
        for path in folder.rglob("*.dll"):
            binaries.append((str(path), str(path.parent.relative_to(SITE))))
openmm_lib = SITE / "OpenMM.libs/lib"
for path in openmm_lib.glob("*.dll"):
    binaries.append((str(path), "Library/lib"))
for path in (openmm_lib / "plugins").glob("*.dll"):
    if "opencl" not in path.name.lower() and "cuda" not in path.name.lower():
        binaries.append((str(path), "Library/lib/plugins"))
rdkit_data = SITE / "rdkit/Data"
if rdkit_data.is_dir():
    datas.append((str(rdkit_data), "Library/share/RDKit/Data"))
runtime_files = {path.name.lower(): path for folder in (VC, UCRT) for path in folder.glob("*.dll")}
binaries += [(str(path), ".") for path in runtime_files.values()]

excluded = ["nicegui", "webview", "pythonnet", "clr", "clr_loader", "PyQt5", "PyQt6", "PySide2", "PySide6",
            "tensorflow", "torch", "jax", "matplotlib", "IPython", "jupyter", "notebook",
            "rdkit.Chem.Draw.IPythonConsole", "rdkit.Chem.PandasTools", "Bio", "prody", "selenium",
            "playwright", "pytest", "sphinx", "sklearn", "skimage", "cv2", "numba", "llvmlite"]
a = Analysis([str(ROOT / "legacy_desktop.py")], pathex=[str(ROOT)], binaries=binaries,
             datas=datas, hiddenimports=hiddenimports, hookspath=[], hooksconfig={},
             runtime_hooks=[str(ROOT / "packaging/legacy/runtime_hook.py")],
             excludes=excluded, noarchive=False)

# The current build computer has a newer VC runtime in System32. Force every
# VC/UCRT occurrence to the signed, older redistributables, including nested
# wheel copies, rather than silently borrowing the build machine's DLLs.
a.binaries = [(destination, str(runtime_files.get(Path(destination).name.lower(), Path(source))), kind)
              for destination, source, kind in a.binaries]
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="InSilicoDockingStudioLegacy",
          debug=False, bootloader_ignore_signals=False, strip=False, upx=False,
          console=False, disable_windowed_traceback=False,
          icon=str(ROOT / "installer/assets/icon.ico"))
coll = COLLECT(exe, a.binaries, a.datas, strip=False, upx=False, name="InSilicoDockingStudioLegacy")

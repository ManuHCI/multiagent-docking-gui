# -*- mode: python ; coding: utf-8 -*-
"""Native macOS bundle. Run build_macos.py on the target Mac architecture."""
import json
import os
from pathlib import Path
import sys

from PyInstaller.utils.hooks import (
    collect_data_files, collect_dynamic_libs, collect_submodules, copy_metadata,
)

ROOT = Path(SPECPATH).resolve()
PREFIX = Path(sys.prefix)
INPUTS = Path(os.environ["INSILICO_MAC_BUILD_INPUTS"]).resolve()
metadata = json.loads((INPUTS / "build-manifest.json").read_text())
os.environ["MATPLOTLIB"] = "false"
datas = [
    (str(ROOT / "assets"), "assets"), (str(ROOT / "LICENSE"), "."),
    (str(INPUTS / "build-manifest.json"), "."),
    (str(INPUTS / "licenses"), "licenses"),
    (str(ROOT / "packaging/macos/READ_ME_FIRST.txt"), "."),
]
binaries = [(str(PREFIX / "bin/vina"), "bin")]
hiddenimports = [
    "app", "forms", "workshop_ui", "study_ui", "core.study", "core.references", "core.project_library", "core.research",
    "core.contract", "core.install", "core.mmcif", "core.pdbqt", "core.desktop_checks", "core.study_checks",
    "core.project", "core.receptor_pdbqt", "core.selftest", "core.settings",
    "core.station", "core.viewer", "core.web", "core.workshop",
    "stages.target", "stages.structure", "stages.receptor", "stages.site",
    "stages.ligands", "stages.dock", "stages.admet", "stages.md", "stations.swissadme",
    "rdkit.Chem.AllChem", "rdkit.Chem.Descriptors", "rdkit.Chem.Crippen",
    "rdkit.Chem.Draw", "rdkit.Chem.MolStandardize.rdMolStandardize",
    "rdkit.Chem.rdMolAlign", "rdkit.Chem.rdFMCS", "rdkit.Chem.rdMolTransforms",
    "scipy.spatial", "scipy.stats", "scipy.optimize", "openmm.app", "openmm.unit", "pdbfixer", "meeko", "gemmi",
    "webview", "webview.platforms.cocoa", "AppKit", "Foundation", "WebKit", "Quartz", "Security",
    "objc", "UniformTypeIdentifiers", "uvicorn.logging", "uvicorn.loops.asyncio",
    "uvicorn.protocols.http.h11_impl", "uvicorn.protocols.websockets.websockets_impl", "uvicorn.lifespan.on",
    "engineio.async_drivers.aiohttp", "engineio.async_drivers.asgi", "socketio.async_client", "certifi", "unittest",
]
datas += collect_data_files("nicegui", include_py_files=True)
hiddenimports += collect_submodules("nicegui", filter=lambda name: not name.startswith("nicegui.testing"))
for package in ("rdkit", "openmm", "pdbfixer", "meeko", "certifi", "webview"):
    datas += collect_data_files(package)
    binaries += collect_dynamic_libs(package)
for package in ("numpy", "scipy"):
    binaries += collect_dynamic_libs(package)
for distribution in ("nicegui", "pywebview", "pdbfixer", "meeko", "gemmi", "openmm", "rdkit", "numpy", "scipy"):
    datas += copy_metadata(distribution)

# Conda scientific dylibs and dynamically loaded OpenMM CPU plugins. GPU
# toolkits are not required for receptor preparation or CPU docking.
for pattern in ("libRDKit*.dylib", "libOpenMM*.dylib"):
    binaries += [(str(path), "lib") for path in (PREFIX / "lib").glob(pattern)]
plugin_dir = PREFIX / "lib/plugins"
plugins = [path for path in plugin_dir.glob("*.dylib")
           if not any(token in path.name.lower() for token in ("cuda", "opencl", "metal", "hip"))]
if not any("cpu" in path.name.lower() for path in plugins):
    raise RuntimeError("The native OpenMM CPU plugin was not found in the build environment.")
binaries += [(str(path), "lib/plugins") for path in plugins]
rdkit_data = PREFIX / "share/RDKit/Data"
if rdkit_data.is_dir():
    datas.append((str(rdkit_data), "share/RDKit/Data"))

excluded = [
    "PyQt5", "PyQt6", "PySide2", "PySide6", "qtpy", "clr", "pythonnet", "clr_loader",
    "tensorflow", "torch", "torchvision", "torchaudio", "jax", "jaxlib",
    "matplotlib", "pandas", "IPython", "ipykernel", "jupyter", "notebook", "jupyterlab", "pytest", "sphinx",
    "rdkit.Chem.Draw.IPythonConsole", "rdkit.Chem.PandasTools", "Bio", "prody", "selenium", "playwright", "redis",
    "pyarrow", "dask", "distributed", "numba", "llvmlite", "xarray", "plotly", "sklearn", "skimage", "cv2",
    "bokeh", "altair", "streamlit", "tables", "h5py", "statsmodels", "sympy", "polars",
    "webview.platforms.winforms", "webview.platforms.edgechromium", "webview.platforms.qt",
    "webview.platforms.gtk", "webview.platforms.android", "webview.platforms.cef",
]
a = Analysis([str(ROOT / "desktop.py")], pathex=[str(ROOT)], binaries=binaries,
             datas=datas, hiddenimports=sorted(set(hiddenimports)), hookspath=[],
             hooksconfig={}, runtime_hooks=[], excludes=excluded, noarchive=False)
pyz = PYZ(a.pure)
identity = os.environ.get("INSILICO_CODESIGN_IDENTITY") or None
entitlements = str(ROOT / "packaging/macos/entitlements.plist")
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="InSilicoTeachingSuite",
          debug=False, bootloader_ignore_signals=False, strip=False, upx=False,
          console=False, argv_emulation=False, target_arch=metadata["architecture"],
          codesign_identity=identity, entitlements_file=entitlements)
coll = COLLECT(exe, a.binaries, a.datas, strip=False, upx=False, name="InSilicoTeachingSuite")
app = BUNDLE(coll, name="In-Silico Docking Studio.app", icon=str(INPUTS / "icon.icns"),
             bundle_identifier="com.akshaisarathi.insilico-docking-studio", version=metadata["version"],
             info_plist={
                 "CFBundleDisplayName": "In-Silico Docking Studio",
                 "CFBundleShortVersionString": metadata["version"],
                 "CFBundleVersion": metadata["version"],
                 "NSPrincipalClass": "NSApplication", "NSHighResolutionCapable": True,
                 "LSMinimumSystemVersion": "15.0",
                 # The native window uses localhost for the loopback UI.
                 "NSAppTransportSecurity": {"NSAllowsLocalNetworking": True},
                 "NSHumanReadableCopyright": "Dr Manu Kumar Shetty and Dr Pratibha Nadig, 2026",
             })


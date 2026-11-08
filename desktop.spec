# -*- mode: python ; coding: utf-8 -*-
"""Relocatable Windows directory bundle; build with build_desktop.py."""

import importlib.util
import os
import sys
from pathlib import Path

from PyInstaller.utils.hooks import (
    collect_data_files,
    collect_dynamic_libs,
    collect_submodules,
    copy_metadata,
)

ROOT = Path(SPECPATH).resolve()
PREFIX = Path(sys.prefix)
TOOLS = ROOT / ".build-tools"
os.environ["MATPLOTLIB"] = "false"

datas = [(str(ROOT / "assets"), "assets"), (str(ROOT / "LICENSE"), "."),
         (str(ROOT / "DESKTOP.md"), ".")]
binaries = [(str(ROOT / "vina.exe"), ".")]
hiddenimports = [
    "app", "forms", "workshop_ui", "study_ui", "core.study", "core.references", "core.project_library",
    "core.contract", "core.install", "core.mmcif", "core.pdbqt",
    "core.project", "core.receptor_pdbqt", "core.selftest", "core.settings",
    "core.station", "core.viewer", "core.web", "core.workshop",
    "stages.target", "stages.structure", "stages.receptor", "stages.site",
    "stages.ligands", "stages.dock", "stages.admet", "stages.md",
    "stations.swissadme",
    "rdkit.Chem.AllChem", "rdkit.Chem.Descriptors", "rdkit.Chem.Crippen",
    "rdkit.Chem.Draw", "rdkit.Chem.MolStandardize.rdMolStandardize",
    "rdkit.Chem.rdMolAlign", "rdkit.Chem.rdFMCS", "rdkit.Chem.rdMolTransforms",
    "scipy.spatial", "scipy.stats", "scipy.optimize",
    "openmm.app", "openmm.unit", "pdbfixer", "meeko",
    "webview", "webview.platforms.winforms", "webview.platforms.edgechromium",
    "clr", "pythonnet", "clr_loader", "uvicorn.logging",
    "uvicorn.loops.asyncio", "uvicorn.protocols.http.h11_impl",
    "uvicorn.protocols.websockets.websockets_impl", "uvicorn.lifespan.on",
    "engineio.async_drivers.aiohttp", "engineio.async_drivers.asgi",
    "socketio.async_client", "certifi", "cgi", "unittest",
]

# NiceGUI loads elements, JavaScript and templates dynamically. Retain its
# package files, including component source files used to locate JS resources.
datas += collect_data_files("nicegui", include_py_files=True)
hiddenimports += collect_submodules("nicegui", filter=lambda name: not name.startswith("nicegui.testing"))

for package in ("rdkit", "openmm", "pdbfixer", "meeko", "certifi", "webview", "pythonnet", "clr_loader"):
    datas += collect_data_files(package)
    binaries += collect_dynamic_libs(package)

# Pip distributions can place their native dependencies beside their package.
# The standard NumPy/SciPy hooks add their own BLAS libraries as well.
for package in ("numpy", "scipy"):
    binaries += collect_dynamic_libs(package)

for distribution in ("nicegui", "pywebview", "pythonnet", "clr_loader", "pdbfixer", "meeko", "openmm", "rdkit", "numpy", "scipy"):
    try:
        datas += copy_metadata(distribution)
    except Exception as error:
        print(f"Metadata unavailable for {distribution}: {error}")

# Conda native libraries live outside site-packages. Explicit entries retain
# their relocatable layout; PyInstaller resolves their dependent DLLs.
library = PREFIX / "Library"
if library.is_dir():
    for pattern in ("RDKit*.dll", "boost*.dll", "libboost*.dll"):
        binaries += [(str(path), "Library/bin") for path in (library / "bin").glob(pattern)]
    binaries += [(str(path), "Library/lib") for path in (library / "lib").glob("OpenMM*.dll")]
    # Preparation runs on CPU. Long MD jobs are exported, so do not ship the
    # CUDA toolkit and GPU kernels with every student's learning app.
    binaries += [(str(path), "Library/lib/plugins") for path in (library / "lib/plugins").glob("*.dll")
                 if "CUDA" not in path.name and "OpenCL" not in path.name]
    rdkit_share = library / "share/RDKit"
    if rdkit_share.is_dir():
        datas.append((str(rdkit_share / "Data"), "Library/share/RDKit/Data"))

# These optional integrations are unrelated to the teaching workflow and can
# otherwise pull a complete scientific IDE or a second GUI toolkit into it.
excluded = [
    "PyQt5", "PyQt6", "PySide2", "PySide6", "qtpy",
    "tensorflow", "torch", "torchvision", "torchaudio", "jax", "jaxlib",
    "matplotlib", "pandas", "IPython", "ipykernel", "jupyter", "notebook",
    "jupyterlab", "pytest", "sphinx", "setuptools",
    "rdkit.Chem.Draw.IPythonConsole", "rdkit.Chem.PandasTools",
    "prody", "Bio", "selenium", "playwright", "redis",
    "pyarrow", "dask", "distributed", "numba", "llvmlite", "xarray",
    "plotly", "pyecharts", "sklearn", "skimage", "cv2", "bokeh", "altair",
    "streamlit", "tables", "h5py", "statsmodels", "sympy", "polars",
    "panel", "pyviz_comms", "intake", "holoviews", "hvplot", "datashader",
    "huggingface_hub", "hf_xet", "transformers", "datasets", "tokenizers",
    "onnx", "onnxruntime", "sentry_sdk", "boto3", "botocore", "s3fs",
    "nbconvert", "nbformat", "jupyter_client", "jupyter_core", "zmq",
    "webview.platforms.qt", "webview.platforms.gtk", "webview.platforms.cocoa",
]

a = Analysis(
    [str(ROOT / "desktop.py")],
    pathex=[str(ROOT), str(TOOLS)],
    binaries=binaries,
    datas=datas,
    hiddenimports=sorted(set(hiddenimports)),
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=excluded,
    noarchive=False,
)
# Conda can contain placeholder API-set DLLs. Windows 10/11 supplies these
# system API sets; one-byte placeholders are not executable libraries.
a.binaries = [entry for entry in a.binaries
              if not (Path(entry[1]).suffix.lower() == ".dll" and Path(entry[1]).stat().st_size < 100)]
pyz = PYZ(a.pure)
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="InSilicoTeachingSuite",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    icon=str(ROOT / "installer/assets/icon.ico"),
    contents_directory="_internal",
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name="InSilicoTeachingSuite",
)

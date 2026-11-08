"""Resource and writable-data locations for source and packaged editions."""
from __future__ import annotations

import os
import hashlib
import shutil
import sys
from pathlib import Path

_dll_handles: list = []


def resource_root() -> Path:
    return Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent.parent))


def data_root() -> Path:
    override = os.environ.get("INSILICO_DATA_DIR", "").strip()
    if override:
        return Path(override).resolve()
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Application Support" / "InSilicoTeachingSuite"
    base = Path(os.environ.get("LOCALAPPDATA", str(Path.home())))
    return base / "InSilicoTeachingSuite"


def workshop_examples() -> Path:
    """Projects rebase their database on open, so templates need writable copies."""
    source = resource_root() / "assets" / "workshop"
    if not getattr(sys, "frozen", False):
        return source
    digest = hashlib.sha256()
    for case in ("egfr", "sv2a"):
        digest.update((source / case / "project.db").read_bytes())
    destination = data_root() / "Examples" / digest.hexdigest()[:16]
    if not (destination / ".ready").exists():
        shutil.copytree(source, destination, dirs_exist_ok=True)
        (destination / ".ready").write_text("Bundled teaching cases copied for this user.\n", encoding="utf-8")
    return destination


def configure() -> None:
    """Find private scientific libraries before any scientific imports."""
    if not getattr(sys, "frozen", False):
        return
    root = resource_root()
    library_dirs = ([root / "bin", root / "lib"] if sys.platform != "win32"
                    else [root / "Library" / "bin", root / "Library" / "lib"])
    os.environ["PATH"] = os.pathsep.join(str(p) for p in library_dirs) + os.pathsep + os.environ.get("PATH", "")
    for directory in (root, root / "Library" / "bin", root / "Library" / "lib",
                      root / "Library" / "lib" / "plugins"):
        if sys.platform == "win32" and directory.is_dir() and hasattr(os, "add_dll_directory"):
            _dll_handles.append(os.add_dll_directory(str(directory)))
    plugin_paths = ([root / "lib" / "plugins", root / "openmm" / "plugins"]
                    if sys.platform != "win32" else [root / "Library" / "lib" / "plugins"])
    for plugins in plugin_paths:
        if plugins.is_dir():
            os.environ["OPENMM_PLUGIN_DIR"] = str(plugins)
            break
    rdkit_paths = ([root / "share" / "RDKit", root / "rdkit"] if sys.platform != "win32"
                   else [root / "Library" / "share" / "RDKit"])
    for rdkit_data in rdkit_paths:
        if (rdkit_data / "Data").is_dir():
            os.environ["RDBASE"] = str(rdkit_data)
            break
    certificates = root / "certifi" / "cacert.pem"
    if certificates.is_file():
        os.environ.setdefault("SSL_CERT_FILE", str(certificates))
        os.environ.setdefault("REQUESTS_CA_BUNDLE", str(certificates))

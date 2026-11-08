"""Locate private legacy scientific libraries before importing the app."""
import os
from pathlib import Path
import sys

root = Path(sys._MEIPASS)
handles = []
for directory in (root, root / "Library/bin", root / "Library/lib"):
    if directory.is_dir() and hasattr(os, "add_dll_directory"):
        handles.append(os.add_dll_directory(str(directory)))
sys._insilico_legacy_dll_handles = handles
os.environ["PATH"] = os.pathsep.join(str(path) for path in
    (root, root / "Library/bin", root / "Library/lib")) + os.pathsep + os.environ.get("PATH", "")
os.environ["OPENMM_PLUGIN_DIR"] = str(root / "Library/lib/plugins")
if (root / "Library/share/RDKit/Data").is_dir():
    os.environ["RDBASE"] = str(root / "Library/share/RDKit")
certificates = root / "certifi/cacert.pem"
if certificates.is_file():
    os.environ.setdefault("SSL_CERT_FILE", str(certificates))
    os.environ.setdefault("REQUESTS_CA_BUNDLE", str(certificates))

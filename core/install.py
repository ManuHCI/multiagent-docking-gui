"""
How to install each bundled tool.

Kept in one place because an error message that names a missing tool
without saying how to get it is only half an error message. The first
version said "reinstall the software", which is sound advice once there
is an installer and useless to anyone running from source - which is
everyone, for now.
"""

from __future__ import annotations
import os

# tool import name -> (friendly name, the command that installs it)
INSTALL_COMMANDS: dict[str, tuple[str, str]] = {
    "rdkit": ("RDKit", "conda install -c conda-forge rdkit"),
    "openmm": ("OpenMM", "conda install -c conda-forge openmm"),
    "pdbfixer": ("PDBFixer", "conda install -c conda-forge pdbfixer"),
    "meeko": ("Meeko", "pip install meeko"),
    "vina": ("AutoDock Vina", "conda install -c conda-forge vina"),
    "numpy": ("NumPy", "conda install -c conda-forge numpy"),
    "py3Dmol": ("py3Dmol", "pip install py3Dmol"),
    "prolif": ("ProLIF", "pip install prolif"),
    "MDAnalysis": ("MDAnalysis", "pip install MDAnalysis"),
    "admet_ai": ("ADMET-AI", "pip install admet-ai"),
    "nicegui": ("NiceGUI", "pip install nicegui"),
}


def how_to_install(*tools: str) -> str:
    """One line per missing tool, ready to paste into Anaconda Prompt."""
    lines = []
    for tool in tools:
        name, command = INSTALL_COMMANDS.get(
            tool, (tool, f"conda install -c conda-forge {tool}"))
        lines.append(command)
    # Collapse conda installs into a single command where possible.
    conda = [l.split()[-1] for l in lines
             if l.startswith("conda install -c conda-forge")]
    pip = [l.split()[-1] for l in lines if l.startswith("pip install")]
    out = []
    if "vina" in tools and os.name == "nt":
        out.append("Download the official Windows Vina executable from "
                   "https://github.com/ccsb-scripps/AutoDock-Vina/releases "
                   "and place it beside app.py as vina.exe")
        conda = [name for name in conda if name != "vina"]
    if conda:
        out.append("conda install -c conda-forge " + " ".join(conda))
    if pip:
        out.append("pip install " + " ".join(pip))
    return "   then   ".join(out)

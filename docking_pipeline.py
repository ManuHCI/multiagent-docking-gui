"""
docking_pipeline.py
===================
AUTO-GENERATED from MultiAgent_Docking.ipynb for use with docking_gui.py.

What was changed vs. the original notebook:
  1. Hardcoded OPENAI_API_KEY line commented out (security).
     Set the key in your environment instead:
         Windows:  set OPENAI_API_KEY=sk-...
  2. IPython-only lines (from IPython.display import / get_ipython()) commented out.
  3. The "## 5. Run the full pipeline" cells removed --
     this file is now an importable module, not a script.

Entry point for the GUI:
    from docking_pipeline import executor_agent, BLACKBOARD, PROJECT_ROOT
    executor_agent("Levetiracetam", ["6CM4", "6A93"])
"""

#!/usr/bin/env python
# coding: utf-8

# # Multi-Agent Molecular Docking Pipeline
# 
# **Author:** Dr. Manu Kumar Shetty, MAMC & Lok Nayak Hospital, New Delhi
# **Pattern:** Planner → Executor → Specialized Agents (mirrors the C1M5 agentic workflow)
# 
# ## What this notebook does
# 
# Given:
# - a **ligand name** (e.g. `"Levetiracetam"`)
# - a **list of receptor PDB IDs** (e.g. `["6CM4", "6A93"]`)
# 
# ...it fully automates:
# 
# 1. 🗂️ **Setup Agent** — create project folder tree + install dependencies
# 2. 💊 **Ligand Agent** — download SDF from PubChem → PDB → PDBQT (Meeko)
# 3. 🧬 **Receptor Agent** — download PDB from RCSB → strip waters/ions/co-crystal ligand → compute grid box center from native ligand → PDBQT (MGLTools or Meeko fallback)
# 4. 🎯 **Docking Agent** — write Vina config → run AutoDock Vina → parse scores
# 5. 📄 **Report Agent** — compile a Markdown report with best pose scores and manuscript-style sentences
# 
# ## Design principle
# 
# The **science is deterministic**. Downloads, cleaning, format conversions, and Vina runs are plain Python functions — not LLM outputs. The LLM only:
# - plans the sequence of steps,
# - routes each step to the right specialized agent,
# - interprets results and writes the final report.
# 
# This keeps docking reproducible while still letting you add natural-language requests like "also dock on 5-HT2A" without rewriting the pipeline.
# 

# ## 1. Imports and LLM client setup
# 
# Same `aisuite` client pattern as your C1M5 assignment.

# In[ ]:


#pip install aisuite

# In[ ]:


#!pip install openai

# In[ ]:


#pip install google-generativeai

# In[ ]:


#pip install anthropic

# In[1]:


import os
# NOTE: this tool uses Ollama (local LLM) and does NOT require an API key.
# If you want to use a cloud provider instead, set the relevant env var
# *outside* this file (never hardcode keys in source code).

# In[2]:


# =============================================================================
# LLM backend configuration -- HARDCODED TO OLLAMA (fully offline)
# =============================================================================
# Uses a local Ollama server at http://localhost:11434
#
# Prerequisite (one-time):
#   1. Install Ollama:   https://ollama.com/download/windows
#   2. Pull a model:     ollama pull llama3.2:3b
#
# To use a different local model, change DEFAULT_MODEL below. Any model pulled
# via `ollama pull <name>` works; the aisuite model-string format is
# "ollama:<name>" -- e.g. "ollama:llama3.2:8b", "ollama:qwen2.5:7b",
# "ollama:phi4-mini", "ollama:gemma3:4b".
# =============================================================================
DEFAULT_MODEL = "ollama:llama3.2:3b"

from aisuite import Client
# Pass a generous timeout so Ollama has time to load the model into memory
# on the first call (cold-start can take 30-120 seconds on CPU-only machines).
CLIENT = Client({"ollama": {"timeout": 600}})

# Quick test
# --- COMMENTED OUT BY GUI EXPORT: say-hello LLM test call at import ---
# r = CLIENT.chat.completions.create(
    # model="openai:gpt-4o-mini",
    # messages=[{"role": "user", "content": "say hello"}],
    # temperature=0,
# )
# print(r.choices[0].message.content)

# In[ ]:




# In[3]:


# =========================
# Standard library
# =========================
import os
import re
import ast
import json
import shutil
import subprocess
from datetime import datetime
from pathlib import Path

# =========================
# Third-party
# =========================
import requests
import pandas as pd
# from IPython.display import Markdown, display  # <-- IPython-only, commented out for GUI

from aisuite import Client

# --- COMMENTED OUT BY GUI EXPORT: duplicate CLIENT init ---
# CLIENT = Client()

# =========================
# Project root — EDIT THIS
# =========================
PROJECT_ROOT = Path(r"E:\Docking\Levetiracetam_Project\MultiAgent_Project")
PROJECT_ROOT.mkdir(parents=True, exist_ok=True)
os.chdir(PROJECT_ROOT)
print("Project root:", PROJECT_ROOT)


# ## 2. Deterministic tools (the actual science)
# 
# These are plain Python functions — no LLM involved. Agents call these to do real work. If any of these fail, the pipeline fails cleanly with a real error message instead of hallucinated "success".
# 
# ### 2.1 Folder setup + dependency install

# In[4]:


def tool_create_folders(root: Path) -> dict:
    """Create the standard docking project folder tree."""
    folders = {
        "ligands":            root / "ligands",
        "receptors_raw":      root / "receptors" / "raw",
        "receptors_prepared": root / "receptors" / "prepared",
        "configs":            root / "configs",
        "outputs":            root / "outputs",
        "reports":            root / "reports",
    }
    for p in folders.values():
        p.mkdir(parents=True, exist_ok=True)
    return {k: str(v) for k, v in folders.items()}


def tool_check_dependencies() -> dict:
    """Check which required tools are available.

    Required Python packages: rdkit, meeko, requests, pandas
    Required external tools:  obabel (Open Babel), vina (AutoDock Vina)
    Optional:                 MGLTools prepare_receptor4.py (for receptor PDBQT)
    """
    report = {}

    # Python packages
    for pkg in ["rdkit", "meeko", "requests", "pandas"]:
        try:
            __import__(pkg)
            report[pkg] = "OK"
        except ImportError:
            report[pkg] = "MISSING — run: pip install " + pkg

    # External binaries
    for binary in ["obabel", "vina"]:
        path = shutil.which(binary)
        report[binary] = path if path else "MISSING — install and add to PATH"

    return report


def tool_install_missing(packages: list) -> str:
    """Install missing Python packages with pip."""
    if not packages:
        return "Nothing to install."
    cmd = ["pip", "install"] + packages
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.stdout[-500:] + "\n" + r.stderr[-500:]


# ### 2.2 Ligand tools — PubChem download + format conversion
# 
# PubChem REST API: `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{NAME}/SDF` gives us a 2D SDF. We then:
# - convert SDF → PDB with Open Babel (adds 3D coordinates + H)
# - convert PDB → PDBQT with Meeko (adds Gasteiger charges + AutoDock atom types)

# In[19]:


PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def tool_download_ligand_from_pubchem(name: str, out_sdf: Path) -> Path:
    """Download a compound by name from PubChem as 3D SDF.

    Tries 3D SDF first (better starting geometry); falls back to 2D.
    """
    out_sdf = Path(out_sdf)
    out_sdf.parent.mkdir(parents=True, exist_ok=True)

    for record_type in ["3d", "2d"]:
        url = f"{PUBCHEM_BASE}/compound/name/{name}/SDF?record_type={record_type}"
        r = requests.get(url, timeout=60)
        if r.status_code == 200 and r.text.strip():
            out_sdf.write_text(r.text)
            print(f"  Downloaded {name} ({record_type.upper()} SDF) -> {out_sdf}")
            return out_sdf

    raise RuntimeError(f"PubChem did not return SDF for '{name}'. "
                       f"Check spelling or try a synonym / CID.")


def tool_sdf_to_pdbqt_ligand(in_sdf: Path, out_pdbqt: Path) -> Path:
    """Prepare ligand PDBQT directly from an SDF file using RDKit + Meeko.
    
    Skips the intermediate PDB file entirely.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation

    in_sdf, out_pdbqt = Path(in_sdf), Path(out_pdbqt)

    # Read the first molecule from the SDF
    suppl = Chem.SDMolSupplier(str(in_sdf), removeHs=False)
    mol = suppl[0]
    if mol is None:
        raise RuntimeError(f"RDKit could not parse {in_sdf}")

    # Add explicit hydrogens with 3D coords
    mol = Chem.AddHs(mol, addCoords=True)

    # Embed 3D coords if the SDF only has 2D
    if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D():
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)

    prep = MoleculePreparation()
    prep.prepare(mol)
    out_pdbqt.write_text(prep.write_pdbqt_string())
    return out_pdbqt





def tool_prepare_ligand(name: str, ligand_dir: Path) -> dict:
    safe = name.lower().replace(" ", "_")
    sdf   = ligand_dir / f"{safe}.sdf"
    pdbqt = ligand_dir / f"{safe}.pdbqt"

    tool_download_ligand_from_pubchem(name, sdf)
    tool_sdf_to_pdbqt_ligand(sdf, pdbqt)

    return {"name": name, "sdf": str(sdf), "pdbqt": str(pdbqt)}


# ### 2.3 Receptor tools — RCSB download + cleaning + grid box + PDBQT
# 
# Pipeline for each PDB ID:
# 1. **Download** the raw `.pdb` from RCSB (`https://files.rcsb.org/download/{ID}.pdb`)
# 2. **Inspect** — list all HETATM residue codes that aren't standard amino acids or water
# 3. **Clean** — remove HOH + all non-standard HETATMs (co-crystal ligand, ions, buffers like PEG, cholesterol OLA, etc.)
# 4. **Grid box center** — compute the centroid of the *native* co-crystal ligand *before* stripping it. This is your binding pocket center.
# 5. **PDBQT** — try MGLTools `prepare_receptor4.py` first (gold standard), fall back to Open Babel if MGLTools isn't installed.

# In[20]:


RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

# Standard residues we always keep
STANDARD_AA = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    "MSE",  # selenomethionine — usually kept
}


def tool_download_receptor_from_rcsb(pdb_id: str, out_pdb: Path) -> Path:
    """Download a raw PDB file from RCSB."""
    pdb_id = pdb_id.upper()
    out_pdb = Path(out_pdb)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)

    r = requests.get(RCSB_PDB_URL.format(pdb_id=pdb_id), timeout=60)
    if r.status_code != 200:
        raise RuntimeError(f"RCSB returned {r.status_code} for {pdb_id}")

    out_pdb.write_text(r.text)
    print(f"  Downloaded {pdb_id} -> {out_pdb}  ({len(r.text):,} bytes)")
    return out_pdb


def tool_inspect_hetatms(pdb_file: Path) -> set:
    """Return set of non-standard HETATM residue codes (excluding water).

    These are candidates for removal: co-crystal ligand, ions, buffers.
    """
    hets = set()
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("HETATM"):
                res = line[17:20].strip()
                if res not in STANDARD_AA and res != "HOH":
                    hets.add(res)
    return hets


def tool_guess_native_ligand(hets: set) -> str | None:
    """Heuristic: the native co-crystal ligand is usually the non-standard
    HETATM that is NOT a common ion / buffer.

    Returns the guessed residue code, or None if ambiguous.
    """
    # Common crystallization additives — NOT the native ligand
    additives = {
        "PEG","PGE","PG4","PE4","EDO","GOL","MPD","DMS","TRS","BOG",
        "OLA","OLB","OLC","CLR","PLM","STE","MYR","LMT",      # lipids / detergents
        "NA","K","CL","CA","MG","ZN","FE","MN","CU","NI","CD", # ions
        "SO4","PO4","NO3","ACT","FMT","CIT","ACE","IOD","BR",  # small anions
    }
    candidates = hets - additives
    if len(candidates) == 1:
        return next(iter(candidates))
    # If multiple: return the one with the longest name (usually the drug-like one)
    if candidates:
        return sorted(candidates, key=len, reverse=True)[0]
    return None


def tool_clean_receptor(in_pdb: Path, out_pdb: Path,
                        remove_resnames: set) -> Path:
    """Strip water + specified HETATM residues. Keep ATOM + kept HETATMs."""
    in_pdb, out_pdb = Path(in_pdb), Path(out_pdb)
    remove = {"HOH", "WAT"} | set(remove_resnames)

    with open(in_pdb) as fin, open(out_pdb, "w") as fout:
        for line in fin:
            rec = line[:6].strip()
            if rec == "HETATM":
                resname = line[17:20].strip()
                if resname in remove:
                    continue
            fout.write(line)
    return out_pdb


def tool_grid_center_from_ligand(raw_pdb: Path, ligand_resname: str,
                                 box_size: float = 22.0) -> dict:
    """Compute docking grid box center as centroid of the native ligand.

    Returns {center_x, center_y, center_z, size_x, size_y, size_z}.
    """
    coords = []
    with open(raw_pdb) as f:
        for line in f:
            if line.startswith("HETATM") and line[17:20].strip() == ligand_resname:
                coords.append((float(line[30:38]),
                               float(line[38:46]),
                               float(line[46:54])))
    if not coords:
        raise RuntimeError(
            f"No atoms found for residue '{ligand_resname}' in {raw_pdb}. "
            f"Cannot compute grid center."
        )
    df = pd.DataFrame(coords, columns=["X", "Y", "Z"])
    return {
        "center_x": round(df["X"].mean(), 3),
        "center_y": round(df["Y"].mean(), 3),
        "center_z": round(df["Z"].mean(), 3),
        "size_x":   box_size,
        "size_y":   box_size,
        "size_z":   box_size,
        "native_ligand": ligand_resname,
        "n_atoms": len(coords),
    }


def tool_pdb_to_pdbqt_receptor(in_pdb: Path, out_pdbqt: Path,
                               mgltools_python: str | None = None,
                               mgltools_prepare_script: str | None = None) -> Path:
    """Convert cleaned receptor PDB -> PDBQT.

    Prefers MGLTools prepare_receptor4.py (gold standard for Vina).
    Falls back to Open Babel if MGLTools paths not given or not found.
    """
    in_pdb, out_pdbqt = Path(in_pdb), Path(out_pdbqt)

    # Try MGLTools first
    if mgltools_python and mgltools_prepare_script \
       and Path(mgltools_python).exists() and Path(mgltools_prepare_script).exists():
        cmd = [
            mgltools_python, mgltools_prepare_script,
            "-r", str(in_pdb),
            "-o", str(out_pdbqt),
        ]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode == 0 and out_pdbqt.exists():
            print(f"  Receptor PDBQT via MGLTools -> {out_pdbqt}")
            return out_pdbqt
        print(f"  MGLTools failed, falling back to obabel.\n{r.stderr[:300]}")

    # Fallback: Open Babel (-xr = rigid receptor)
    cmd = ["obabel", str(in_pdb), "-O", str(out_pdbqt), "-xr"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not out_pdbqt.exists():
        raise RuntimeError(f"obabel PDB->PDBQT failed:\n{r.stderr}")
    print(f"  Receptor PDBQT via obabel -> {out_pdbqt}")
    return out_pdbqt


def tool_prepare_receptor(pdb_id: str, receptor_raw_dir: Path,
                          receptor_prepared_dir: Path,
                          mgltools_python: str | None = None,
                          mgltools_prepare_script: str | None = None) -> dict:
    """Full receptor pipeline for one PDB ID."""
    pdb_id = pdb_id.upper()
    raw        = receptor_raw_dir / f"{pdb_id}.pdb"
    clean      = receptor_prepared_dir / f"{pdb_id}_clean.pdb"
    pdbqt      = receptor_prepared_dir / f"{pdb_id}.pdbqt"

    tool_download_receptor_from_rcsb(pdb_id, raw)
    hets = tool_inspect_hetatms(raw)
    native = tool_guess_native_ligand(hets)

    # Remove everything non-standard EXCEPT we'll also strip the native ligand
    # so the binding pocket is empty for docking
    to_remove = set(hets)  # strip all non-standard HETATMs + waters
    tool_clean_receptor(raw, clean, to_remove)

    grid = (tool_grid_center_from_ligand(raw, native) if native
            else {"center_x": None, "center_y": None, "center_z": None,
                  "size_x": 22, "size_y": 22, "size_z": 22,
                  "native_ligand": None, "n_atoms": 0})

    tool_pdb_to_pdbqt_receptor(clean, pdbqt,
                               mgltools_python, mgltools_prepare_script)

    return {
        "pdb_id": pdb_id,
        "raw":    str(raw),
        "clean":  str(clean),
        "pdbqt":  str(pdbqt),
        "hetatms_found":   sorted(hets),
        "hetatms_removed": sorted(to_remove),
        "native_ligand":   native,
        "grid":            grid,
    }


# ### 2.4 Docking tools — Vina config, run, parse
# 
# AutoDock Vina reads a plain-text config file and writes the docked poses to a `.pdbqt`. Scores for each pose are embedded as `REMARK VINA RESULT` lines — we parse those.

# In[21]:


def tool_write_vina_config(receptor_pdbqt, ligand_pdbqt, grid, config_path,
                           exhaustiveness: int = 24,
                           num_modes: int = 20,
                           energy_range: int = 4) -> Path:
    """Write a Vina config file. Paths use forward slashes to avoid Windows backslash issues."""
    config_path = Path(config_path)
    config_path.parent.mkdir(parents=True, exist_ok=True)

    # Vina is picky on Windows — forward slashes, no quotes, absolute paths
    recp = str(Path(receptor_pdbqt).resolve()).replace("\\", "/")
    lig  = str(Path(ligand_pdbqt).resolve()).replace("\\", "/")

    text = f"""receptor = {recp}
ligand = {lig}

center_x = {grid['center_x']}
center_y = {grid['center_y']}
center_z = {grid['center_z']}

size_x = {grid['size_x']}
size_y = {grid['size_y']}
size_z = {grid['size_z']}

exhaustiveness = {exhaustiveness}
num_modes = {num_modes}
energy_range = {energy_range}
"""
    config_path.write_text(text)
    return config_path


def tool_run_vina(config_path: Path, out_pdbqt: Path, log_path: Path) -> dict:
    """Run AutoDock Vina. Works with both Vina 1.1.2 and 1.2.x.
    
    Vina 1.2 removed the --log option, so we capture stdout to the log file ourselves.
    """
    cmd = ["vina",
           "--config", str(config_path),
           "--out",    str(out_pdbqt)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    
    # Write our own log file (captures what used to go to --log)
    Path(log_path).write_text(
        f"=== STDOUT ===\n{r.stdout}\n\n=== STDERR ===\n{r.stderr}\n"
    )
    
    return {
        "returncode": r.returncode,
        "stdout":     r.stdout[-2000:],
        "stderr":     r.stderr[-2000:],
        "out_pdbqt":  str(out_pdbqt),
        "log":        str(log_path),
    }


def tool_parse_vina_scores(out_pdbqt: Path) -> list:
    """Parse all VINA RESULT lines. Returns list of dicts with pose scores."""
    poses = []
    with open(out_pdbqt) as f:
        for line in f:
            if "VINA RESULT" in line:
                parts = line.split()
                # Format: REMARK VINA RESULT:  affinity  rmsd_lb  rmsd_ub
                try:
                    poses.append({
                        "affinity_kcal_per_mol": float(parts[3]),
                        "rmsd_lower_bound":      float(parts[4]),
                        "rmsd_upper_bound":      float(parts[5]),
                    })
                except (ValueError, IndexError):
                    continue
    return poses


def tool_pdbqt_to_pdb(in_pdbqt: Path, out_pdb: Path) -> Path:
    """Convert docked pose PDBQT -> PDB for visualization (PyMOL / Discovery Studio)."""
    in_pdbqt, out_pdb = Path(in_pdbqt), Path(out_pdb)
    cmd = ["obabel", str(in_pdbqt), "-O", str(out_pdb)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"obabel PDBQT->PDB failed:\n{r.stderr}")
    return out_pdb


def tool_dock_pair(ligand_info: dict, receptor_info: dict,
                   configs_dir: Path, outputs_dir: Path) -> dict:
    """Dock one ligand against one receptor. Returns full result dict."""
    lig_name   = ligand_info["name"].lower().replace(" ", "_")
    pdb_id     = receptor_info["pdb_id"]
    tag        = f"{pdb_id}_{lig_name}"

    cfg    = configs_dir / f"{tag}_config.txt"
    out_q  = outputs_dir / f"{tag}_out.pdbqt"
    out_p  = outputs_dir / f"{tag}_out.pdb"
    log    = outputs_dir / f"{tag}_log.txt"

    if receptor_info["grid"]["center_x"] is None:
        return {"tag": tag, "status": "skipped",
                "reason": "no native ligand found — grid center unknown"}

    tool_write_vina_config(receptor_info["pdbqt"], ligand_info["pdbqt"],
                           receptor_info["grid"], cfg)
    run = tool_run_vina(cfg, out_q, log)

    if run["returncode"] != 0 or not Path(out_q).exists():
        return {"tag": tag, "status": "failed", "run": run,
                "config": str(cfg)}

    poses = tool_parse_vina_scores(out_q)
    try:
        tool_pdbqt_to_pdb(out_q, out_p)
        pdb_vis = str(out_p)
    except Exception as e:
        pdb_vis = f"conversion failed: {e}"

    return {
        "tag": tag,
        "status": "success",
        "config":      str(cfg),
        "out_pdbqt":   str(out_q),
        "out_pdb":     pdb_vis,
        "log":         str(log),
        "poses":       poses,
        "best_score":  min(p["affinity_kcal_per_mol"] for p in poses) if poses else None,
        "n_poses":     len(poses),
    }


# ## 3. Agents
# 
# Each agent is a thin LLM wrapper around the deterministic tools. The pattern matches your C1M5 structure:
# - **system_prompt** defines the role
# - **messages** = `[system, user]`
# - LLM call with `temperature` tuned per agent
# 
# The agents also print progress, because docking is slow and silent pipelines are scary.

# ### 3.1 Setup Agent
# 
# Creates folders and verifies that `obabel` and `vina` are on PATH. If Python packages are missing, it installs them.

# In[22]:


def setup_agent(task: str, model: str = DEFAULT_MODEL) -> str:
    """Setup agent: creates folder tree, checks & installs dependencies."""
    print("==================================")
    print("🗂️  Setup Agent")
    print("==================================")

    # --- Do the actual work deterministically ---
    folders = tool_create_folders(PROJECT_ROOT)
    deps    = tool_check_dependencies()

    missing_pkgs = [k for k, v in deps.items()
                    if v.startswith("MISSING") and k in {"rdkit", "meeko", "requests", "pandas"}]
    install_log = tool_install_missing(missing_pkgs) if missing_pkgs else "All Python packages OK."

    missing_bin = [k for k in ("obabel", "vina") if deps[k].startswith("MISSING")]

    # --- Let the LLM summarize the result in the agent style ---
    system_prompt = """You are a setup agent for a molecular docking pipeline.
Your job: confirm the workspace is ready and clearly flag any missing binaries
(obabel, vina, MGLTools) that the user must install manually. Be concise.
Format your reply as a short checklist."""

    user_prompt = f"""
Task: {task}

Folders created:
{json.dumps(folders, indent=2)}

Dependency check:
{json.dumps(deps, indent=2)}

Python install log (if any):
{install_log}

Missing external binaries (user must install manually): {missing_bin}

Write a short checklist-style status report for the user.
"""

    response = CLIENT.chat.completions.create(
        model=model,
        messages=[{"role": "system", "content": system_prompt},
                  {"role": "user",   "content": user_prompt}],
        temperature=0.3,
    )
    return response.choices[0].message.content


# ### 3.2 Ligand Agent
# 
# Takes a compound name, downloads from PubChem, converts through SDF → PDB → PDBQT.
# 
# **Tip:** If PubChem returns nothing, check spelling or use a PubChem CID. The agent reports this clearly instead of pretending to succeed.

# In[23]:


def ligand_agent(task: str, model: str = DEFAULT_MODEL) -> dict:
    """Ligand agent: PubChem name -> SDF -> PDB -> PDBQT."""
    print("==================================")
    print("💊 Ligand Agent")
    print("==================================")

    # Extract the compound name from the task using the LLM
    extract_prompt = f"""Extract ONLY the chemical compound name from the task
below. Return just the name, nothing else, no quotes, no explanation.

Task: {task}"""
    r = CLIENT.chat.completions.create(
        model=model,
        messages=[{"role": "user", "content": extract_prompt}],
        temperature=0,
    )
    compound = r.choices[0].message.content.strip().strip('"').strip("'")
    print(f"  Target compound: {compound}")

    # Deterministic pipeline
    try:
        result = tool_prepare_ligand(compound, PROJECT_ROOT / "ligands")
        result["status"] = "success"
        print(f"  ✅ Ligand ready: {result['pdbqt']}")
    except Exception as e:
        result = {"status": "failed", "name": compound, "error": str(e)}
        print(f"  ❌ Ligand prep failed: {e}")

    return result


# ### 3.3 Receptor Agent
# 
# Takes a list of PDB IDs. For each:
# - downloads from RCSB
# - auto-detects the native co-crystal ligand (so we get a real binding pocket center)
# - strips waters, ions, buffers, and native ligand
# - computes grid box center from native ligand coordinates
# - generates receptor PDBQT
# 
# **Heads up:** If MGLTools isn't installed, we fall back to Open Babel. MGLTools is strongly preferred for Vina — set the two paths below if you have it.

# In[24]:


# ---- OPTIONAL: set these if you have MGLTools installed ----
MGLTOOLS_PYTHON = r"C:\Program Files (x86)\MGLTools-1.5.7\python.exe"
MGLTOOLS_PREPARE = r"C:\Program Files (x86)\MGLTools-1.5.7\Lib\site-packages\AutoDockTools\Utilities24\prepare_receptor4.py"


def receptor_agent(task: str, model: str = DEFAULT_MODEL) -> list:
    """Receptor agent: list of PDB IDs -> prepared receptors with grid boxes."""
    print("==================================")
    print("🧬 Receptor Agent")
    print("==================================")

    # Ask LLM to extract a clean list of PDB IDs from the task
    extract_prompt = f"""Extract ALL 4-character PDB IDs from the task below.
Return ONLY a valid Python list of uppercase strings, nothing else.
Example: ["6CM4", "6A93"]

Task: {task}"""
    r = CLIENT.chat.completions.create(
        model=model,
        messages=[{"role": "user", "content": extract_prompt}],
        temperature=0,
    )
    raw = r.choices[0].message.content.strip()
    raw = re.sub(r"^```(?:python)?\n?", "", raw).rstrip("`").strip()
    try:
        pdb_ids = ast.literal_eval(raw)
    except Exception:
        # Fallback: regex
        pdb_ids = re.findall(r"\b[0-9][A-Z0-9]{3}\b", raw.upper())
    print(f"  PDB IDs to prepare: {pdb_ids}")

    results = []
    for pid in pdb_ids:
        print(f"\n  --- {pid} ---")
        try:
            info = tool_prepare_receptor(
                pid,
                PROJECT_ROOT / "receptors" / "raw",
                PROJECT_ROOT / "receptors" / "prepared",
                MGLTOOLS_PYTHON if Path(MGLTOOLS_PYTHON).exists() else None,
                MGLTOOLS_PREPARE if Path(MGLTOOLS_PREPARE).exists() else None,
            )
            info["status"] = "success"
            print(f"  ✅ {pid}: native={info['native_ligand']}, "
                  f"center=({info['grid']['center_x']}, {info['grid']['center_y']}, {info['grid']['center_z']})")
        except Exception as e:
            info = {"pdb_id": pid, "status": "failed", "error": str(e)}
            print(f"  ❌ {pid} failed: {e}")
        results.append(info)

    return results


# In[25]:


import subprocess
# --- COMMENTED OUT BY GUI EXPORT: vina version check at import ---
# r = subprocess.run(["vina", "--version"], capture_output=True, text=True)
# print(r.stdout)
# print(r.stderr)

# ### 3.4 Docking Agent
# 
# Takes the prepared ligand + all prepared receptors. Runs Vina for each pair. Parses all poses. Returns a table of results.

# In[26]:


def docking_agent(ligand_info: dict, receptor_infos: list,
                  model: str = DEFAULT_MODEL) -> list:
    """Docking agent: dock one ligand against every prepared receptor."""
    print("==================================")
    print("🎯 Docking Agent")
    print("==================================")

    if ligand_info.get("status") != "success":
        print("  ❌ Ligand is not ready, cannot dock.")
        return []

    results = []
    for rec in receptor_infos:
        if rec.get("status") != "success":
            print(f"  ⏭️  Skipping {rec.get('pdb_id')} (not prepared)")
            continue

        print(f"\n  Docking {ligand_info['name']} → {rec['pdb_id']} ...")
        res = tool_dock_pair(ligand_info, rec,
                             PROJECT_ROOT / "configs",
                             PROJECT_ROOT / "outputs")
        if res["status"] == "success":
            print(f"  ✅ Best score: {res['best_score']} kcal/mol "
                  f"({res['n_poses']} poses)")
        else:
            print(f"  ❌ {res['status']}: {res.get('reason') or res.get('run', {}).get('stderr', '')[:200]}")
        results.append(res)

    return results


# ### 3.5 Report Agent
# 
# Takes all the intermediate results and writes a Markdown research report — same role as your `writer_agent` + `editor_agent` combined, scoped to docking output.

# In[27]:


def report_agent(ligand_info: dict, receptor_infos: list, docking_results: list,
                 model: str = DEFAULT_MODEL) -> str:
    """Report agent: compile results into a manuscript-style Markdown report."""
    print("==================================")
    print("📄 Report Agent")
    print("==================================")

    # Build a compact, LLM-friendly summary
    summary = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "ligand": {
            "name":   ligand_info.get("name"),
            "status": ligand_info.get("status"),
            "pdbqt":  ligand_info.get("pdbqt"),
        },
        "receptors": [
            {"pdb_id": r["pdb_id"],
             "status": r["status"],
             "native_ligand": r.get("native_ligand"),
             "grid": r.get("grid"),
             "hetatms_removed": r.get("hetatms_removed", [])}
            for r in receptor_infos
        ],
        "docking": [
            {"tag": d["tag"], "status": d["status"],
             "best_score": d.get("best_score"),
             "n_poses":    d.get("n_poses"),
             "top_3_poses": d.get("poses", [])[:3]}
            for d in docking_results
        ],
    }

    system_prompt = """You are a scientific writing agent for a molecular docking pipeline.
Produce a clean, publication-oriented Markdown report with the following sections:

1. # Summary (1–2 sentences stating what was docked and the best receptor hit)
2. ## Ligand preparation (compound, PubChem source, final PDBQT path)
3. ## Receptor preparation (one subsection per PDB; list native ligand, grid center, removed HETATMs)
4. ## Docking results (a Markdown table: receptor | best score (kcal/mol) | n poses | status)
5. ## Interpretation — for each successful dock, one manuscript-style sentence using this score banding:
   > -4.5 weak | -4.5 to -5.5 mild | -5.5 to -6.5 moderate | -6.5 to -7.5 interesting | < -7.5 strong
6. ## Caveats (3 bullets max: docking ≠ pharmacology; needs in vitro; MGLTools vs obabel receptor prep affects charges)

Be precise with numbers. Do NOT invent scores. Use ONLY the data provided."""

    user_prompt = f"""Here is the pipeline output as JSON:

```json
{json.dumps(summary, indent=2, default=str)}
```

Write the Markdown report."""

    response = CLIENT.chat.completions.create(
        model=model,
        messages=[{"role": "system", "content": system_prompt},
                  {"role": "user",   "content": user_prompt}],
        temperature=0.4,
    )
    report_md = response.choices[0].message.content

    # -------------------------------------------------------------------------
    # Embed figures (if visualization_agent has already run) at the end of the
    # report. Each figure is a relative path that renders inline on GitHub,
    # VS Code, and most Markdown viewers.
    # -------------------------------------------------------------------------
    figures = BLACKBOARD.get("figures") if isinstance(BLACKBOARD.get("figures"), dict) else None
    if figures:
        figure_labels = {
            "heatmap":           "Per-pose affinity heatmap across receptors",
            "bar_chart":         "Best docking score per receptor",
            "pose_distribution": "Distribution of docking scores per receptor",
            "rmsd_vs_affinity":  "RMSD vs affinity (pose clustering)",
        }
        figure_md = ["\n\n---\n\n## Figures\n"]
        for key, caption in figure_labels.items():
            path = figures.get(key)
            if path and Path(path).exists():
                # Use path relative to the report's own folder so links don't break
                try:
                    rel = Path(path).relative_to(PROJECT_ROOT / "reports")
                except ValueError:
                    rel = Path(path)
                figure_md.append(f"### {caption}\n")
                figure_md.append(f"![{key}]({rel.as_posix()})\n")
        if len(figure_md) > 1:  # only append if at least one figure exists
            report_md += "\n".join(figure_md)

    # Save to disk
    report_path = PROJECT_ROOT / "reports" / f"docking_report_{datetime.now():%Y%m%d_%H%M%S}.md"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(report_md, encoding="utf-8")
    print(f"  📄 Report saved: {report_path}")

    return report_md


# In[28]:


# ==========================================================================
# Tool: visualization of docking results
# ==========================================================================
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


# Score-band colors (your original manuscript bands)
SCORE_BANDS = [
    (-4.5, float("inf"),  "#d62728", "weak (>-4.5)"),
    (-5.5, -4.5,          "#ff7f0e", "mild (-4.5 to -5.5)"),
    (-6.5, -5.5,          "#ffbb00", "moderate (-5.5 to -6.5)"),
    (-7.5, -6.5,          "#2ca02c", "interesting (-6.5 to -7.5)"),
    (float("-inf"), -7.5, "#1f77b4", "strong (<-7.5)"),
]


def _band_color(score):
    if score is None:
        return "#cccccc"
    for lo, hi, color, _ in SCORE_BANDS:
        if lo <= score < hi:
            return color
    return "#cccccc"


def tool_plot_heatmap(docking_results, ligand_name, out_path):
    """Heatmap: ligand x receptor best affinities.
    
    With a single ligand this is a 1-row strip; becomes more informative
    when you add reference ligands (haloperidol, risperidone) later.
    """
    successful = [d for d in docking_results if d.get("status") == "success"]
    if not successful:
        print("  No successful dockings — nothing to plot.")
        return None

    receptors = [d["tag"].split("_")[0] for d in successful]
    scores = np.array([[d["best_score"] for d in successful]])

    fig, ax = plt.subplots(figsize=(max(6, len(receptors) * 1.5), 2.5))
    im = ax.imshow(scores, cmap="RdYlGn_r", aspect="auto", vmin=-10, vmax=-4)

    ax.set_xticks(range(len(receptors)))
    ax.set_xticklabels(receptors, fontsize=11)
    ax.set_yticks([0])
    ax.set_yticklabels([ligand_name], fontsize=11)

    # Annotate each cell with the score
    for i, score in enumerate(scores[0]):
        ax.text(i, 0, f"{score:.2f}", ha="center", va="center",
                color="black", fontsize=11, fontweight="bold")

    cbar = plt.colorbar(im, ax=ax, orientation="horizontal",
                        pad=0.3, shrink=0.6)
    cbar.set_label("Binding affinity (kcal/mol) — more negative = stronger",
                   fontsize=10)
    ax.set_title(f"Docking affinity heatmap — {ligand_name}",
                 fontsize=12, pad=12)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  📊 Heatmap saved: {out_path}")
    return out_path


def tool_plot_best_affinity_bar(docking_results, ligand_name, out_path):
    """Bar chart of best affinity per receptor, colored by score band."""
    successful = [d for d in docking_results if d.get("status") == "success"]
    if not successful:
        return None

    receptors = [d["tag"].split("_")[0] for d in successful]
    scores = [d["best_score"] for d in successful]
    colors = [_band_color(s) for s in scores]

    fig, ax = plt.subplots(figsize=(max(6, len(receptors) * 1.2), 5))
    bars = ax.bar(receptors, scores, color=colors, edgecolor="black",
                  linewidth=1.2)

    # Label bars with the numeric score
    for bar, score in zip(bars, scores):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() - 0.15,
                f"{score:.2f}",
                ha="center", va="top",
                fontsize=11, fontweight="bold", color="white")

    # Horizontal dashed lines at each band threshold
    for lo, hi, color, label in SCORE_BANDS:
        if lo != float("-inf") and lo > -10:
            ax.axhline(lo, linestyle="--", color=color, alpha=0.5, linewidth=1)

    ax.set_ylabel("Binding affinity (kcal/mol)", fontsize=12)
    ax.set_title(f"Best docking affinity: {ligand_name} vs receptors",
                 fontsize=13, pad=10)
    ax.invert_yaxis()  # more negative (stronger) at the top
    ax.grid(axis="y", alpha=0.3)

    # Legend of score bands
    handles = [mpatches.Patch(color=c, label=lbl)
               for _, _, c, lbl in SCORE_BANDS]
    ax.legend(handles=handles, loc="lower right", fontsize=9, title="Score band")

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  📊 Bar chart saved: {out_path}")
    return out_path


def tool_plot_pose_distribution(docking_results, ligand_name, out_path):
    """Strip plot: every pose's affinity per receptor.
    
    Tight clustering near the best score = confident docking.
    Wide spread = the ligand doesn't have a clear preferred pose.
    """
    successful = [d for d in docking_results if d.get("status") == "success"]
    if not successful:
        return None

    fig, ax = plt.subplots(figsize=(max(6, len(successful) * 1.5), 6))

    for i, d in enumerate(successful):
        receptor = d["tag"].split("_")[0]
        affinities = [p["affinity_kcal_per_mol"] for p in d["poses"]]
        jitter = np.random.uniform(-0.15, 0.15, size=len(affinities))
        x = np.full(len(affinities), i) + jitter

        ax.scatter(x, affinities, s=60, alpha=0.6,
                   color=_band_color(d["best_score"]),
                   edgecolor="black", linewidth=0.5)
        # Highlight the best pose
        ax.scatter([i], [min(affinities)], s=180, marker="*",
                   color="red", edgecolor="black", linewidth=1, zorder=5,
                   label="Best pose" if i == 0 else "")

    ax.set_xticks(range(len(successful)))
    ax.set_xticklabels([d["tag"].split("_")[0] for d in successful], fontsize=11)
    ax.set_ylabel("Affinity (kcal/mol)", fontsize=12)
    ax.set_title(f"Pose distribution — {ligand_name}\n(all 20 Vina poses per receptor)",
                 fontsize=13, pad=10)
    ax.invert_yaxis()
    ax.grid(axis="y", alpha=0.3)
    ax.legend(loc="lower right", fontsize=10)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  📊 Pose distribution saved: {out_path}")
    return out_path


def tool_plot_rmsd_vs_affinity(docking_results, ligand_name, out_path):
    """Scatter: RMSD (upper bound) vs affinity, one subplot per receptor.
    
    Bottom-left cluster = tight, reproducible best pose.
    Bottom-right spread = strong score but pose uncertainty.
    """
    successful = [d for d in docking_results if d.get("status") == "success"]
    if not successful:
        return None

    n = len(successful)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5), squeeze=False)

    for ax, d in zip(axes[0], successful):
        receptor = d["tag"].split("_")[0]
        affinities = [p["affinity_kcal_per_mol"] for p in d["poses"]]
        rmsds = [p["rmsd_upper_bound"] for p in d["poses"]]

        ax.scatter(rmsds, affinities, s=80, alpha=0.7,
                   color=_band_color(d["best_score"]),
                   edgecolor="black", linewidth=0.5)
        ax.scatter([rmsds[0]], [affinities[0]], s=250, marker="*",
                   color="red", edgecolor="black", linewidth=1.2, zorder=5,
                   label=f"Best: {affinities[0]:.2f} kcal/mol")

        ax.set_xlabel("RMSD upper bound (Å)", fontsize=11)
        ax.set_ylabel("Affinity (kcal/mol)", fontsize=11)
        ax.set_title(f"{receptor}", fontsize=12)
        ax.invert_yaxis()
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9, loc="upper right")

    fig.suptitle(f"Pose quality — {ligand_name}", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  📊 RMSD vs affinity saved: {out_path}")
    return out_path

# In[29]:


def visualization_agent(ligand_info: dict, docking_results: list,
                        model: str = DEFAULT_MODEL) -> dict:
    """Visualization agent: makes publication-ready figures from docking results."""
    print("==================================")
    print("📊 Visualization Agent")
    print("==================================")

    figures_dir = PROJECT_ROOT / "reports" / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    ligand_name = ligand_info.get("name", "ligand")
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    outputs = {}
    try:
        outputs["heatmap"] = tool_plot_heatmap(
            docking_results, ligand_name,
            figures_dir / f"heatmap_{stamp}.png")
        outputs["bar_chart"] = tool_plot_best_affinity_bar(
            docking_results, ligand_name,
            figures_dir / f"best_affinity_{stamp}.png")
        outputs["pose_distribution"] = tool_plot_pose_distribution(
            docking_results, ligand_name,
            figures_dir / f"pose_distribution_{stamp}.png")
        outputs["rmsd_vs_affinity"] = tool_plot_rmsd_vs_affinity(
            docking_results, ligand_name,
            figures_dir / f"rmsd_vs_affinity_{stamp}.png")
    except Exception as e:
        print(f"  ❌ Visualization failed: {e}")
        outputs["error"] = str(e)

    return outputs

# ## 4. Planner and Executor
# 
# This is the **same pattern as your C1M5 `planner_agent` + `executor_agent`**, adapted for docking. The planner produces a list of steps; the executor decides which specialized agent handles each step.
# 
# For docking, the step sequence is always the same (setup → ligand → receptor → dock → report), so the planner mostly just confirms the workflow — but keeping it as an agent lets you naturally extend the pipeline later (e.g. "also run MD simulation" or "also compare with haloperidol").

# In[30]:


def planner_agent(ligand_name: str, pdb_ids: list,
                  model: str = DEFAULT_MODEL) -> list:
    """Produce an ordered plan of steps for the docking workflow."""
    user_prompt = f"""
You are a planning agent for a molecular docking pipeline.

🧠 Available specialized agents:
- setup_agent     → creates folders, verifies obabel / vina / MGLTools, installs Python deps
- ligand_agent    → downloads a compound from PubChem by name, prepares SDF/PDB/PDBQT
- receptor_agent  → downloads PDB IDs from RCSB, cleans, computes grid center, writes PDBQT
- docking_agent   → runs AutoDock Vina (ligand × each receptor)
- visualization_agent → makes heatmap + bar chart + pose distribution + RMSD plots from docking results
- report_agent    → writes the final Markdown report

🎯 Ligand: "{ligand_name}"
🎯 Receptors: {pdb_ids}

Write the plan as a valid Python list of strings. Each string must be an atomic step
assignable to exactly one of the above agents. Do NOT include installation details,
visualization steps, MD simulation, or anything outside these agents' capabilities.
Return ONLY the Python list, no markdown fences, no explanation.
"""
    response = CLIENT.chat.completions.create(
        model=model,
        messages=[{"role": "user", "content": user_prompt}],
        temperature=1,
    )
    steps_str = response.choices[0].message.content.strip()
    steps_str = re.sub(r"^```(?:python)?\n?", "", steps_str).rstrip("`").strip()
    return ast.literal_eval(steps_str)


# In[31]:


# Shared blackboard — each agent writes its result here so later agents can read it
BLACKBOARD = {
    "ligand":    None,
    "receptors": None,
    "docking":   None,
    "figures":   None,   # ← NEW
    "report":    None,
}


def clean_json_block(raw: str) -> str:
    raw = raw.strip()
    if raw.startswith("```"):
        raw = re.sub(r"^```(?:json)?\n?", "", raw)
        raw = re.sub(r"\n?```$", "", raw)
    return raw.strip()


AGENT_REGISTRY = {
    "setup_agent":    setup_agent,
    "ligand_agent":   ligand_agent,
    "receptor_agent": receptor_agent,
    "docking_agent":  docking_agent,   # takes (ligand_info, receptor_infos)
    "visualization_agent": visualization_agent,
    "report_agent":   report_agent,    # takes (ligand_info, receptor_infos, docking_results)
}


def executor_agent(ligand_name: str, pdb_ids: list,
                   model: str = DEFAULT_MODEL) -> dict:
    """Orchestrate the full docking pipeline."""
    print("==================================")
    print("🚀 Executor — Multi-Agent Docking")
    print("==================================")

    plan = planner_agent(ligand_name, pdb_ids)
    print("\n📋 Plan:")
    for i, step in enumerate(plan, 1):
        print(f"  {i}. {step}")

    history = []

    for i, step in enumerate(plan, 1):
        # Route the step to the right agent
        route_prompt = f"""You are an execution manager for a molecular docking pipeline.

Available agents:
- setup_agent
- ligand_agent
- receptor_agent
- docking_agent
- visualization_agent
- report_agent

Given the instruction, return ONLY a JSON object:
{{"agent": "<one of the above>", "task": "<the cleaned task string>"}}

Instruction: "{step}"
"""
        r = CLIENT.chat.completions.create(
            model=model,
            messages=[{"role": "user", "content": route_prompt}],
            temperature=0,
        )
        info = json.loads(clean_json_block(r.choices[0].message.content))
        agent_name, task = info["agent"], info["task"]

        print(f"\n🛠️  Step {i}: {agent_name} — {task}")

        if agent_name == "setup_agent":
            out = setup_agent(task)
        elif agent_name == "ligand_agent":
            out = ligand_agent(f"{task}. Ligand: {ligand_name}")
            BLACKBOARD["ligand"] = out
        elif agent_name == "receptor_agent":
            out = receptor_agent(f"{task}. PDB IDs: {pdb_ids}")
            BLACKBOARD["receptors"] = out
        elif agent_name == "docking_agent":
            if not BLACKBOARD["ligand"] or not BLACKBOARD["receptors"]:
                out = "⚠️ Skipped: ligand or receptors not ready"
            else:
                out = docking_agent(BLACKBOARD["ligand"], BLACKBOARD["receptors"])
                BLACKBOARD["docking"] = out
        elif agent_name == "visualization_agent":                         # ← NEW block
            if not BLACKBOARD["docking"]:                                  # ← NEW
                out = "⚠️ Skipped: no docking results to visualize"       # ← NEW
            else:                                                          # ← NEW
                out = visualization_agent(BLACKBOARD["ligand"],            # ← NEW
                                          BLACKBOARD["docking"])           # ← NEW
                BLACKBOARD["figures"] = out                                # ← NEW
        elif agent_name == "report_agent":
            if not BLACKBOARD["docking"]:
                out = "⚠️ Skipped: no docking results to report"
            else:
                out = report_agent(BLACKBOARD["ligand"],
                                   BLACKBOARD["receptors"],
                                   BLACKBOARD["docking"])
                BLACKBOARD["report"] = out
        else:
            out = f"⚠️ Unknown agent: {agent_name}"

        history.append({"step": step, "agent": agent_name, "output": out})

    # -------------------------------------------------------------------------
    # SAFETY NET -- ensures visualization and report always run when there are
    # docking results, regardless of whether the planner included them.
    # Small local LLMs sometimes omit or duplicate steps; this fills the gaps.
    # -------------------------------------------------------------------------
    if BLACKBOARD.get("docking") and not BLACKBOARD.get("figures"):
        print("\n🛡️  Safety net: running visualization_agent (planner omitted it)")
        try:
            BLACKBOARD["figures"] = visualization_agent(
                BLACKBOARD["ligand"], BLACKBOARD["docking"])
            history.append({"step": "[safety-net] generate figures",
                            "agent": "visualization_agent",
                            "output": BLACKBOARD["figures"]})
        except Exception as e:
            print(f"   ❌ Safety-net visualization failed: {e}")

    if BLACKBOARD.get("docking") and not BLACKBOARD.get("report"):
        print("\n🛡️  Safety net: running report_agent (planner omitted it)")
        try:
            BLACKBOARD["report"] = report_agent(
                BLACKBOARD["ligand"],
                BLACKBOARD["receptors"],
                BLACKBOARD["docking"])
            history.append({"step": "[safety-net] generate report",
                            "agent": "report_agent",
                            "output": BLACKBOARD["report"]})
        except Exception as e:
            print(f"   ❌ Safety-net report failed: {e}")

    return {"plan": plan, "history": history, "blackboard": BLACKBOARD}

"""Verified source links, a citable methods catalogue and per-run records.

The catalogue covers available workflow components, including optional ones.
Project provenance distinguishes the tools and sources actually used in a run.
No request is made until a user explicitly resolves a drug name or CID.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import re
import threading
import time
from urllib.parse import quote, urlencode


_REFERENCES = [
    {"id": "vina12", "title": "AutoDock Vina 1.2", "kind": "Docking",
     "citation": "Eberhardt J, Santos-Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. J Chem Inf Model. 2021;61:3891–3898.",
     "url": "https://doi.org/10.1021/acs.jcim.1c00203"},
    {"id": "vina", "title": "Original AutoDock Vina method", "kind": "Docking",
     "citation": "Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. J Comput Chem. 2010;31:455–461.",
     "url": "https://doi.org/10.1002/jcc.21334"},
    {"id": "vina-guide", "title": "Docking setup and search-space guidance", "kind": "Docking",
     "citation": "AutoDock Vina documentation. Basic docking, search-space selection, preparation and output interpretation.",
     "url": "https://autodock-vina.readthedocs.io/en/latest/docking_basic.html"},
    {"id": "autodock-protocol", "title": "Protein–ligand docking protocol", "kind": "Docking",
     "citation": "Forli S et al. Computational protein–ligand docking and virtual drug screening with the AutoDock suite. Nat Protoc. 2016;11:905–919.",
     "url": "https://doi.org/10.1038/nprot.2016.051"},
    {"id": "meeko", "title": "Meeko docking-input preparation", "kind": "Molecule preparation",
     "citation": "Santos-Martins D et al. Meeko: molecule parameterization and software interoperability for docking and beyond. J Chem Inf Model. 2025;65:13045–13050. The installed software version is recorded with each run.",
     "url": "https://doi.org/10.1021/acs.jcim.5c02271"},
    {"id": "meeko-guide", "title": "Meeko input preparation documentation", "kind": "Molecule preparation",
     "citation": "Forli Laboratory. Meeko software documentation: molecule preparation, parameterisation and PDBQT input/output.",
     "url": "https://meeko.readthedocs.io/"},
    {"id": "rdkit", "title": "RDKit chemistry and descriptors", "kind": "Molecule preparation",
     "citation": "RDKit: Open-source cheminformatics. The installed software version is recorded with each run; version-specific releases are archived on Zenodo.",
     "url": "https://www.rdkit.org/"},
    {"id": "rdkit-book", "title": "RDKit standardisation, stereochemistry and force fields", "kind": "Molecule preparation",
     "citation": "RDKit documentation. The RDKit Book: chemical representation, standardisation, ETKDG conformers and MMFF/UFF molecular mechanics.",
     "url": "https://www.rdkit.org/docs/RDKit_Book.html"},
    {"id": "etkdg", "title": "ETKDG conformer generation", "kind": "Molecule preparation",
     "citation": "Riniker S, Landrum GA. Better Informed Distance Geometry: Using What We Know To Improve Conformation Generation. J Chem Inf Model. 2015;55:2562–2574.",
     "url": "https://doi.org/10.1021/acs.jcim.5b00654"},
    {"id": "etkdgv3", "title": "Small-ring and macrocycle conformers", "kind": "Molecule preparation",
     "citation": "Wang S, Witek J, Landrum GA, Riniker S. Improving Conformer Generation for Small Rings and Macrocycles Based on Distance Geometry and Experimental Torsional-Angle Preferences. J Chem Inf Model. 2020;60:2044–2058.",
     "url": "https://doi.org/10.1021/acs.jcim.0c00025"},
    {"id": "smiles", "title": "SMILES chemical-structure notation", "kind": "Molecule preparation",
     "citation": "Daylight Chemical Information Systems. SMILES theory manual: atoms, bonds, stereochemistry and chemical-structure notation.",
     "url": "https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html"},
    {"id": "pdbfixer", "title": "PDBFixer protein repair", "kind": "Protein preparation",
     "citation": "OpenMM project. PDBFixer software and manual: missing atoms, hydrogens, residues and structure preparation. The installed version is recorded with each run.",
     "url": "https://github.com/openmm/pdbfixer"},
    {"id": "openmm8", "title": "OpenMM protein-preparation runtime", "kind": "Protein preparation",
     "citation": "Eastman P et al. OpenMM 8: Molecular Dynamics Simulation with Machine Learning Potentials. J Phys Chem B. 2024;128:109–116.",
     "url": "https://doi.org/10.1021/acs.jpcb.3c06662"},
    {"id": "pdb", "title": "Protein Data Bank", "kind": "Protein data",
     "citation": "Berman HM et al. The Protein Data Bank. Nucleic Acids Res. 2000;28:235–242. Cite the selected structure's original publication as well as the database.",
     "url": "https://doi.org/10.1093/nar/28.1.235"},
    {"id": "pdb-downloads", "title": "Protein structure download formats", "kind": "Protein data",
     "citation": "RCSB Protein Data Bank. File Download Services. PDBx/mmCIF is the archive format; legacy PDB files are available only for supported entries.",
     "url": "https://www.rcsb.org/docs/programmatic-access/file-download-services"},
    {"id": "uniprot", "title": "UniProt target information", "kind": "Protein data",
     "citation": "The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. Nucleic Acids Res. 2025;53:D609–D617.",
     "url": "https://doi.org/10.1093/nar/gkae1010"},
    {"id": "alphafold", "title": "AlphaFold predicted structures", "kind": "Protein data",
     "citation": "Varadi M et al. AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences. Nucleic Acids Res. 2024;52:D368–D375.",
     "url": "https://doi.org/10.1093/nar/gkad1011"},
    {"id": "opentargets", "title": "Open Targets target evidence", "kind": "Target evidence",
     "citation": "Ochoa D et al. The next-generation Open Targets Platform: reimagined, redesigned, rebuilt. Nucleic Acids Res. 2023;51:D1353–D1359.",
     "url": "https://doi.org/10.1093/nar/gkac1046"},
    {"id": "chembl", "title": "ChEMBL bioactivity data", "kind": "Drug data",
     "citation": "Zdrazil B et al. The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods. Nucleic Acids Res. 2024;52:D1180–D1192.",
     "url": "https://doi.org/10.1093/nar/gkad1004"},
    {"id": "pubchem", "title": "PubChem drug structures", "kind": "Drug data",
     "citation": "Kim S et al. PubChem 2023 update. Nucleic Acids Res. 2023;51:D1373–D1380. Compound identifiers and retrieval URLs accompany resolved drugs.",
     "url": "https://doi.org/10.1093/nar/gkac956"},
    {"id": "pubchem-api", "title": "PubChem name/CID lookup and structure downloads", "kind": "Drug data",
     "citation": "National Center for Biotechnology Information. PubChem PUG REST documentation: compound properties, SMILES and SDF structures.",
     "url": "https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest"},
    {"id": "lipinski", "title": "Lipinski rule of five", "kind": "ADME descriptors",
     "citation": "Lipinski CA, Lombardo F, Dominy BW, Feeney PJ. Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings. Adv Drug Deliv Rev. 2001;46:3–26.",
     "url": "https://doi.org/10.1016/S0169-409X(00)00129-0"},
    {"id": "veber", "title": "Veber molecular-property rules", "kind": "ADME descriptors",
     "citation": "Veber DF et al. Molecular properties that influence the oral bioavailability of drug candidates. J Med Chem. 2002;45:2615–2623.",
     "url": "https://doi.org/10.1021/jm020017n"},
    {"id": "crippen", "title": "Calculated logP and molar refractivity", "kind": "ADME descriptors",
     "citation": "Wildman SA, Crippen GM. Prediction of physicochemical parameters by atomic contributions. J Chem Inf Comput Sci. 1999;39:868–873.",
     "url": "https://doi.org/10.1021/ci990307l"},
    {"id": "tpsa", "title": "Topological polar surface area", "kind": "ADME descriptors",
     "citation": "Ertl P, Rohde B, Selzer P. Fast calculation of molecular polar surface area as a sum of fragment-based contributions and its application to the prediction of drug transport properties. J Med Chem. 2000;43:3714–3717.",
     "url": "https://doi.org/10.1021/jm000942e"},
    {"id": "swissadme", "title": "SwissADME external predictions", "kind": "Optional external service",
     "citation": "Daina A, Michielin O, Zoete V. SwissADME: a free web tool to evaluate pharmacokinetics, drug-likeness and medicinal chemistry friendliness of small molecules. Sci Rep. 2017;7:42717.",
     "url": "https://doi.org/10.1038/srep42717"},
    {"id": "openmm7", "title": "OpenMM optional MD export", "kind": "Optional simulation",
     "citation": "Eastman P et al. OpenMM 7: Rapid development of high performance algorithms for molecular dynamics. PLoS Comput Biol. 2017;13:e1005659. MD packages are exported for separate execution.",
     "url": "https://doi.org/10.1371/journal.pcbi.1005659"},
    {"id": "openmm-guide", "title": "Optional MD force fields and simulation settings", "kind": "Optional simulation",
     "citation": "OpenMM User Guide. Force fields, water models, integration, constraints and simulation settings. Exported job.json records the selected models and parameters; cite those models if a job is actually executed.",
     "url": "https://docs.openmm.org/latest/userguide/"},
    {"id": "3dmol", "title": "3Dmol.js molecular viewer", "kind": "Visualisation",
     "citation": "Rego N, Koes D. 3Dmol.js: molecular visualization with WebGL. Bioinformatics. 2015;31:1322–1324.",
     "url": "https://doi.org/10.1093/bioinformatics/btu829"},
    {"id": "numpy", "title": "NumPy numerical calculations", "kind": "Software infrastructure",
     "citation": "Harris CR et al. Array programming with NumPy. Nature. 2020;585:357–362.",
     "url": "https://doi.org/10.1038/s41586-020-2649-2"},
    {"id": "scipy", "title": "SciPy scientific calculations", "kind": "Software infrastructure",
     "citation": "Virtanen P et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nat Methods. 2020;17:261–272.",
     "url": "https://doi.org/10.1038/s41592-019-0686-2"},
    {"id": "python", "title": "Python runtime", "kind": "Software infrastructure",
     "citation": "Python Software Foundation. Python programming language. The bundled interpreter version appears in the installation diagnostics.",
     "url": "https://www.python.org/"},
    {"id": "nicegui", "title": "NiceGUI application interface", "kind": "Software infrastructure",
     "citation": "Zauberzeug. NiceGUI: Python-based user interface framework.",
     "url": "https://nicegui.io/"},
    {"id": "pywebview", "title": "pywebview desktop window", "kind": "Software infrastructure",
     "citation": "pywebview contributors. pywebview: a native GUI window for web-based application interfaces.",
     "url": "https://pywebview.flowrl.com/"},
]


DOWNLOAD_LINKS = [
    {"title": "RCSB Protein Data Bank", "url": "https://www.rcsb.org/",
     "description": "Find experimental protein structures and their original publications."},
    {"title": "Protein structure batch downloads", "url": "https://www.rcsb.org/downloads",
     "description": "Download one or several structures in PDBx/mmCIF or supported legacy PDB format."},
    {"title": "AlphaFold Protein Structure Database", "url": "https://alphafold.ebi.ac.uk/",
     "description": "Find predicted models; inspect confidence and model limitations before docking."},
    {"title": "UniProt protein information", "url": "https://www.uniprot.org/",
     "description": "Find protein identifiers, sequences, functions and links to structures."},
    {"title": "PubChem drug and compound structures", "url": "https://pubchem.ncbi.nlm.nih.gov/",
     "description": "Find a compound by name or CID, then download SMILES or SDF. Confirm stereochemistry and salt form."},
    {"title": "ChEMBL molecules and bioactivity", "url": "https://www.ebi.ac.uk/chembl/",
     "description": "Find compound structures and experimentally reported bioactivity records."},
    {"title": "Open Targets evidence", "url": "https://platform.opentargets.org/",
     "description": "Review target evidence and existing drug-target information."},
    {"title": "SwissADME", "url": "https://www.swissadme.ch/",
     "description": "Optional external prediction service, opened and submitted manually."},
]


def references() -> list[dict[str, str]]:
    """Return independent catalogue records so UI state cannot alter citations."""
    return deepcopy(_REFERENCES)


def structure_links(pdb_id: str) -> list[dict[str, str]]:
    """Return exact-entry records and primary coordinate downloads."""
    identifier = str(pdb_id).strip().upper()
    if not re.fullmatch(r"(?:[1-9][A-Z0-9]{3}|PDB_[A-Z0-9]{8})", identifier):
        return []
    return [
        {"label": "Protein record and original publication", "title": "Protein record and original publication",
         "url": f"https://www.rcsb.org/structure/{identifier}"},
        {"label": "Download protein mmCIF", "title": "Download protein mmCIF",
         "url": f"https://files.rcsb.org/download/{identifier}.cif"},
        {"label": "Download protein PDB (when available)", "title": "Download protein PDB (when available)",
         "url": f"https://files.rcsb.org/download/{identifier}.pdb"},
    ]


def drug_links(name: str, smiles: str = "", cid: str | int = "") -> list[dict[str, str]]:
    """Point to the exact CID/SMILES when known, rather than guessing from a label."""
    clean_cid = str(cid).strip()
    if re.fullmatch(r"[1-9][0-9]*", clean_cid):
        record = f"https://pubchem.ncbi.nlm.nih.gov/compound/{clean_cid}"
        sdf = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{clean_cid}/SDF?"
    elif smiles:
        record = "https://pubchem.ncbi.nlm.nih.gov/#query=" + quote(smiles, safe="")
        sdf = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/SDF?" + urlencode({"smiles": smiles}) + "&"
    elif str(name).strip():
        term = quote(str(name).strip(), safe="")
        record = "https://pubchem.ncbi.nlm.nih.gov/#query=" + term
        sdf = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{term}/SDF?"
    else:
        return []
    return [
        {"label": "PubChem compound record", "title": "PubChem compound record", "url": record},
        {"label": "Download drug 2D SDF", "title": "Download drug 2D SDF", "url": sdf + "record_type=2d"},
        {"label": "Download drug 3D SDF (when available)", "title": "Download drug 3D SDF (when available)",
         "url": sdf + "record_type=3d"},
    ]


_PUBCHEM_LOCK = threading.Lock()
_LAST_PUBCHEM_REQUEST = 0.0


def resolve_drug(query: str, cache_dir: str | Path) -> dict:
    """Resolve one explicit name/CID, preserving available stereochemical SMILES.

    This uses one PubChem properties request. Calls are serialised and paced
    to at most four requests per second; errors never silently substitute a
    different compound or infer activity from its name.
    """
    from .web import Web

    term = str(query).strip()
    if not term or len(term) > 200 or any(ord(char) < 32 for char in term):
        raise ValueError("Enter a drug name or a positive PubChem CID.")
    numeric = re.fullmatch(r"(?:CID\s*[:#]?\s*)?([0-9]+)", term, re.IGNORECASE)
    namespace = "cid" if numeric else "name"
    value = numeric.group(1) if numeric else term
    if numeric and int(value) <= 0:
        raise ValueError("A PubChem CID must be a positive whole number.")
    endpoint = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
                + namespace + "/" + quote(value, safe="") + "/property/SMILES,Title/JSON")
    global _LAST_PUBCHEM_REQUEST
    with _PUBCHEM_LOCK:
        pause = 0.25 - (time.monotonic() - _LAST_PUBCHEM_REQUEST)
        if pause > 0:
            time.sleep(pause)
        _LAST_PUBCHEM_REQUEST = time.monotonic()
        response = Web(Path(cache_dir)).get(endpoint, service="PubChem", retries=0)
    if not response.ok:
        raise ValueError(response.problem + (" " + response.suggestion if response.suggestion else ""))
    payload = response.json()
    try:
        rows = payload["PropertyTable"]["Properties"]
        if not isinstance(rows, list) or not rows:
            raise KeyError("properties")
    except (KeyError, TypeError):
        raise ValueError("PubChem returned no usable compound structure. Use a specific CID or enter a SMILES structure.") from None
    if len(rows) != 1:
        raise ValueError("PubChem returned several compound records. Use a specific CID or an exact chemical name.")
    try:
        row = rows[0]
        if not isinstance(row, dict):
            raise TypeError("compound record")
        compound_id = int(row["CID"])
        # New PubChem naming first, then compatible older cached responses.
        structure = row.get("SMILES") or row.get("IsomericSMILES") or row.get("CanonicalSMILES") or row.get("ConnectivitySMILES")
        if compound_id <= 0 or not isinstance(structure, str) or not structure:
            raise KeyError("structure")
    except (KeyError, TypeError, IndexError, ValueError):
        raise ValueError("PubChem returned no usable compound structure. Use a specific CID or enter a SMILES structure.") from None
    from rdkit import Chem
    if Chem.MolFromSmiles(structure) is None:
        raise ValueError("The downloaded chemical structure could not be read. Use a verified SMILES structure.")
    name = str(row.get("Title") or term)
    return {"name": name, "smiles": structure, "cid": compound_id,
            "query": term, "source": "PubChem", "source_url": endpoint,
            "retrieved": response.retrieved, "from_cache": response.from_cache,
            "structure_links": drug_links(name, structure, compound_id)}


def _project_record(project_or_path):
    """Read owned folder handles, and leave borrowed Project handles open."""
    if isinstance(project_or_path, (str, Path)):
        from .project import Project
        owned = Project.open(Path(project_or_path))
        try:
            return owned.name, owned.provenance()
        finally:
            owned.close()
    return str(getattr(project_or_path, "name", "Docking target")), project_or_path.provenance()


def render_references(projects=None) -> str:
    """Export actual reproducibility records followed by the complete catalogue."""
    if projects is None:
        selected = []
    elif hasattr(projects, "provenance") or isinstance(projects, (str, Path)):
        selected = [projects]
    else:
        selected = list(projects)
    lines = ["IN-SILICO DOCKING STUDIO — METHODS, SOURCES AND REFERENCES", "",
             "Study records identify the tools and data actually used. The catalogue also includes optional workflow components.", ""]
    for project in selected:
        try:
            name, records = _project_record(project)
        except Exception as exc:
            lines.append("PROJECT: " + str(getattr(project, "name", project)))
            lines.append("Run record unavailable: " + str(exc))
            continue
        lines.extend(["PROJECT: " + str(name), ""])
        if not records:
            lines.append("No tool-use record has been saved for this project yet.")
        for record in records:
            lines.append(str(record.get("tool", "Tool")) + " | version: " + str(record.get("version", "not recorded")))
            for label, key in [("Stage", "stage_id"), ("Recorded", "timestamp"), ("Source release", "source_release")]:
                if record.get(key):
                    lines.append(label + ": " + str(record[key]))
            if record.get("seed") is not None:
                lines.append("Random seed: " + str(record["seed"]))
            if record.get("parameters"):
                parameters = record["parameters"]
                lines.append("Parameters: " + (parameters if isinstance(parameters, str) else json.dumps(parameters, ensure_ascii=False, sort_keys=True, default=str)))
            if record.get("citation"):
                lines.append("Recorded citation: " + str(record["citation"]))
            lines.append("")
    lines.extend(["PROTEIN AND DRUG DOWNLOAD SOURCES", ""])
    for link in DOWNLOAD_LINKS:
        lines.extend([link["title"] + ": " + link["url"], link["description"], ""])
    lines.extend(["REFERENCE CATALOGUE", ""])
    for index, reference in enumerate(_REFERENCES, 1):
        lines.extend([f"{index}. {reference['title']} [{reference['kind']}]",
                      reference["citation"], reference["url"], ""])
    lines.extend(["INTERPRETATION", "",
                  "Each target needs its own structure preparation and search box. Raw docking scores across different proteins are not a direct potency or selectivity comparison.",
                  "A reference-pose recovery threshold, a small comparison set or a descriptor range is a workflow check, not validation of affinity, clinical benefit or safety.",
                  "The built-in geometric pocket detector and descriptor-range comparison are local workflow methods; they are not external validated prediction services."])
    return "\n".join(lines) + "\n"

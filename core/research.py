"""Plain language study reports and portable, selected-pair research packages."""
from __future__ import annotations

import hashlib
import html
import io
import json
import math
from pathlib import Path
import zipfile

from . import contacts as contact_analysis
from . import viewer
from .runtime import resource_root

SWISSADME_URL = "https://www.swissadme.ch/index.php"


def _complex_pdb(receptor_pdb: str, ligand_pdb: str) -> str:
    """
    The protein and the docked pose in one file, for a 3D viewer.

    The export already carries receptor.pdbqt and recorded_pose.pdbqt,
    but PDBQT is a docking format that viewers read poorly, and neither
    file holds both molecules -- so a student had no way to look at the
    drug sitting in the pocket, which is the one result of a docking run
    worth looking at.
    """
    lines, serial = [], 0
    for line in receptor_pdb.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            serial += 1
            lines.append(f"{line[:6]}{serial:5d}{line[11:]}")
        elif line.startswith(("TER", "HELIX", "SHEET", "SSBOND", "LINK")):
            lines.append(line)
    lines.append("TER")
    # Vina names every ligand atom after its element alone, so a drug
    # arrives as "C, C, C, N, O, C ...". Inside one residue those are
    # duplicate atom names, and a PDB reader is entitled to drop the
    # repeats -- OpenMM silently discarded 26 of erlotinib's 29 atoms
    # from an earlier version of this file, which would have shown the
    # student a broken molecule. Numbering them per element fixes it.
    seen: dict[str, int] = {}
    for line in ligand_pdb.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        serial += 1
        element = (line[76:78].strip() or line[12:16].strip()[:1]).upper()
        seen[element] = seen.get(element, 0) + 1
        atom = f"{element.capitalize()}{seen[element]}"[:4]
        # A one- or two-letter element starts in column 14, not 13.
        name = atom.ljust(4) if len(atom) == 4 else f" {atom:<3}"
        lines.append(f"HETATM{serial:5d} {name}{line[16:]}"[:80])
    lines.append("END")
    return "\n".join(lines) + "\n"


def _chimerax_script(reference: str) -> str:
    """
    A ChimeraX session that opens on the binding site.

    Hydrogen bonds are drawn by ChimeraX's own `hbonds`, which checks
    donor-hydrogen-acceptor geometry. contacts.txt in the same folder
    lists distances only and says so. The two will not always agree,
    and where they differ ChimeraX is the one to believe.
    """
    resname = (reference.split()[0].upper() if reference.strip() else "")
    crystal = ""
    if resname and resname != "LIG":
        crystal = f"""
# The molecule from the crystal structure, in green, for comparison.
# Orange is where the computer predicted the drug would sit; green is
# where the experiment found it. If your re-docking check passed, the
# two should sit almost on top of each other.
show :{resname}
style :{resname} stick
color :{resname} green
"""
    return f"""# Binding site view -- In-Silico Docking Studio
#
# HOW TO USE THIS
#   1. Install ChimeraX (free for academic use) from
#      https://www.rbvi.ucsf.edu/chimerax/download.html
#   2. Start ChimeraX, then choose File > Open and pick this file.
#      Keep complex.pdb in the same folder.
#
# Every line below is one ChimeraX command, and each is explained. If
# any single line fails on your version, delete it and open the file
# again; the rest will still work.

# Load the protein with the docked drug already inside it.
open complex.pdb

# Start from a clean, plain white ribbon of the protein.
hide atoms
show cartoons
color white

# The docked drug, in orange.
show :LIG
style :LIG stick
color :LIG orange
{crystal}
# Every amino acid with an atom within 4.5 angstroms of the drug.
# These are the residues the drug is touching. Coloured by element:
# red is oxygen, blue is nitrogen, yellow is sulphur.
show :LIG :<4.5
style :LIG :<4.5 stick
color :LIG :<4.5 byhetero

# Hydrogen bonds between the drug and the protein, drawn as red dashes
# with the distance written on them. "restrict cross" means only bonds
# with one atom in the drug and one in the protein, not bonds inside
# the drug itself. ChimeraX checks the angles here, which the
# contacts.txt list in this folder does not.
hbonds :LIG restrict cross reveal true color red dashes 6 showDist true log true

# Name tags on the contacting residues. Delete this line if the view
# gets crowded.
label :LIG :<4.5 residues

# Point the camera at the drug.
view :LIG

# WHAT TO LOOK FOR
#   Which amino acids line the pocket, and are they the ones the
#   literature names for this target?
#   Is the drug buried, or sitting on the surface? A pose on the
#   outside of the protein is usually a wrong answer.
#   Are the red dashes going to sensible places -- an oxygen or
#   nitrogen at both ends?
#   Does the oily part of the drug face the oily part of the pocket?
#
# The H-bond list is written in the ChimeraX Log window (Tools >
# General > Log) because of "log true" above.
"""


def _chimera_script(reference: str) -> str:
    """The same view for the retired UCSF Chimera, whose syntax differs."""
    resname = (reference.split()[0].upper() if reference.strip() else "")
    crystal = (f"display :{resname}\ncolor green :{resname}\n"
               if resname and resname != "LIG" else "")
    return f"""# Binding site view -- for the older UCSF Chimera
#
# Use this only if you have the retired UCSF Chimera. If you have
# ChimeraX, which is the current program, use view_in_chimerax.cxc
# instead -- it is better explained and better tested.
#
# In Chimera: File > Open, then pick this file. Keep complex.pdb beside it.

open complex.pdb
~display
ribbon
color white
display :LIG
represent stick :LIG
color orange :LIG
{crystal}display :LIG zr<4.5
represent stick :LIG zr<4.5
findhbond makePseudobonds true color red lineWidth 2 log true
focus :LIG
"""


def contact_files(receptor_pdb: str, ligand_pdb: str, reference: str,
                  target: str = "", drug: str = "",
                  score: float | None = None) -> dict[str, bytes]:
    """complex.pdb, the two viewer scripts and the written contact list."""
    residues = contact_analysis.find(receptor_pdb, ligand_pdb,
                                     exclude=reference)
    return {
        "complex.pdb": _complex_pdb(receptor_pdb, ligand_pdb).encode("utf-8"),
        "view_in_chimerax.cxc": _chimerax_script(reference).encode("utf-8"),
        "view_in_chimera.com": _chimera_script(reference).encode("utf-8"),
        "contacts.txt": contact_analysis.report(
            residues, target, drug, score).encode("utf-8"),
    }


def preparation_review(project) -> list[dict]:
    r = viewer.values(project, "03_receptor")
    if not r:
        return []
    rows = [
        ("Hydrogen atoms added", r.get("hydrogens_added"), "These small atoms help describe chemical bonding and hydrogen bonds."),
        ("Missing atoms added", r.get("missing_atoms_added"), "Missing atoms in existing amino acids were modelled. Their positions are estimates."),
        ("Water molecules removed", r.get("waters_removed"), "These water molecules are absent from the prepared receptor. Important pocket waters may need a separate review."),
        ("Missing amino acids modelled", r.get("gaps_repaired"), "These amino acids were modelled rather than experimentally observed."),
        ("Protein segments left missing", len(r["gaps_left_open"]) if "gaps_left_open" in r else None, "These missing protein segments remain missing. Check whether they are near the drug pocket."),
    ]
    return [{"change": name, "value": "Not recorded" if value is None else value, "meaning": meaning}
            for name, value, meaning in rows]


def drug_review(drug: dict) -> str:
    changes = drug.get("changes", [])
    lines = list(changes) if changes else ["No salt or charge changes were recorded."]
    details = drug.get("preparation", {})
    if details:
        lines += [f"Hydrogen atoms added: {details['hydrogens_added']}. A 3D shape was built using {details['shape_method']}.",
                  "Shape relaxation using " + details["relaxation_method"] + (" finished within the set limit." if details["relaxation_converged"] else " reached its iteration limit. Review the starting geometry.")]
    else:
        lines.append("The saved file contains a prepared 3D structure. Detailed shape relaxation settings were not recorded by this earlier version.")
    lines += ["Docking can search other shapes.", "The prepared charge is a workflow choice. It may differ from the form present in the body."]
    return "\n".join(lines)


def _display(study, text) -> str:
    text = str(text)
    for drug in sorted(study.data["drugs"], key=lambda d: len(d["internal_name"]), reverse=True):
        text = text.replace(drug["internal_name"], drug["name"])
    return text


def report_text(study) -> str:
    from .study import results_rows, target_project
    from .workshop import docking_checks
    out = ["In Silico Docking Studio", "Docking report", study.data["name"], "Study status: " + study.data.get("status", "Not recorded"), "",
           "What this workflow does", "The software prepares protein and drug structures, then searches for ways each drug could fit in the selected protein pocket. Completed results are listed below.",
           "A docking score is a model estimate. A more negative number is favoured by this scoring method. It is not a measured binding strength or a recommended dose.", ""]
    for target in study.data["targets"]:
        out += ["Protein: " + target["name"], "PDB entry: " + (target.get("pdb_id") or "User supplied structure")]
        project = target_project(study, target["id"])
        try:
            out += ["", "Preparation changes"]
            for row in preparation_review(project):
                out.append(f"{row['change']}: {row['value']}. {row['meaning']}")
            receptor = viewer.values(project, "03_receptor")
            for label, key in (("Kept components", "kept"), ("Removed components", "removed")):
                if key in receptor:
                    out.append(label + ": " + (", ".join(str(x.get("label", "")) + " " + str(x.get("what", "")) for x in receptor[key]) or "None recorded"))
            for drug in viewer.values(project, "05_ligands").get("prepared", []):
                out += ["", "Prepared drug: " + _display(study, drug["name"]), _display(study, drug_review(drug)),
                        "Prepared SMILES: " + drug.get("smiles", "")]
            out += ["", "Docking results"]
            for row in results_rows(study):
                if row["target_id"] == target["id"]:
                    out.append(f"{row['drug']}: {row['status']}." + (f" Score: {row['score_kcal_mol']} kcal/mol." if row.get("score_kcal_mol") != "" else ""))
                    if row.get("error"):
                        out.append(row["error"])
            record = project.get_step("06_dock")
            dock = viewer.values(project, "06_dock")
            if dock:
                checks = docking_checks(dock, json.loads(record["control_json"] or "{}"))
                for name, check in checks.items():
                    label = {"baseline": "Score comparison", "re-docking": "Known pose recovery", "size": "Effect of molecular size"}.get(name, name)
                    state = {"pass": "Recorded criterion met", "fail": "Needs review", "warn": "Needs review", "not_applicable": "Not checked"}.get(check.get("verdict"), "Not checked")
                    out += ["", f"{label}: {state}", _display(study, check.get("plain_verdict", ""))]
                    if name == "baseline":
                        out.append("Comparison molecules are not proven nonbinders. A failed score comparison does not overturn experimental evidence for a known drug and target.")
                    if check.get("detail"):
                        out.append("Recorded numbers: " + _display(study, json.dumps(check["detail"], ensure_ascii=False)))
                out += ["", "Run settings: " + json.dumps(dock.get("settings", study.data.get("settings", {})), ensure_ascii=False),
                        "Search box: " + json.dumps(dock.get("box", {}), ensure_ascii=False)]
            else:
                out.append("No current completed docking result. Earlier results, if present, need updating.")
        finally:
            project.close()
        out.append("")
    out += ["Drug behaviour in the body", "The local drug descriptors are calculated chemical properties. SwissADME can provide additional ADME predictions.",
            "No SwissADME prediction was obtained automatically by this report. Paste the exported SMILES into " + SWISSADME_URL,
            "SwissADME predictions do not establish clinical pharmacokinetics, dose or safety.", "",
            "What to report", "Report the protein entry, prepared drug identity, search box, settings, scores, pose recovery and all comparison outcomes.",
            "Treat new docking results as hypotheses for further testing. A good score or a stable MD trajectory alone does not establish binding or activity."]
    return "\n".join(out) + "\n"


def report_html(study) -> bytes:
    return ("<!doctype html><html lang='en'><meta charset='utf-8'><title>Docking report</title>"
            "<style>body{max-width:900px;margin:40px auto;padding:20px;font:16px/1.65 system-ui;color:#21354b}pre{white-space:pre-wrap;overflow-wrap:anywhere}"
            "@media print{body{margin:0}}</style><pre>" + html.escape(report_text(study)) + "</pre></html>").encode("utf-8")


def swissadme_input(study) -> bytes:
    from .study import target_project
    prepared = {}
    # Use the prepared identities consistently; before preparation use explicit inputs.
    for target in study.data["targets"]:
        project = target_project(study, target["id"])
        try:
            for drug in viewer.values(project, "05_ligands").get("prepared", []):
                prepared.setdefault(drug["name"], drug["smiles"])
        finally:
            project.close()
    return ("\n".join(prepared.get(d["internal_name"], d["smiles"]) + " " + d["name"].replace("\n", " ").replace("\r", " ")
                       for d in study.data["drugs"]) + "\n").encode("utf-8")


def _file(project, filename: str) -> Path:
    path = Path(filename).resolve()
    if not filename or not path.is_relative_to(project.folder.resolve()) or not path.is_file():
        raise ValueError("A current structure file is missing or outside this project. Prepare and dock this pair again.")
    return path


def _opening(title: str, md: bool) -> str:
    """The first thing a student reads. Plain words, no jargon unexplained."""
    if md:
        return (
            "# " + title + "\n\n"
            "## What this notebook does\n\n"
            "Docking gave you a **still photograph**: one picture of the drug\n"
            "sitting in the pocket. It cannot tell you whether the drug would\n"
            "stay there. Molecular dynamics (MD) is the moving picture. It\n"
            "gives every atom of the protein, the drug and the surrounding\n"
            "water a push, works out the forces between them, and takes one\n"
            "tiny step forward in time. Then it does that again, two\n"
            "femtoseconds at a time, millions of times over.\n\n"
            "The question MD answers is: **does the pose survive?** A drug\n"
            "that slides out of the pocket in the first few nanoseconds was\n"
            "never really bound, whatever score docking gave it. A drug whose\n"
            "position settles and whose hydrogen bonds hold for most of the\n"
            "run is a real prediction.\n\n"
            "## Before you start\n\n"
            "Choose **Runtime > Change runtime type > GPU**. On a GPU this\n"
            "takes hours; on the ordinary processor it can take days for the\n"
            "same work.\n\n"
            "## What this notebook is not\n\n"
            "It is a starting protocol, not a finished study. Read `job.json`\n"
            "and `README.txt` to see every setting that was chosen for you.\n"
            "A simulation that runs without crashing is not evidence that the\n"
            "drug binds; it is evidence that the model did not fall apart.\n"
        )
    return (
        "# " + title + "\n\n"
        "## What this notebook does\n\n"
        "Molecular docking takes one protein and one drug molecule and asks\n"
        "two questions. **Where does the drug sit, and in what shape?** That\n"
        "answer is called the *pose*. **How tightly does it hold?** That\n"
        "answer is called the *score*, and it is printed in kilocalories per\n"
        "mole (kcal/mol), where a more negative number means the program\n"
        "thinks the binding is stronger.\n\n"
        "The program does not do an experiment. It takes a still photograph\n"
        "of the protein, treats your drug as a bendable object, tries a very\n"
        "large number of ways of fitting one into the other, and keeps the\n"
        "best few.\n\n"
        "## What is already decided for you\n\n"
        "Everything in the ZIP came from your own run on the desktop: the\n"
        "prepared protein, the prepared drug, the search box, and the search\n"
        "settings. This notebook re-runs that same job here so that you can\n"
        "see the code and change it. It does **not** download a new protein\n"
        "or a new drug.\n\n"
        "## One warning before you read any score\n\n"
        "A ten-fold change in binding strength is worth only about **1.4\n"
        "kcal/mol**, and this program's own error is about **2.5 kcal/mol**.\n"
        "So a score of -8.4 against -8.1 is not a difference; it is noise.\n"
        "Read `Docking_report.txt` in the ZIP for the checks that tell you\n"
        "how far your particular run can be trusted. Those checks are not\n"
        "repeated here.\n"
    )


def _install_text(md: bool) -> str:
    tools = ("OpenMM, which does the physics, along with the force fields\n"
             "that say how strongly each pair of atoms pulls on the other"
             if md else
             "AutoDock Vina, the docking program itself")
    return (
        "## Step 2. Install the scientific software\n\n"
        "Google Colab gives you a plain Python. It does not have " + tools
        + ".\n\nThis cell installs them. It does so inside a **separate**\n"
        "Python 3.11 of its own, in a folder under `/content`, instead of\n"
        "changing the Python that Colab itself runs on. That is deliberate:\n"
        "installing scientific packages over Colab's own Python usually ends\n"
        "with a demand to restart the runtime, and then you lose your\n"
        "uploaded files and start again.\n\n"
        "`micromamba` is a small program that fetches packages from\n"
        "conda-forge, a public store of scientific software.\n\n"
        "**This is the slow cell** -- usually three to six minutes. There is\n"
        "nothing to do but wait. A lot of text will scroll past; that is\n"
        "normal.\n"
    )


def _code_text(script: str, md: bool) -> str:
    """Show the student the actual science, instead of hiding it."""
    if md:
        return (
            "## Step 3. The code, and what it does\n\n"
            "The cell below prints `" + script + "` so you can read it before\n"
            "you run it. Follow it from the top. In order, it:\n\n"
            "1. **Checks the files have not changed** since the desktop wrote\n"
            "   them, by comparing a fingerprint of each one.\n"
            "2. **Builds the system**: your protein, the docked drug, a box of\n"
            "   water around them, and enough sodium and chloride to reach the\n"
            "   0.15 molar salt of body fluid.\n"
            "3. **Assigns a force field** -- AMBER ff14SB for the protein and\n"
            "   OpenFF for the drug. A force field is the set of rules saying\n"
            "   how stiff each bond is and how strongly each pair of atoms\n"
            "   attracts or repels.\n"
            "4. **Minimises the energy**, which pushes apart any atoms that\n"
            "   were placed too close together.\n"
            "5. **Warms the system up** to 310 K, body temperature, holding\n"
            "   the volume fixed, then lets the pressure settle.\n"
            "6. **Runs the production simulation**, saving the positions every\n"
            "   so often so you have a trajectory to look at afterwards.\n\n"
            "Change the code if you want to. That is the whole reason it is\n"
            "here rather than hidden inside the application.\n"
        )
    return (
        "## Step 3. The code, and what it does\n\n"
        "The cell below prints `" + script + "` so you can read it before you\n"
        "run it. This is the part worth your attention -- it is short, and\n"
        "the five lines in the middle are the whole of docking.\n\n"
        "```python\n"
        "vina = Vina(sf_name=\"vina\", cpu=0, seed=...)   # the program\n"
        "vina.set_receptor(\"receptor.pdbqt\")            # your protein\n"
        "vina.set_ligand_from_file(\"ligand.pdbqt\")      # your drug\n"
        "vina.compute_vina_maps(center=..., box_size=...)  # where to search\n"
        "vina.dock(exhaustiveness=..., n_poses=...)     # try, score, repeat\n"
        "```\n\n"
        "Taking those in turn:\n\n"
        "- **PDBQT** is the file format AutoDock uses. It is an ordinary PDB\n"
        "  coordinate file with two extra columns: the partial electrical\n"
        "  charge on each atom, and its atom type.\n"
        "- **compute_vina_maps** builds the search box. Nothing outside the\n"
        "  box is ever considered, which is what makes the search possible at\n"
        "  all -- searching a whole protein surface would take months.\n"
        "- **exhaustiveness** is how many independent attempts the program\n"
        "  makes, each starting from a random position. The default is 8.\n"
        "  Raise it to 32 if your molecule has many rotatable bonds.\n"
        "- **seed** is what makes the run repeatable. The search uses random\n"
        "  numbers, so the same seed gives the same answer twice. Change the\n"
        "  seed, run again, and compare: if the two scores differ by more\n"
        "  than about 0.2 kcal/mol, the search has not settled and you should\n"
        "  raise exhaustiveness.\n\n"
        "Try changing a number and running it again. That experiment teaches\n"
        "more about docking than reading about it does.\n"
    )


def _results_text(md: bool) -> str:
    if md:
        return (
            "## Step 4. Take your results away\n\n"
            "Colab deletes everything when the session closes, so download\n"
            "before you leave. The ZIP contains the trajectory, the final\n"
            "structure and a log of temperature, pressure and energy.\n\n"
            "**What to look at first.** Plot how far the drug has moved from\n"
            "its starting pose against time. If that distance rises and keeps\n"
            "rising, the drug left the pocket and the docking pose was wrong.\n"
            "If it rises for a nanosecond or two and then stays flat, the pose\n"
            "is holding, which is the result you were testing for.\n"
        )
    return (
        "## Step 4. Take your results away, and read them honestly\n\n"
        "Colab deletes everything when the session closes, so download now.\n\n"
        "**What is in the ZIP**\n\n"
        "- `docked.pdbqt` -- the predicted poses, best first. Several poses,\n"
        "  not one. If the top two score within 0.5 kcal/mol of each other but\n"
        "  sit in different places, the program is telling you it cannot\n"
        "  decide, and you should report that rather than pick one.\n"
        "- `results.json` -- the scores and the settings used. The\n"
        "  `energies_kcal_mol` list has one row per pose; the first number in\n"
        "  each row is the score.\n\n"
        "**Before you write anything down**\n\n"
        "You may fairly write: *the molecule is predicted to occupy this\n"
        "pocket, making contact with these residues*. Open `contacts.txt` and\n"
        "`view_in_chimerax.cxc` from the ZIP to see exactly which residues\n"
        "those are.\n\n"
        "You may **not** write that the score is an affinity, or that a\n"
        "compound scoring -8.4 binds more strongly than one scoring -8.1.\n"
        "With an error bar of 2.5 kcal/mol, a printed score of -8.0 is\n"
        "consistent with a real potency anywhere between about 40 nanomolar\n"
        "and 133 micromolar. That is a three-thousand-fold range.\n\n"
        "Docking is a way of deciding what to test. It is not a way of\n"
        "avoiding the test.\n"
    )


def notebook(title: str, script: str, md: bool = False) -> bytes:
    def cell(kind, source):
        return {"cell_type": kind, "id": hashlib.sha256(source.encode("utf-8")).hexdigest()[:12], "metadata": {}, "source": source.splitlines(keepends=True),
                **({"execution_count": None, "outputs": []} if kind == "code" else {})}
    packages = (["openmm>=8.1,<9", "openmmforcefields", "openff-toolkit-base", "openff-forcefields", "ambertools", "rdkit"]
                if md else ["vina=1.2.7"])
    install = ("import subprocess, sys, os\nfrom pathlib import Path\n"
               "tools = Path('/content/studio_tools')\nmicromamba = tools / 'micromamba'\nscience = tools / '" + ("md" if md else "docking") + "'\n"
               "environment = os.environ.copy()\nenvironment['MAMBA_ROOT_PREFIX'] = str(tools / 'root')\n"
               "subprocess.run([sys.executable, 'install_environment.py', '--platform', 'linux-64', '--destination', str(micromamba)], check=True)\n"
               "subprocess.run([str(micromamba), 'create', '--yes', '--prefix', str(science), '--override-channels', '--channel', 'conda-forge', 'python=3.11', "
               + ", ".join(repr(p) for p in packages) + "], env=environment, check=True)\n")
    # Every code cell gets a plain-English cell above it saying what it
    # does and why. An earlier version of this notebook explained only
    # how to operate it -- "run the cells in order" -- and the science
    # itself ran inside a subprocess the student never saw. They could
    # execute the whole thing and learn nothing from it.
    cells = [cell("markdown", _opening(title, md))]
    cells += [cell("markdown",
                   "## Step 1. Bring in your files\n\n"
                   "The ZIP you downloaded from Docking Studio holds everything this\n"
                   "notebook needs: your protein, your drug, the code, and the exact\n"
                   "settings used on the desktop. Nothing is downloaded from the\n"
                   "internet, so you are working on your own molecules and not on an\n"
                   "example.\n\n"
                   "Run the cell below. A **Choose Files** button appears. Pick the ZIP.\n"
                   "The cell unpacks it and prints `README.txt`, which names the\n"
                   "protein and drug so you can check you uploaded the right one.\n"),
              cell("code", "from google.colab import files\nfrom pathlib import Path\nimport zipfile, os, tempfile\nuploaded = files.upload()\narchives = [name for name in uploaded if name.endswith('.zip')]\nassert len(archives) == 1, 'Upload exactly one Docking Studio ZIP.'\nwork = Path(tempfile.mkdtemp(prefix='studio_job_', dir='/content')).resolve()\nwith zipfile.ZipFile(archives[0]) as archive:\n    for member in archive.infolist():\n        assert (work / member.filename).resolve().is_relative_to(work), 'Unsafe archive path'\n    archive.extractall(work)\nos.chdir(work)\nprint(Path('README.txt').read_text())\n")]
    if md:
        cells += [cell("markdown",
                       "### If your protein sits in a membrane\n\n"
                       "A receptor such as SV2A does not float in water. It sits in the\n"
                       "fatty membrane of the cell, and simulating it in water alone\n"
                       "would be simulating something that does not exist. Building a\n"
                       "membrane is a separate job, done in CHARMM-GUI or a similar\n"
                       "tool, and this notebook will stop until you supply the result.\n\n"
                       "For a protein that really is soluble, this cell does nothing.\n"),
                  cell("code", "import json\nif json.loads(Path('job.json').read_text())['system_type'] == 'membrane':\n    print('Upload the reviewed membrane_system.xml and matching membrane.pdb described in README.txt.')\n    files.upload()\n    assert all(Path(name).is_file() for name in ['membrane_system.xml', 'membrane.pdb']), 'Both prepared membrane files are required.'\n")]
    cells += [cell("markdown", _install_text(md)),
              cell("code", install),
              cell("markdown", _code_text(script, md)),
              cell("code", "print(Path('" + script + "').read_text())\n"),
              cell("markdown",
                   "### Now run it\n\n"
                   "This hands the script to the scientific environment built above.\n"
                   + ("Expect anything from several minutes to several hours, depending\n"
                      "on how many nanoseconds you asked for and whether you have a GPU.\n"
                      if md else
                      "It usually takes under a minute. The output appears below the cell\n"
                      "as it runs.\n")),
              cell("code", "import subprocess\nsubprocess.run([str(micromamba), 'run', '--prefix', str(science), 'python', '-u', '" + script + "'], env=environment, check=True)\n"),
              cell("markdown", _results_text(md)),
              cell("code", "import shutil\nfrom google.colab import files\nshutil.make_archive('/content/studio_results', 'zip', 'output')\nfiles.download('/content/studio_results.zip')\n")]
    return json.dumps({"nbformat": 4, "nbformat_minor": 5, "metadata": {"kernelspec": {"name": "python3", "display_name": "Python 3"}, "colab": {"provenance": []}}, "cells": cells}, indent=2).encode("utf-8")


def _zip(files: dict[str, bytes]) -> bytes:
    stream = io.BytesIO()
    with zipfile.ZipFile(stream, "w", zipfile.ZIP_DEFLATED) as archive:
        for name, data in files.items():
            archive.writestr(name, data)
    return stream.getvalue()


def export_pair(study, target_id: str, drug_id: str, kind: str = "docking", *, system_type: str = "water", nanoseconds: float = 1.0) -> bytes:
    """Export exactly the selected current pair; never fall back to a different drug."""
    from .study import target_project, results_rows
    target = next((t for t in study.data["targets"] if t["id"] == target_id), None)
    drug = next((d for d in study.data["drugs"] if d["id"] == drug_id), None)
    if target is None or drug is None:
        raise ValueError("Select a protein and a drug from this study.")
    pair = next(r for r in results_rows(study) if r["target_id"] == target_id and r["drug_id"] == drug_id)
    if pair["status"] != "completed":
        raise ValueError("This pair needs a current completed docking run before export.")
    if kind not in {"docking", "md"} or system_type not in {"water", "membrane"}:
        raise ValueError("Choose docking or MD and a supported system type.")
    if kind == "md":
        try:
            nanoseconds = float(nanoseconds)
        except (TypeError, ValueError):
            raise ValueError("Enter an MD duration greater than zero and no more than 1000 ns.") from None
        if not math.isfinite(nanoseconds) or not 0 < nanoseconds <= 1000:
            raise ValueError("Enter an MD duration greater than zero and no more than 1000 ns.")
    project = target_project(study, target_id)
    try:
        dock = viewer.values(project, "06_dock")
        result = next(r for r in dock["results"] if r["name"] == drug["internal_name"] and "score" in r)
        pose = _file(project, result["file"])
        ligand_input = _file(project, result.get("input_file") or str(pose.with_name(pose.name.removesuffix("_docked.pdbqt") + ".pdbqt")))
        receptor = _file(project, dock["receptor"])
        from stages.dock import DEFAULT_SEED, DEFAULT_POSES
        settings = dock.get("settings")
        if not settings:
            # Historical settings must come from this target's own run provenance.
            runs = [p for p in project.provenance() if p.get("stage_id") == "06_dock" and p.get("tool") == "AutoDock Vina"]
            if not runs:
                raise ValueError("This historical run has no recorded search settings. Run docking again before exporting code.")
            run = runs[-1]
            parameters = run.get("parameters", {})
            if isinstance(parameters, str):
                parameters = json.loads(parameters)
            settings = {"exhaustiveness": parameters.get("exhaustiveness"), "seed": run.get("seed", DEFAULT_SEED), "n_poses": DEFAULT_POSES, "vina_version": run.get("version", "unknown")}
        if not settings.get("exhaustiveness"):
            raise ValueError("Search effort was not recorded. Run docking again before exporting code.")
        job = {"target": target["name"], "pdb_id": target.get("pdb_id", ""), "drug": drug["name"], "smiles": result["smiles"],
               "recorded_score_kcal_mol": result["score"], "box": dock["box"], "settings": settings,
               "checks": dock.get("checks", {}), "redock_rmsd_angstrom": dock.get("redock_rmsd"),
               "purpose": "Selected pair rerun; comparison checks are preserved, not rerun."}
        files = {"receptor.pdbqt": receptor.read_bytes(), "ligand.pdbqt": ligand_input.read_bytes(), "recorded_pose.pdbqt": pose.read_bytes(),
                 "install_environment.py": (resource_root() / "assets/templates/bootstrap_micromamba.py").read_bytes(),
                 "Docking_report.txt": report_text(study).encode("utf-8")}
        # Something the student can actually look at. A missing or
        # unreadable prepared structure must not lose them the whole
        # export, so the viewer files are best-effort.
        try:
            from stages.md import pdbqt_pose_to_pdb
            receptor_step = viewer.values(project, "03_receptor")
            prepared_pdb = _file(project, receptor_step.get("prepared_file", "")
                                 ).read_text(errors="replace")
            files.update(contact_files(
                prepared_pdb, pdbqt_pose_to_pdb(pose.read_text(errors="replace")),
                receptor_step.get("ligand") or "",
                target["name"], drug["name"], result.get("score")))
        except (ValueError, OSError, KeyError):
            pass
        if kind == "docking":
            files["run_docking.py"] = (resource_root() / "assets/templates/run_docking.py").read_bytes()
            files["Run_Docking.ipynb"] = notebook("Dock your selected protein and drug", "run_docking.py")
            files["README.txt"] = ("Selected protein: " + target["name"] + "\nSelected drug: " + drug["name"] + "\n\n"
                "Upload Run_Docking.ipynb to Google Colab. Run its cells and upload this complete ZIP when asked.\n"
                "The code uses the exact prepared PDBQT inputs, recorded box and run settings. No new protein or drug is downloaded.\n"
                "It runs AutoDock Vina through its Python interface. It does not repeat preparation or the comparison controls.\n"
                "The desktop uses the Vina command line interface. Differences in processors, platforms and interfaces can affect results.\n"
                "The original checks remain in job.json and Docking_report.txt. Interpret scores with those checks.\n"
                "Local Linux use: install vina==1.2.7, then run python run_docking.py in this folder.\n"
                "\n"
                "TO LOOK AT THE RESULT\n"
                "complex.pdb holds the protein and the docked drug together in one file.\n"
                "Open view_in_chimerax.cxc in ChimeraX (free for academic use, from\n"
                "https://www.rbvi.ucsf.edu/chimerax/download.html) to see the drug in its\n"
                "pocket with the contacting amino acids labelled and the hydrogen bonds\n"
                "drawn. view_in_chimera.com does the same for the older UCSF Chimera.\n"
                "contacts.txt lists the same contacts as text, with distances, for anyone\n"
                "who has neither program installed.\n"
                "\n"
                "Contacts in contacts.txt are distances only. A hydrogen bond also needs the\n"
                "right angle, which ChimeraX checks and that list does not. Where the two\n"
                "disagree, believe ChimeraX.\n").encode("utf-8")
        else:
            files.update(_md_files(project, pose, drug, target, system_type, nanoseconds))
            job.update(system_type=system_type, nanoseconds=nanoseconds, timestep_fs=2.0, temperature_K=310,
                       salt_M=0.15, padding_nm=1.0, seed=2026, protein_forcefield="AMBER ff14SB", ligand_forcefield="OpenFF 2.2.1 with AM1-BCC",
                       purpose="MD starting protocol, not an executed simulation")
            job["water_model"] = "TIP3P"
            if system_type == "membrane":
                job.update(protein_forcefield="Defined in the supplied membrane System; record builder details",
                           ligand_forcefield="Defined in the supplied membrane System; record builder details",
                           water_model="Defined in the supplied membrane System", salt_M=None, padding_nm=None)
        job["input_sha256"] = {name: hashlib.sha256(data).hexdigest() for name, data in files.items() if name.endswith((".pdbqt", ".pdb", ".sdf"))}
        files["job.json"] = json.dumps(job, indent=2, ensure_ascii=False).encode("utf-8")
        return _zip(files)
    finally:
        project.close()


def _md_files(project, pose, drug, target, system_type, nanoseconds):
    from meeko import PDBQTMolecule, RDKitMolCreate
    from rdkit import Chem
    molecules = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(pose.read_text(encoding="utf-8"), poses_to_read=1, skip_typing=True))
    if len(molecules) != 1 or molecules[0] is None or molecules[0].GetNumConformers() != 1:
        raise ValueError("The docked drug cannot be reconstructed with its bond orders. Run docking again with this version.")
    mol = Chem.AddHs(molecules[0], addCoords=True)
    mol.SetProp("_Name", drug["name"])
    ligand_sdf = (Chem.MolToMolBlock(mol) + "\n$$$$\n").encode("utf-8")
    receptor = viewer.values(project, "03_receptor")
    prepared = _file(project, receptor["prepared_file"]).read_text(encoding="utf-8")
    raw = viewer.local_text(project, viewer.values(project, "02_structure").get("chosen", {}).get("local_file", ""))
    import re
    if system_type == "water" and re.search(r"\bSV2[ABC]\b|MEMBRANE PROTEIN", target["name"] + "\n" + raw[:1500], re.IGNORECASE):
        raise ValueError("This target is a membrane protein. Choose membrane mode and prepare a lipid system before MD.")
    reference = receptor.get("ligand", "").split()
    def keep(line):
        if not line.startswith(("ATOM  ", "HETATM", "TER", "SSBOND", "LINK")):
            return False
        if len(reference) >= 2 and line.startswith("HETATM"):
            label = line[17:20].strip() + " " + line[21:22] + line[22:27].strip()
            return label != " ".join(reference[:2])
        return True
    protein = "\n".join(line for line in prepared.splitlines() if keep(line)) + "\nEND\n"
    instruction = ("Membrane mode: build a lipid membrane system for this selected complex using a suitable membrane builder, such as CHARMM-GUI (https://www.charmm-gui.org/).\n"
                   "Export an OpenMM System as membrane_system.xml and its matching topology and coordinates as membrane.pdb into this folder.\n"
                   "The system must contain the selected drug as residue LIG, with correct parameters, solvent, ions and a periodic box.\n"
                   "Review protein identity, drug chemistry, membrane orientation and equilibration. The runner checks particle counts and the selected drug's element counts, but these checks cannot establish correct system construction.\n"
                   "The script stops until those files are provided. It does not simulate SV2A in water alone.\n" if system_type == "membrane" else
                   "Water mode is for a soluble protein with a standard, parameterisable drug. It is not suitable for membrane proteins such as SV2A.\n"
                   "Review the retained components and missing protein regions first. Unknown cofactors need appropriate parameters; the script stops on unsupported chemistry.\n")
    return {"protein.pdb": protein.encode("utf-8"), "ligand.sdf": ligand_sdf,
            "run_md.py": (resource_root() / "assets/templates/run_md.py").read_bytes(),
            "Run_MD.ipynb": notebook("MD for your selected protein and drug", "run_md.py", md=True),
            "README.txt": (f"Selected protein: {target['name']}\nSelected drug: {drug['name']}\nRequested production duration: {nanoseconds:g} ns\n\n" + instruction +
                "\nThis is downloadable code. No MD simulation has been run by the desktop app.\n"
                "ligand.sdf contains the docked coordinates and chemical bond orders recovered from Meeko mapping.\n"
                "Upload Run_MD.ipynb to Colab and follow its cells. Upload this ZIP first, then install its separate scientific environment. No runtime restart is needed.\n"
                "In membrane mode, the notebook also asks you to upload the two prepared membrane files before running.\n"
                "The default soluble protocol uses AMBER ff14SB, TIP3P water, OpenFF 2.2.1, AM1-BCC charges, 310 K and a 2 fs timestep.\n"
                "It minimises energy, runs 100 ps of NVT equilibration and 100 ps of NPT equilibration, then production. These durations require review for your system.\n"
                "Membrane mode preserves the supplied System's force fields and barostat. Its settings must be documented separately.\n"
                "Results and checkpoints are saved in output. Download them before a temporary Colab session ends.\n"
                "Resume requires the saved system, topology and checkpoint on a compatible OpenMM platform.\n"
                "A stable trajectory alone does not prove binding. Review trajectories, replicates, contacts and force field assumptions.\n"
                "References: https://docs.openmm.org/ ; https://github.com/openmm/openmmforcefields ; https://docs.openforcefield.org/\n").encode("utf-8")}

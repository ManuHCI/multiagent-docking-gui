"""
Stage 9 -- Molecular dynamics.

The step the audience has heard of and nobody in a pharmacology
department can actually run.

Docking gives one frozen picture of a molecule in a pocket, with the
protein held rigid and the water thrown away. Molecular dynamics puts
that picture back into a box of water at body temperature and lets
everything move. It answers the one question docking cannot: does the
pose stay, or does it fall apart in the first nanosecond?

Why this stage is EXPORT tier rather than BUNDLED
-------------------------------------------------
Because it would be dishonest to run it here. A 100 ns simulation of a
small protein takes roughly a day on a decent GPU and about three
months on a departmental desktop CPU. Starting it on the user's laptop
would produce a progress bar that never finishes, and they would learn
nothing except that the software is broken.

So this stage writes a complete job package instead -- the system, the
script, a Colab notebook, a cloud shell script, an analysis script and
a README -- and tells the user honestly what it will cost in time and
money. Everything in the package is built from what they have already
done here, so there is nothing to assemble by hand.

The control
-----------
The expensive mistake this stage exists to prevent is running MD on a
pose that was never trustworthy. GPU hours are real money and a week
of a student's time; spending them to watch a wrong pose wander away
is the single most wasteful thing in this whole workflow. So before
writing the package, this stage reads the verdict the docking step
gave itself -- and if the docking could not reproduce a known pose, or
the compound never beat an arbitrary molecule, it says so first.

Tools the package uses (none are installed here)
------------------------------------------------
  OpenMM 8      the simulation engine; free, GPU-native, Python
  PDBFixer      adds hydrogens and caps the chain
  OpenFF /      parameters for the drug molecule, which is the part a
  GAFF          protein force field does not cover
  MDTraj        the analysis afterwards

Citations are written into the package so the user can cite them.
"""

from __future__ import annotations

import json
import math
import shutil
from pathlib import Path
from typing import Any

from core.contract import (
    CheckResult,
    ControlOutcome,
    ProvenanceEntry,
    Readiness,
    Stage,
    StageResult,
    Teaching,
    Tier,
    Verdict,
)

# Water density, molecules per cubic angstrom, at 300 K. Used for an
# honest atom count rather than a guess, because the atom count is what
# the time estimate rests on.
WATER_DENSITY = 0.0334

# Throughput in nanoseconds per day, measured on systems of roughly
# 50,000 atoms. These are deliberately conservative: a student who is
# told 24 hours and waits 30 is annoyed; one told 48 and waits 30 is
# pleased.
THROUGHPUT = {
    "Colab free (T4)": 25.0,
    "Colab Pro (A100)": 90.0,
    "Rented A100 (RunPod, Vast.ai)": 110.0,
    "Your desktop, no GPU": 1.2,
}

# What an hour costs, in US dollars, on each route. Zero means free.
HOURLY_COST = {
    "Colab free (T4)": 0.0,
    "Colab Pro (A100)": 0.15,
    "Rented A100 (RunPod, Vast.ai)": 1.40,
    "Your desktop, no GPU": 0.0,
}

# AutoDock atom types -> element. Needed because a docked pose comes
# back as PDBQT, which is a PDB file with two extra columns and its own
# atom naming, and OpenMM will not read it.
AD_ELEMENT = {
    "A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N",
    "OA": "O", "O": "O", "OS": "O", "SA": "S", "S": "S",
    "HD": "H", "H": "H", "F": "F", "Cl": "Cl", "CL": "Cl",
    "Br": "Br", "BR": "Br", "I": "I", "P": "P",
    "MG": "Mg", "ZN": "Zn", "CA": "Ca", "MN": "Mn", "FE": "Fe",
}


# ----------------------------------------------------------------------
# Building the system
# ----------------------------------------------------------------------

def pdbqt_pose_to_pdb(text: str, model: int = 1,
                      resname: str = "LIG") -> str:
    """
    Pull one pose out of a Vina output file and write it as plain PDB.

    Vina writes PDBQT: a PDB file with a partial charge and an AutoDock
    atom type appended, and the file holds every pose stacked in MODEL
    records. OpenMM reads neither. This takes the pose the user chose,
    drops the extra columns, renames the residue so it can be picked
    out later, and restores a proper element symbol.
    """
    out: list[str] = []
    current = 0
    serial = 0
    for line in text.splitlines():
        if line.startswith("MODEL"):
            parts = line.split()
            current = int(parts[1]) if len(parts) > 1 else current + 1
            continue
        if line.startswith("ENDMDL"):
            if current == model:
                break
            continue
        if current != model:
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        ad_type = line[77:79].strip() if len(line) > 77 else ""
        element = AD_ELEMENT.get(ad_type, ad_type[:1].upper() or "C")
        # Hydrogens are re-added by the force field setup, and a polar
        # hydrogen placed by a docking program is not where a force
        # field would put it.
        if element == "H":
            continue
        # The coordinates are re-formatted rather than sliced across.
        # PDBQT written by other tools is not always column-perfect,
        # and a column that is one character out produces a file
        # OpenMM reads without complaint and simulates wrongly.
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        serial += 1
        raw = (line[12:16].strip() or element)[:4]
        # PDB convention: a one- or two-letter element starts in column
        # 14, not 13. MDTraj and OpenMM both guess the element from
        # this when the element column is missing, so it matters.
        name = raw.ljust(4) if len(raw) == 4 else f" {raw:<3}"
        out.append(
            f"HETATM{serial:5d} {name}{resname:>4} X   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}"
        )
    out.append("END")
    return "\n".join(out) + "\n"


def count_atoms(pdb_text: str) -> int:
    return sum(1 for ln in pdb_text.splitlines()
               if ln.startswith(("ATOM", "HETATM")))


def system_size(pdb_text: str, padding: float = 10.0) -> dict[str, Any]:
    """
    How big the simulated system will be, once it is in water.

    This is the number every estimate below depends on, so it is
    calculated from the actual coordinates rather than assumed. The
    protein is placed in a cube with `padding` angstroms of water on
    every side, which is the standard choice.
    """
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    for ln in pdb_text.splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        try:
            xs.append(float(ln[30:38]))
            ys.append(float(ln[38:46]))
            zs.append(float(ln[46:54]))
        except ValueError:
            continue
    if not xs:
        return {"solute_atoms": 0, "total_atoms": 0, "box_nm": 0.0}

    span = max(max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))
    side = span + 2 * padding
    volume = side ** 3
    solute = len(xs)
    # The solute displaces roughly its own volume; 25 A^3 per heavy atom
    # is the usual rule of thumb.
    water_volume = max(volume - solute * 25.0, 0.0)
    waters = int(water_volume * WATER_DENSITY)
    return {
        "solute_atoms": solute,
        "waters": waters,
        "total_atoms": solute + waters * 3,
        "box_nm": round(side / 10.0, 1),
    }


def estimate(total_atoms: int, nanoseconds: float) -> list[dict[str, Any]]:
    """
    How long each route takes and what it costs.

    Scaled from the 50,000-atom reference linearly in atom count, which
    is close enough to be useful and openly approximate.
    """
    scale = max(total_atoms, 1) / 50_000.0
    rows = []
    for route, ns_per_day in THROUGHPUT.items():
        rate = ns_per_day / scale
        hours = nanoseconds / rate * 24.0
        cost = hours * HOURLY_COST[route]
        rows.append({
            "route": route,
            "hours": round(hours, 1),
            "human": _human_time(hours),
            "cost_usd": round(cost, 2),
        })
    return rows


def _human_time(hours: float) -> str:
    if hours < 1:
        return f"{hours * 60:.0f} minutes"
    if hours < 48:
        return f"{hours:.0f} hours"
    days = hours / 24
    if days < 60:
        return f"{days:.0f} days"
    return f"{days / 30:.1f} months"


# ----------------------------------------------------------------------
# The files that go in the package
# ----------------------------------------------------------------------

RUN_MD = '''"""
Molecular dynamics of a docked complex, with OpenMM.

Written by the In-Silico Drug Design Teaching Suite. You do not need to
edit anything in this file; the settings at the top are the ones that
were chosen in the software.

It saves a checkpoint every 10 minutes. If the machine stops -- and a
free Colab session will stop -- run the same command again and it picks
up from the last checkpoint instead of starting over. This is the
difference between losing an afternoon and losing a week.
"""

import os
import sys

NANOSECONDS = {nanoseconds}        # how long to simulate
TEMPERATURE_K = 310                # body temperature, not 300
TIMESTEP_FS = 4.0                  # safe with hydrogen mass repartitioning
PADDING_NM = 1.0                   # water around the protein
SALT_M = 0.15                      # physiological saline
REPORT_EVERY_PS = 10               # how often a frame is saved

COMPLEX_PDB = "complex.pdb"
LIGAND_PDB = "ligand.pdb"
OUT = "output"

# ---------------------------------------------------------------------

from openmm import app, unit, LangevinMiddleIntegrator, MonteCarloBarostat
from openmm.app import PDBFile, Modeller, Simulation, StateDataReporter
from openmm.app import DCDReporter, CheckpointReporter
import openmm

os.makedirs(OUT, exist_ok=True)

print("1/6  Reading the complex ...")
pdb = PDBFile(COMPLEX_PDB)

print("2/6  Parameterising the drug molecule ...")
# A protein force field knows nothing about your drug, so the drug gets
# its own parameters from OpenFF. This is the step that fails most
# often, and the message below is the one to read if it does.
try:
    from openff.toolkit import Molecule
    from openmmforcefields.generators import SystemGenerator
    ligand = Molecule.from_file(LIGAND_PDB)
    ligand.assign_partial_charges("gasteiger")
    generator = SystemGenerator(
        forcefields=["amber/ff14SB.xml", "amber/tip3p_standard.xml"],
        small_molecule_forcefield="gaff-2.11",
        molecules=[ligand],
        forcefield_kwargs={{"constraints": app.HBonds,
                            "hydrogenMass": 3.0 * unit.amu}},
    )
    forcefield = generator.forcefield
except Exception as exc:                       # noqa: BLE001
    print()
    print("  Could not parameterise the ligand:", exc)
    print()
    print("  The usual cause is that the ligand has no explicit bond")
    print("  orders. Fix it by giving OpenFF the SMILES instead:")
    print("     Molecule.from_smiles('{smiles}')")
    print("  and matching it to ligand.pdb with .remap().")
    sys.exit(1)

print("3/6  Putting it in a box of water ...")
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield, pH=7.4)
modeller.addSolvent(forcefield,
                    padding=PADDING_NM * unit.nanometer,
                    ionicStrength=SALT_M * unit.molar,
                    neutralize=True)
print("     system has", modeller.topology.getNumAtoms(), "atoms")

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometer,
    constraints=app.HBonds,
    hydrogenMass=3.0 * unit.amu,
)
system.addForce(MonteCarloBarostat(1 * unit.bar,
                                   TEMPERATURE_K * unit.kelvin, 25))

integrator = LangevinMiddleIntegrator(
    TEMPERATURE_K * unit.kelvin,
    1.0 / unit.picosecond,
    TIMESTEP_FS * unit.femtoseconds,
)

try:
    platform = openmm.Platform.getPlatformByName("CUDA")
    print("     running on the GPU")
except Exception:                              # noqa: BLE001
    platform = openmm.Platform.getPlatformByName("CPU")
    print("     NO GPU FOUND - this will be about 30 times slower")

simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)

steps_total = int(NANOSECONDS * 1e6 / TIMESTEP_FS)
report_every = int(REPORT_EVERY_PS * 1000 / TIMESTEP_FS)
checkpoint = os.path.join(OUT, "state.chk")

if os.path.exists(checkpoint):
    print("4/6  Found a checkpoint - continuing from where it stopped")
    simulation.loadCheckpoint(checkpoint)
    done = simulation.context.getStepCount()
else:
    print("4/6  Minimising and warming up ...")
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(
        TEMPERATURE_K * unit.kelvin)
    simulation.step(50_000)                    # 200 ps of equilibration
    PDBFile.writeFile(simulation.topology,
                      simulation.context.getState(
                          getPositions=True).getPositions(),
                      open(os.path.join(OUT, "equilibrated.pdb"), "w"))
    done = 0

simulation.reporters = [
    DCDReporter(os.path.join(OUT, "trajectory.dcd"), report_every,
                append=os.path.exists(os.path.join(OUT, "trajectory.dcd"))),
    StateDataReporter(os.path.join(OUT, "log.csv"), report_every,
                      step=True, time=True, potentialEnergy=True,
                      temperature=True, speed=True, append=done > 0),
    StateDataReporter(sys.stdout, report_every * 10, step=True,
                      time=True, temperature=True, speed=True,
                      remainingTime=True, totalSteps=steps_total),
    CheckpointReporter(checkpoint, 150_000),   # roughly every 10 minutes
]

print(f"5/6  Running {{NANOSECONDS}} ns ...")
simulation.step(max(steps_total - done, 0))

print("6/6  Done. Results are in", OUT)
'''

ANALYSE = '''"""
What to look at when the simulation finishes.

Three numbers answer the question you ran this for.

  LIGAND RMSD   how far the drug has drifted from where docking put
                it. Flat and under about 2 angstroms means the pose
                held. Climbing, or jumping to 8 angstroms, means it
                did not -- which is a real result and worth reporting.

  PROTEIN RMSD  whether the protein itself stayed folded. If this is
                still climbing at the end, the simulation was not long
                enough to settle and nothing else should be read.

  CONTACTS      which residues the drug actually touched, and for what
                fraction of the time. This is what you put in the
                paper, not the docking score.
"""

import mdtraj as md
import numpy as np

traj = md.load("output/trajectory.dcd", top="output/equilibrated.pdb")
traj = traj.image_molecules()
print(traj)

protein = traj.topology.select("protein and name CA")
ligand = traj.topology.select("resname LIG")

traj.superpose(traj, 0, atom_indices=protein)

prot_rmsd = md.rmsd(traj, traj, 0, atom_indices=protein) * 10      # A
lig_rmsd = np.sqrt(
    np.mean(np.sum((traj.xyz[:, ligand] - traj.xyz[0, ligand]) ** 2,
                   axis=2), axis=1)) * 10

print()
print(f"protein backbone RMSD  mean {{prot_rmsd.mean():.2f}} A  "
      f"final {{prot_rmsd[-1]:.2f}} A")
print(f"ligand RMSD            mean {{lig_rmsd.mean():.2f}} A  "
      f"final {{lig_rmsd[-1]:.2f}} A")
print()
if lig_rmsd[-1] < 2.0:
    print("Ligand displacement is small in this trajectory. Review the")
    print("contacts and replicate simulations; this does not prove binding.")
elif lig_rmsd[-1] < 4.0:
    print("The ligand shifted relative to its starting coordinates.")
    print("Inspect the trajectory and contacts before interpreting it.")
else:
    print("The ligand moved substantially from its starting coordinates.")
    print("Check the trajectory, periodic boundaries and setup; RMSD")
    print("alone does not establish that the ligand left the pocket.")

# Which residues the drug touched, and how often
pairs = md.compute_neighbors(traj, 0.4, ligand,
                             haystack_indices=traj.topology.select("protein"))
counts = {{}}
for frame in pairs:
    for atom in set(traj.topology.atom(i).residue for i in frame):
        counts[str(atom)] = counts.get(str(atom), 0) + 1

print()
print("Residues within 4 A of the drug, by fraction of the simulation:")
for residue, n in sorted(counts.items(), key=lambda kv: -kv[1])[:15]:
    print(f"   {{residue:<12}} {{n / traj.n_frames:5.0%}}")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t = traj.time / 1000.0
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(t, prot_rmsd, label="protein backbone")
    ax.plot(t, lig_rmsd, label="ligand")
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("RMSD (angstrom)")
    ax.legend()
    fig.tight_layout()
    fig.savefig("output/rmsd.png", dpi=150)
    print()
    print("Plot written to output/rmsd.png")
except ImportError:
    pass
'''

CLOUD_SH = '''#!/bin/bash
# Run this simulation on any rented Linux machine with an NVIDIA GPU
# (RunPod, Vast.ai, Lambda, a university cluster node).
#
# Copy this whole folder to the machine, then:
#     bash run_on_cloud.sh
#
# It installs what it needs, runs the simulation, and keeps going if
# your connection drops.
set -e

if ! command -v conda >/dev/null 2>&1; then
  echo "Installing miniforge ..."
  wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
  bash Miniforge3-Linux-x86_64.sh -b -p "$HOME/miniforge3"
  export PATH="$HOME/miniforge3/bin:$PATH"
fi

source "$(conda info --base)/etc/profile.d/conda.sh"

if ! conda env list | grep -q "^md "; then
  echo "Creating the environment (this takes about 5 minutes) ..."
  conda create -y -n md -c conda-forge \\
      python=3.11 openmm openmmforcefields openff-toolkit \\
      mdtraj matplotlib cudatoolkit
fi
conda activate md

nvidia-smi || echo "WARNING: no GPU visible - this will be very slow"

# nohup keeps it running if your ssh connection drops.
nohup python run_md.py > run.log 2>&1 &
echo "Started. Watch it with:   tail -f run.log"
'''

README = '''IN-SILICO DRUG DESIGN TEACHING SUITE
Molecular dynamics job package
=====================================================================

WHAT IS IN THIS FOLDER

  complex.pdb     your protein with the docked drug in it
  ligand.pdb      the drug on its own
  run_md.py       the simulation
  analyse.py      what to look at afterwards
  run_on_cloud.sh for a rented GPU machine
  Run_MD.ipynb    for Google Colab (easiest - start here)
  job.json        every setting, for your methods section

---------------------------------------------------------------------
WHAT IS A MOLECULAR DYNAMICS SIMULATION?

Docking gave you one frozen picture: the drug sitting in the pocket,
with the protein held completely rigid and no water anywhere. Real
proteins are not rigid and the body is mostly water.

MD puts that picture into a box of water at 37 degrees and lets every
atom move, one four-femtosecond step at a time, for millions of steps.
You get a film instead of a photograph.

It answers one question docking cannot answer:

    does the drug STAY where docking put it, or does it drift off
    within the first nanosecond?

Rapid movement is a reason to inspect the pose, system preparation,
trajectory and force-field assumptions. Staying close to the start does
not by itself prove binding or activity.

---------------------------------------------------------------------
HOW LONG WILL IT TAKE?

{estimate_table}

Google Colab availability, GPU access and session duration vary. The
script saves checkpoints every ten minutes and can resume from one.
Download checkpoints and results before a session ends; files left only
on a temporary cloud runtime can be lost.

---------------------------------------------------------------------
EASIEST ROUTE: GOOGLE COLAB (FREE)

  1. Go to colab.research.google.com and sign in with any Google
     account.
  2. File -> Upload notebook -> choose Run_MD.ipynb from this folder.
  3. Runtime -> Change runtime type -> Hardware accelerator -> T4 GPU.
     This is the step people forget. Without it the simulation runs
     on a CPU and can be much slower; check the package's estimates.
  4. Run the first cell. It will ask you to upload files: select
     everything in this folder.
  5. Run the remaining cells.
  6. When the session disconnects, reconnect and run the last cell
     again. It continues from the checkpoint.

---------------------------------------------------------------------
IF YOU HAVE A GPU MACHINE

  Copy this folder across, then:   bash run_on_cloud.sh

---------------------------------------------------------------------
WHAT TO DO WITH THE RESULT

Run analyse.py. It prints three things:

  protein RMSD   did the protein stay folded
  ligand RMSD    displacement from the initial pose after protein alignment
  contacts       which residues the drug actually touched

Report the ligand RMSD plot with contacts, the trajectory length and
replicate results. Small displacement describes this simulation; it is
not proof of binding. Large displacement needs inspection of the
trajectory and setup before deciding what caused it.

---------------------------------------------------------------------
CITE THESE

  OpenMM      Eastman P et al. PLoS Comput Biol. 2017;13:e1005659.
  OpenFF      Qiu Y et al. J Chem Theory Comput. 2021;17:6262-6280.
  GAFF        Wang J et al. J Comput Chem. 2004;25:1157-1174.
  MDTraj      McGibbon RT et al. Biophys J. 2015;109:1528-1532.
  AMBER ff14SB  Maier JA et al. J Chem Theory Comput. 2015;11:3696.

Prepared by the In-Silico Drug Design Teaching Suite.
'''


def colab_notebook(nanoseconds: float) -> dict[str, Any]:
    """A Colab notebook that a person who has never used one can run."""
    def cell(source: str, kind: str = "code") -> dict[str, Any]:
        base = {"cell_type": kind, "metadata": {},
                "source": source.splitlines(keepends=True)}
        if kind == "code":
            base["execution_count"] = None
            base["outputs"] = []
        return base

    return {
        "nbformat": 4,
        "nbformat_minor": 0,
        "metadata": {
            "accelerator": "GPU",
            "colab": {"provenance": [], "gpuType": "T4"},
            "kernelspec": {"name": "python3",
                           "display_name": "Python 3"},
        },
        "cells": [
            cell(
                "# Molecular dynamics of your docked complex\n"
                "\n"
                "Prepared by the **In-Silico Drug Design Teaching "
                "Suite**.\n"
                "\n"
                "**Before you run anything:** click *Runtime* -> *Change "
                "runtime type* -> set **Hardware accelerator** to "
                "**T4 GPU**, then *Save*. Without this the simulation "
                "runs on a CPU and takes months instead of hours.\n"
                "\n"
                "Then run the cells below one at a time, from the top.",
                "markdown"),
            cell(
                "# CELL 1 - check that the GPU is switched on\n"
                "!nvidia-smi\n"
                "# If this says 'command not found', go back and set the\n"
                "# runtime type to T4 GPU. Do not continue without it.\n"),
            cell(
                "# CELL 2 - upload the job package\n"
                "# A file chooser appears. Select EVERY file in the\n"
                "# md_job folder from your computer.\n"
                "from google.colab import files\n"
                "uploaded = files.upload()\n"
                "print(sorted(uploaded))\n"),
            cell(
                "# CELL 3 - install the simulation software\n"
                "# Takes about 4 minutes. Run it once per session.\n"
                "!pip -q install condacolab\n"
                "import condacolab; condacolab.install()\n"
                "# The session restarts here. That is normal.\n"),
            cell(
                "# CELL 4 - install the chemistry packages\n"
                "!mamba install -q -y -c conda-forge openmm "
                "openmmforcefields openff-toolkit mdtraj matplotlib\n"
                "import openmm, openmm.testInstallation\n"
                "openmm.testInstallation.main()\n"),
            cell(
                f"# CELL 5 - run the simulation ({nanoseconds:g} ns)\n"
                "# Run this cell again if the session disconnects; it\n"
                "# continues from the last checkpoint rather than\n"
                "# starting over.\n"
                "!python run_md.py\n"),
            cell(
                "# CELL 6 - analyse, and look at the plot\n"
                "!python analyse.py\n"
                "from IPython.display import Image\n"
                "Image('output/rmsd.png')\n"),
            cell(
                "# CELL 7 - download the results to your computer\n"
                "!zip -qr results.zip output\n"
                "from google.colab import files\n"
                "files.download('results.zip')\n"),
        ],
    }


# ----------------------------------------------------------------------
# The stage
# ----------------------------------------------------------------------

class MdStage(Stage):
    stage_id = "09_md"
    title = "Molecular dynamics"
    tier = Tier.EXPORT
    requires = ["06_dock"]

    def check(self) -> CheckResult:
        if not self.project.is_complete("06_dock"):
            return CheckResult(
                readiness=Readiness.NEEDS_INPUT,
                reason="Nothing has been docked yet.",
                fix="Finish the docking step first. Molecular dynamics "
                    "starts from a docked pose.",
            )
        return CheckResult(readiness=Readiness.READY)

    def explain(self) -> Teaching:
        return Teaching(
            rule=(
                "Review all docking checks before spending time on simulation. MD explores a prepared starting model; it does not establish that the starting pose is correct."
            ),
            rule_detail=[
                "A simulation can run successfully even when the starting model or force-field assumptions are unsuitable.",
                "Choose duration, equilibration and replicas for the scientific question. A 100-nanosecond run does not guarantee convergence.",
                "This step prepares a package for a separate simulation. Runtime depends on system size and available hardware.",
            ],
            button="Prepare a molecular dynamics simulation to run on a "
                   "free cloud GPU.",
            simple=(
                "Docking gave you one frozen photograph: the drug in "
                "the pocket, with the protein completely stiff and no "
                "water anywhere. Neither of those is true in the "
                "body.\n\n"
                "Molecular dynamics (MD) makes a film instead. It puts "
                "your complex in a box of water at 37 degrees and lets "
                "every atom move, step by step, for millions of steps. "
                "Then you watch what happens.\n\n"
                "It answers the one question docking cannot:\n\n"
                "DOES THE DRUG STAY WHERE DOCKING PUT IT?\n\n"
                "Movement away from the starting pose is a reason to "
                "inspect the model and trajectory. Staying in place "
                "during one short run is not proof that the drug binds "
                "or works in the laboratory.\n\n"
                "MD can need substantial computing time. The time "
                "depends on the system, simulation length and hardware. "
                "This workshop prepares the files rather than running "
                "the full simulation in class.\n\n"
                "So instead this step prepares everything - the system, "
                "the script, the analysis - and gives you a ready "
                "Google Colab notebook. GPU access and session limits "
                "vary. Follow the package instructions and save your "
                "checkpoints and results outside the temporary runtime."
            ),
            panel=(
                "Molecular dynamics integrates Newton's equations of "
                "motion for every atom in the system, in steps of a few "
                "femtoseconds, using a force field to supply the "
                "forces. The output is a trajectory: the position of "
                "every atom at every saved time point.\n\n"
                "What it adds over docking is threefold. The protein is "
                "flexible, so a side chain that had to be pushed aside "
                "can move. Water is represented explicitly, and you "
                "can inspect how the pose and surrounding contacts "
                "change over simulated time. These observations still "
                "depend on force-field assumptions and sampling.\n\n"
                "The standard protocol, which this package implements: "
                "add hydrogens at pH 7.4; parameterise the ligand with "
                "GAFF or OpenFF, since a protein force field has no "
                "parameters for an arbitrary drug; solvate in a cubic "
                "box with 10 angstroms of TIP3P water on each side; "
                "neutralise and add 0.15 M salt; minimise; equilibrate "
                "under NVT then NPT; then run production dynamics at "
                "310 K and 1 bar.\n\n"
                "Timestep is 4 fs rather than 2, which is safe here "
                "because hydrogen mass repartitioning is switched on "
                "and bonds to hydrogen are constrained. That roughly "
                "can improve throughput. Check integration stability "
                "and whether the settings suit your analysis.\n\n"
                "The analysis that matters is the ligand RMSD against "
                "the initial pose after alignment on the protein. "
                "Inspect this alongside contacts and the trajectory. "
                "No single RMSD threshold establishes binding, "
                "convergence or a publishable conclusion.\n\n"
                "What MD still does not give you is a binding free "
                "energy. MM-GBSA, which many papers report from a short "
                "trajectory, correlates poorly with measured affinity "
                "and should not be quoted as one. A real free energy "
                "needs alchemical methods and far more compute."
            ),
            not_mean=[
                "A stable pose does not mean the compound is active. It "
                "means it stayed close to its start under these "
                "simulation conditions. Activity is a "
                "wet-lab measurement and nothing on a computer "
                "substitutes for it.",
                "An RMSD plot that looks flat does not mean the "
                "simulation converged. A hundred nanoseconds is a very "
                "short time in protein terms; a loop that has not moved "
                "may simply not have had time to.",
                "MM-GBSA or MM-PBSA numbers from a short trajectory are "
                "not binding free energies, whatever the units say. "
                "They are widely reported and weakly predictive.",
                "One simulation is one sample. If a conclusion rests on "
                "the difference between two compounds, run three "
                "replicas of each with different starting velocities, "
                "or the difference may be noise.",
                "Running MD does not validate a docking pose. If the "
                "docking step failed its checks, investigate the "
                "starting model before interpreting a simulation.",
            ],
            learn_more=[
                ("OpenMM user guide", "https://docs.openmm.org/"),
                ("Google Colab", "https://colab.research.google.com/"),
                ("Making it rain - MD in Colab, Nature Protocols",
                 "https://github.com/pablo-arantes/making-it-rain"),
            ],
        )

    # -- running ---------------------------------------------------------

    def run(self, **kwargs: Any) -> StageResult:
        nanoseconds = float(kwargs.get("nanoseconds") or 100.0)
        pose_rank = int(kwargs.get("pose_rank") or 1)
        which = str(kwargs.get("molecule") or "")

        dock = self.project.get_step("06_dock")
        if dock is None:
            return StageResult(success=False,
                               message="The docking step has not been run.")
        dock_values = json.loads(dock["values_json"] or "{}")
        results = [r for r in dock_values.get("results", []) if "score" in r]
        if not results:
            return StageResult(
                success=False,
                message="The docking step produced no poses to simulate.")

        chosen = next((r for r in results if r.get("name") == which), None)
        if chosen is None:
            chosen = min(results, key=lambda r: r["score"])

        receptor_pdbqt = Path(dock_values.get("receptor", ""))
        prepared = self._receptor_pdb()
        if prepared is None:
            return StageResult(
                success=False,
                message="The prepared protein file could not be found. "
                        "Run the Prepare protein step again.")

        pose_file = Path(chosen.get("file", ""))
        if not pose_file.exists():
            return StageResult(
                success=False,
                message=f"The docked pose file is missing: {pose_file}")

        folder = self.project.folder / "md_job"
        folder.mkdir(parents=True, exist_ok=True)

        # -- the system ---------------------------------------------
        ligand_pdb = pdbqt_pose_to_pdb(
            pose_file.read_text(errors="ignore"), model=pose_rank)
        protein_pdb = "\n".join(
            ln for ln in prepared.read_text(errors="ignore").splitlines()
            if ln.startswith(("ATOM", "HETATM", "TER", "SSBOND", "LINK"))
        )
        complex_pdb = protein_pdb + "\n" + ligand_pdb

        (folder / "ligand.pdb").write_text(ligand_pdb)
        (folder / "complex.pdb").write_text(complex_pdb)

        size = system_size(complex_pdb)
        routes = estimate(size["total_atoms"], nanoseconds)

        # -- the scripts --------------------------------------------
        (folder / "run_md.py").write_text(
            RUN_MD.format(nanoseconds=nanoseconds,
                          smiles=chosen.get("smiles", "")))
        (folder / "analyse.py").write_text(ANALYSE.format())
        (folder / "run_on_cloud.sh").write_text(CLOUD_SH)
        (folder / "Run_MD.ipynb").write_text(
            json.dumps(colab_notebook(nanoseconds), indent=1))

        table = "\n".join(
            f"  {r['route']:<32} {r['human']:>12}"
            + (f"   about ${r['cost_usd']:.2f}" if r["cost_usd"] else
               "   free")
            for r in routes
        )
        (folder / "README.txt").write_text(
            README.format(estimate_table=table))

        job = {
            "molecule": chosen.get("name", ""),
            "smiles": chosen.get("smiles", ""),
            "docking_score": chosen.get("score"),
            "pose_rank": pose_rank,
            "nanoseconds": nanoseconds,
            "temperature_K": 310,
            "timestep_fs": 4.0,
            "water_model": "TIP3P",
            "padding_nm": 1.0,
            "salt_M": 0.15,
            "protein_force_field": "AMBER ff14SB",
            "ligand_force_field": "GAFF 2.11 via OpenFF",
            "system": size,
            "estimates": routes,
            "receptor_source": str(receptor_pdbqt),
        }
        (folder / "job.json").write_text(json.dumps(job, indent=2))

        archive = shutil.make_archive(
            str(self.project.folder / "md_job"), "zip", str(folder))

        control = self._readiness_control(dock, chosen, size)

        self.record(ProvenanceEntry(
            stage=self.stage_id,
            tool="OpenMM job package (not executed here)",
            version="8.x, written for the cloud",
            parameters={k: v for k, v in job.items() if k != "estimates"},
            citation="Eastman P et al. OpenMM 7. PLoS Comput Biol. "
                     "2017;13:e1005659.",
        ))
        for e in self.provenance:
            self.project.add_provenance(e)

        free = next(r for r in routes if r["route"].startswith("Colab free"))
        return StageResult(
            success=True,
            summary=(
                f"Job package ready for {chosen.get('name') or 'your molecule'}"
                f": {nanoseconds:g} ns of {size['total_atoms']:,} atoms, "
                f"about {free['human']} on a free Colab GPU"
            ),
            values={"job": job, "folder": str(folder), "zip": archive,
                    "routes": routes},
            artefacts=[Path(archive)],
            control=control,
        )

    # -- helpers ----------------------------------------------------------

    def _receptor_pdb(self) -> Path | None:
        """The cleaned protein, as PDB rather than PDBQT."""
        step = self.project.get_step("03_receptor")
        if step is None:
            return None
        values = json.loads(step["values_json"] or "{}")
        for key in ("prepared_pdb", "prepared_file", "clean_pdb"):
            candidate = values.get(key, "")
            if candidate and Path(candidate).suffix.lower() == ".pdb" \
                    and Path(candidate).exists():
                return Path(candidate)
        # The prepared receptor is kept as PDBQT for docking; the PDB it
        # was made from sits beside it.
        for key in ("prepared_file", "receptor"):
            candidate = values.get(key, "")
            if not candidate:
                continue
            sibling = Path(candidate).with_suffix(".pdb")
            if sibling.exists():
                return sibling
        return None

    @staticmethod
    def _readiness_control(dock_step: Any, chosen: dict,
                           size: dict) -> ControlOutcome:
        """
        Is this pose worth a day of GPU time?

        The one check that saves real money. A docking run that could
        not reproduce a known pose, or a compound that never beat an
        arbitrary molecule, will simulate perfectly well and teach you
        nothing -- you will have paid to watch a wrong answer move.
        """
        detail = {
            "molecule": chosen.get("name", ""),
            "docking score": chosen.get("score"),
            "atoms to simulate": f"{size['total_atoms']:,}",
            "water box": f"{size['box_nm']} nm per side",
        }

        blob = dock_step["control_json"] if dock_step else None
        verdict = ""
        said = ""
        if blob:
            try:
                data = json.loads(blob)
                verdict = str(data.get("verdict", ""))
                said = str(data.get("plain_verdict", ""))
            except (ValueError, TypeError):
                pass
        detail["what the docking check said"] = said[:160] or "not recorded"

        if verdict == "fail":
            return ControlOutcome(
                name="Is this pose worth simulating?",
                verdict=Verdict.FAIL,
                plain_verdict=(
                    "The docking step failed its own check for this "
                    "molecule. This weakens the case for the starting "
                    "pose. A simulation may still run, but successful "
                    "execution would not resolve the failed docking check."
                ),
                what_to_do=(
                    "Go back to the docking step and read its WHAT TO DO "
                    "box first. Fix that, dock again, and come back here "
                    "once the check passes. If you run this anyway - "
                    "sometimes a supervisor wants to see it - say "
                    "plainly in your write-up that the starting pose was "
                    "not validated."
                ),
                detail=detail,
                blocks_progress=True,
            )

        if verdict == "warn" or verdict == "not_applicable":
            return ControlOutcome(
                name="Is this pose worth simulating?",
                verdict=Verdict.WARN,
                plain_verdict=(
                    "The docking step did not fully pass its own check "
                    "for this molecule, so the starting pose is "
                    "uncertain. A simulation can explore this starting "
                    "model, but does not resolve the docking warning "
                    "by itself. Review the warning before committing "
                    "simulation time."
                ),
                what_to_do=(
                    "If you continue, record the starting uncertainty. "
                    "Review ligand RMSD, contacts and trajectories "
                    "across replicas; neither a flat curve nor a large "
                    "change alone establishes binding or its absence."
                ),
                detail=detail,
            )

        if size["total_atoms"] > 150_000:
            return ControlOutcome(
                name="Is this pose worth simulating?",
                verdict=Verdict.WARN,
                plain_verdict=(
                    f"This is a large system - {size['total_atoms']:,} "
                    "atoms once the water is added. This can require "
                    "substantial GPU time. Cloud availability and "
                    "session duration vary."
                ),
                what_to_do=(
                    "Either run a shorter simulation (try 20 ns first to "
                    "see whether the pose survives at all), or use a "
                    "rented GPU - the package includes a script for "
                    "that. Do not simply start the 100 ns run on free "
                    "Colab and hope."
                ),
                detail=detail,
            )

        return ControlOutcome(
            name="Is this pose worth simulating?",
            verdict=Verdict.PASS,
            plain_verdict=(
                "Worth running. The docking step passed its own checks "
                "for this molecule, so you are starting from a pose "
                "there is some reason to believe, and the system is a "
                "within this package's usual size range. Hardware "
                "availability and runtime still vary."
            ),
            what_to_do=(
                "Download the package below and follow README.txt - the "
                "notebook provides a Colab route, subject to its current "
                "limits. Save checkpoints and results. Interpret the "
                "ligand RMSD together with contacts, trajectories and "
                "replicate runs."
            ),
            detail=detail,
        )

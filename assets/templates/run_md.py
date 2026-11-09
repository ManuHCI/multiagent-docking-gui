"""Reviewed starting protocol for one exported protein and docked drug.

Water mode uses explicit solvent. Membrane mode requires a separately built,
parameterised lipid system. Run from the folder containing job.json.
"""
from collections import Counter
import hashlib
import json
import math
from pathlib import Path


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smoke", action="store_true", help="Only 10 NVT, 10 NPT and 10 production steps to check execution. Not a scientific trajectory.")
    args = parser.parse_args()
    import numpy as np
    import openmm as mm
    from openmm import app, unit
    job = json.loads(Path("job.json").read_text(encoding="utf-8"))
    for name, digest in job["input_sha256"].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != digest:
            raise ValueError("An exported input changed: " + name)
    duration = float(job["nanoseconds"])
    if not math.isfinite(duration) or not 0 < duration <= 1000:
        raise ValueError("Duration must be greater than zero and no more than 1000 ns.")
    steps = 10 if args.smoke else int(duration * 1_000_000 / 2.0)
    equilibration_steps = 10 if args.smoke else 50_000
    out = Path("output")
    out.mkdir(exist_ok=True)
    signature = hashlib.sha256(Path("job.json").read_bytes() + Path(__file__).read_bytes() + str(args.smoke).encode()).hexdigest()
    checkpoint = out / "state.chk"
    restarting = checkpoint.exists()
    if restarting:
        if json.loads((out / "run.json").read_text())["signature"] != signature:
            raise ValueError("The job or code changed. Keep the old output and use a new output folder for a new protocol.")
        pdb = app.PDBFile(str(out / "topology.pdb"))
        system = mm.XmlSerializer.deserialize((out / "system.xml").read_text())
        topology, positions = pdb.topology, pdb.positions
    elif job["system_type"] == "membrane":
        if not Path("membrane_system.xml").is_file() or not Path("membrane.pdb").is_file():
            raise ValueError("Membrane proteins need lipids, water, ions and appropriate parameters. Read README.txt and provide membrane_system.xml and membrane.pdb. No simulation has started.")
        pdb = app.PDBFile("membrane.pdb")
        system = mm.XmlSerializer.deserialize(Path("membrane_system.xml").read_text())
        job["external_membrane_inputs_sha256"] = {name: hashlib.sha256(Path(name).read_bytes()).hexdigest()
                                                   for name in ("membrane_system.xml", "membrane.pdb")}
        topology, positions = pdb.topology, pdb.positions
        if system.getNumParticles() != topology.getNumAtoms():
            raise ValueError("The membrane System and topology have different particle counts.")
        from rdkit import Chem
        reference = next(iter(Chem.SDMolSupplier("ligand.sdf", removeHs=False)))
        if reference is None:
            raise ValueError("Cannot read the selected drug identity.")
        expected = Counter(a.GetSymbol() for a in reference.GetAtoms() if a.GetAtomicNum() > 1)
        actual = Counter(a.element.symbol for a in topology.atoms() if a.residue.name == "LIG" and a.element is not None and a.element.atomic_number > 1)
        if expected != actual:
            raise ValueError("The membrane system must contain exactly the selected drug as residue LIG. Review its chemistry and parameters.")
        if not system.usesPeriodicBoundaryConditions():
            raise ValueError("The membrane system must use periodic boundaries.")
        if topology.getPeriodicBoxVectors() is None:
            topology.setPeriodicBoxVectors(system.getDefaultPeriodicBoxVectors())
    elif job["system_type"] == "water":
        from rdkit import Chem
        from openff.toolkit import Molecule
        from openmmforcefields.generators import SMIRNOFFTemplateGenerator
        rdmol = next(iter(Chem.SDMolSupplier("ligand.sdf", removeHs=False)))
        if rdmol is None:
            raise ValueError("Cannot read the docked ligand with its chemical bond information.")
        ligand = Molecule.from_rdkit(rdmol, allow_undefined_stereo=False)
        ligand.assign_partial_charges("am1bcc")
        template = SMIRNOFFTemplateGenerator(molecules=[ligand], forcefield="openff-2.2.1")
        forcefield = app.ForceField("amber14/protein.ff14SB.xml", "amber14/tip3p.xml")
        forcefield.registerTemplateGenerator(template.generator)
        pdb = app.PDBFile("protein.pdb")
        modeller = app.Modeller(pdb.topology, pdb.positions)
        # Unknown retained cofactors deliberately fail unless parameters exist.
        modeller.addHydrogens(forcefield, pH=7.4)
        ligand_topology = ligand.to_topology().to_openmm()
        for residue in ligand_topology.residues():
            residue.name = "LIG"
        modeller.add(ligand_topology, np.array(rdmol.GetConformer().GetPositions()) * unit.angstrom)
        modeller.addSolvent(forcefield, model="tip3p", padding=1.0 * unit.nanometer,
                            ionicStrength=0.15 * unit.molar, neutralize=True)
        topology, positions = modeller.topology, modeller.positions
        system = forcefield.createSystem(topology, nonbondedMethod=app.PME,
                                          nonbondedCutoff=1.0 * unit.nanometer, constraints=app.HBonds)
        system.addForce(mm.MonteCarloBarostat(1 * unit.bar, 310 * unit.kelvin, 25))
    else:
        raise ValueError("Choose water or membrane mode and review the preparation.")

    barostats = [f for f in system.getForces() if hasattr(f, "getFrequency") and "Barostat" in type(f).__name__]
    if len(barostats) != 1:
        raise ValueError("The reviewed system must contain exactly one appropriate pressure barostat.")
    barostat = barostats[0]
    frequency = barostat.getFrequency()
    if not restarting:
        barostat.setDefaultTemperature(310 * unit.kelvin)
        barostat.setFrequency(0)
    integrator = mm.LangevinMiddleIntegrator(310 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    integrator.setRandomNumberSeed(int(job["seed"]))
    try:
        platform = mm.Platform.getPlatformByName("CUDA")
        properties = {"Precision": "mixed"}
        simulation = app.Simulation(topology, system, integrator, platform, properties)
    except Exception as error:
        print("CUDA could not start; trying the CPU:", error)
        platform = mm.Platform.getPlatformByName("CPU")
        properties = {}
        integrator = mm.LangevinMiddleIntegrator(310 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
        integrator.setRandomNumberSeed(int(job["seed"]))
        simulation = app.Simulation(topology, system, integrator, platform, properties)
    simulation.context.setPositions(positions)
    if restarting:
        simulation.loadCheckpoint(str(checkpoint))
    else:
        print(f"Minimising energy, then running {equilibration_steps} NVT and {equilibration_steps} NPT steps.")
        simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(310 * unit.kelvin, int(job["seed"]))
        simulation.step(equilibration_steps)
        barostat.setFrequency(frequency or 25)
        simulation.context.reinitialize(preserveState=True)
        simulation.step(equilibration_steps)
        simulation.currentStep = 0
        simulation.context.setTime(0 * unit.picoseconds)
        state = simulation.context.getState(getPositions=True)
        with (out / "topology.pdb").open("w") as stream:
            app.PDBFile.writeFile(topology, state.getPositions(), stream)
        (out / "system.xml").write_text(mm.XmlSerializer.serialize(system))
        (out / "run.json").write_text(json.dumps({"signature": signature, "openmm_version": mm.version.version,
              "platform": platform.getName(), "job": job, "barostat": type(barostat).__name__,
              "smoke_test_only": args.smoke, "production_steps": steps, "equilibration_steps_per_ensemble": equilibration_steps}, indent=2))
        simulation.saveCheckpoint(str(checkpoint))
    remaining = max(0, steps - simulation.currentStep)
    if remaining:
        simulation.reporters.append(app.DCDReporter(str(out / "trajectory.dcd"), 500, append=(out / "trajectory.dcd").exists()))
        simulation.reporters.append(app.StateDataReporter(str(out / "energy.csv"), 500, step=True, time=True,
                          potentialEnergy=True, temperature=True, density=True, append=(out / "energy.csv").exists()))
        simulation.reporters.append(app.CheckpointReporter(str(checkpoint), 5000))
        print(f"Running {remaining} production steps on {platform.getName()}.")
        simulation.step(remaining)
    simulation.saveCheckpoint(str(checkpoint))
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    if not math.isfinite(energy):
        raise RuntimeError("The simulation produced nonfinite energy. Review the setup and trajectory.")
    with (out / "final.pdb").open("w") as stream:
        app.PDBFile.writeFile(topology, state.getPositions(), stream)
    print("Finished. Download output, including the checkpoint, system and topology. Review the trajectory and replicates before interpreting stability.")


if __name__ == "__main__":
    main()

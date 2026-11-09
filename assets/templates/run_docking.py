"""Rerun the selected prepared pair using recorded settings and Vina 1.2.7."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import re


def main():
    from vina import Vina
    job = json.loads(Path("job.json").read_text(encoding="utf-8"))
    for name, digest in job["input_sha256"].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != digest:
            raise ValueError("An exported input changed: " + name)
    settings, box = job["settings"], job["box"]
    recorded_version = settings.get("vina_version", "")
    installed_version = importlib.metadata.version("vina")
    semantic_version = re.search(r"\b\d+\.\d+\.\d+\b", recorded_version)
    if semantic_version and installed_version != semantic_version.group():
        raise ValueError(f"Recorded Vina version: {recorded_version}. Installed: {installed_version}. Install the matching version before rerunning.")
    notice = ("The desktop command line build reports " + (recorded_version or "an unknown version")
              + ". This rerun uses the Vina " + installed_version + " Python interface. Results can differ.")
    print(notice)
    output = Path("output")
    output.mkdir(exist_ok=True)
    vina = Vina(sf_name="vina", cpu=0, seed=int(settings["seed"]))
    vina.set_receptor("receptor.pdbqt")
    vina.set_ligand_from_file("ligand.pdbqt")
    vina.compute_vina_maps(center=[box["center_" + a] for a in "xyz"],
                           box_size=[box["size_" + a] for a in "xyz"])
    vina.dock(exhaustiveness=int(settings["exhaustiveness"]), n_poses=int(settings["n_poses"]))
    vina.write_poses(str(output / "docked.pdbqt"), n_poses=int(settings["n_poses"]), overwrite=True)
    energies = vina.energies(n_poses=int(settings["n_poses"]))
    results = {"target": job["target"], "drug": job["drug"], "vina_version": installed_version,
               "settings": settings, "box": box, "energies_kcal_mol": energies.tolist(), "rerun_notice": notice,
               "scope": "Selected pair rerun. This does not rerun or validate the comparison checks."}
    (output / "results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()

"""Run one fresh docking search using only a distributed build's resources.

This verifies the packaged executable and a real search, not docking accuracy.
Reference resources are opened read-only; outputs are temporary.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sqlite3
import subprocess
import tempfile
import time
import traceback


def verify(resources: Path, destination: Path) -> int:
    started = time.monotonic()
    report = {"passed": False, "mode": "fresh docking search", "seed": 20261107,
              "exhaustiveness": 1, "cpu": 1, "requested_poses": 1,
              "interpretation": "Executable/search verification; not evidence of potency or accuracy."}
    try:
        resources = resources.resolve()
        destination = destination.resolve()
        destination.parent.mkdir(parents=True, exist_ok=True)
        case = resources / "assets" / "workshop" / "egfr"
        receptor = case / "steps" / "06_dock" / "receptor.pdbqt"
        ligand = case / "steps" / "06_dock" / "Erlotinib.pdbqt"
        executable = resources / "vina.exe"
        for path in (executable, receptor, ligand, case / "project.db"):
            if not path.is_file():
                raise FileNotFoundError(f"The distributed resource is missing: {path}")
        with sqlite3.connect((case / "project.db").as_uri() + "?mode=ro", uri=True) as connection:
            row = connection.execute("SELECT values_json FROM steps WHERE stage_id=?",
                                     ("04_site",)).fetchone()
        if row is None:
            raise ValueError("The saved binding-site box is missing.")
        box = json.loads(row[0])["box"]
        report["input_sha256"] = {
            "vina.exe": hashlib.sha256(executable.read_bytes()).hexdigest(),
            "receptor.pdbqt": hashlib.sha256(receptor.read_bytes()).hexdigest(),
            "Erlotinib.pdbqt": hashlib.sha256(ligand.read_bytes()).hexdigest(),
        }
        with tempfile.TemporaryDirectory(prefix="live_docking_", dir=destination.parent) as temporary:
            output = Path(temporary) / "erlotinib_fresh_docked.pdbqt"
            command = [str(executable), "--receptor", str(receptor), "--ligand", str(ligand),
                       "--out", str(output), "--exhaustiveness", "1", "--num_modes", "1",
                       "--cpu", "1", "--seed", "20261107"]
            for key in ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z"):
                command.extend(["--" + key, str(box[key])])
            search = subprocess.run(command, capture_output=True, text=True, timeout=60,
                                    creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            if search.returncode != 0 or not output.is_file():
                raise RuntimeError((search.stdout + search.stderr)[-2500:] or
                                   f"Vina exited with {search.returncode} and no poses.")
            text = output.read_text(encoding="utf-8")
            scores = [float(value) for value in
                      re.findall(r"REMARK VINA RESULT:\s*([-+\d.]+)", text)]
            if len(scores) != 1 or not math.isfinite(scores[0]):
                raise ValueError("The fresh docking output did not contain one finite-scored pose.")
            if not any(line.startswith(("ATOM", "HETATM")) for line in text.splitlines()):
                raise ValueError("The fresh docking pose contains no atoms.")
            report.update(passed=True, score_kcal_mol=scores[0], poses=len(scores),
                          output_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                          output_bytes=output.stat().st_size)
    except Exception as error:
        report.update(error=f"{type(error).__name__}: {error}", traceback=traceback.format_exc())
    report["seconds"] = round(time.monotonic() - started, 3)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in report.items() if key not in {"traceback", "input_sha256"}}))
    return 0 if report["passed"] else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resources", type=Path, required=True,
                        help="Application resource directory, normally the installed _internal folder.")
    parser.add_argument("--output", type=Path, required=True, help="JSON acceptance report.")
    arguments = parser.parse_args()
    return verify(arguments.resources, arguments.output)


if __name__ == "__main__":
    raise SystemExit(main())

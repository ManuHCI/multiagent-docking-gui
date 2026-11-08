"""Real offline acceptance of the multi-target/multi-drug engine."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import tempfile
import time
import traceback
import zipfile

from .runtime import resource_root


def run_study_checks(output_path: str) -> int:
    started = time.monotonic()
    report = {"passed": False, "mode": "real offline two-target by two-drug docking",
              "interpretation": "Execution and result-identity checks, not validation of affinity or potency."}
    try:
        from .study import (DrugSpec, TargetSpec, create_study, export_csv, export_package,
                            open_study, results_rows, run_docking, run_preparation, target_project)
        examples = resource_root() / "assets" / "workshop"
        hashes = {case: hashlib.sha256((examples / case / "project.db").read_bytes()).hexdigest()
                  for case in ("egfr", "sv2a")}
        parent = os.environ.get("INSILICO_DATA_DIR") or None
        if parent:
            Path(parent).mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="batch_acceptance_", dir=parent) as temporary:
            study = create_study(Path(temporary), "Two-target acceptance", [
                TargetSpec("EGFR acceptance", pdb_id="1M17", reference_project=str(examples / "egfr")),
                TargetSpec("SV2A acceptance", pdb_id="9OKG", reference_project=str(examples / "sv2a")),
            ], [DrugSpec("Ethanol acceptance", "CCO"), DrugSpec("Acetamide acceptance", "CC(N)=O")])
            run_preparation(study)
            if any(target["status"] != "prepared" for target in study.data["targets"]):
                raise AssertionError("Target preparation did not finish: " + json.dumps(study.data["targets"]))
            boxes = []
            for target in study.data["targets"]:
                project = target_project(study, target["id"])
                try:
                    boxes.append(json.loads(project.get_step("04_site")["values_json"])["box"])
                finally:
                    project.close()
            if boxes[0] == boxes[1]:
                raise AssertionError("Different proteins unexpectedly share an identical search box")
            run_docking(study, exhaustiveness=1, n_decoys=1)
            rows = results_rows(open_study(study.folder))
            if len(rows) != 4 or len({(row["target_id"], row["drug_id"]) for row in rows}) != 4:
                raise AssertionError("The two-by-two pairing matrix is incomplete")
            if any(row["status"] != "completed" or not math.isfinite(float(row["score_kcal_mol"])) for row in rows):
                raise AssertionError("A docking pair failed: " + json.dumps(rows))
            exported = list(csv.DictReader(io.StringIO(export_csv(study).decode("utf-8-sig"))))
            if len(exported) != 4:
                raise AssertionError("The CSV omitted a target/drug pair")
            archive = export_package(study)
            with zipfile.ZipFile(archive) as package:
                names = package.namelist()
                if not all(name in names for name in ("results.csv", "references.txt", "methods_and_checks.txt")):
                    raise AssertionError("The study package is missing sources or results")
                if "AutoDock Vina" not in package.read("references.txt").decode("utf-8"):
                    raise AssertionError("The exported references omit the docking tool")
            report.update(passed=True, targets=2, drugs=2, completed_pairs=4,
                          target_specific_boxes=boxes, csv_pairs=len(exported), exported_files=len(names),
                          checks=[{"target": r["target"], "drug": r["drug"], "score_kcal_mol": r["score_kcal_mol"],
                                   "scientific_status": r["scientific_status"]} for r in rows])
        if any(hashlib.sha256((examples / case / "project.db").read_bytes()).hexdigest() != digest
               for case, digest in hashes.items()):
            raise AssertionError("A bundled reference project was modified")
        report["reference_files_unchanged"] = True
    except Exception as exc:
        report.update(passed=False, error=f"{type(exc).__name__}: {exc}", traceback=traceback.format_exc())
    report["seconds"] = round(time.monotonic() - started, 3)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).write_text(json.dumps(report, indent=2), encoding="utf-8")
    return 0 if report["passed"] else 1

"""Persistent, target-specific docking studies using the existing scientific stages.

Every target owns a separate Project, prepared receptor, box and set of poses.
The study matrix retains failed and pending pairs rather than dropping them.
"""
from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import io
import json
import math
from pathlib import Path
import re
import shutil
import tempfile
import threading
from typing import Any, Callable
from uuid import uuid4
import zipfile

from .contract import ControlOutcome, ProvenanceEntry, StageResult, Verdict
from .project import Project

Progress = Callable[[str], None] | None
_locks: dict[str, threading.Lock] = {}
_lock_guard = threading.Lock()
PAIR_COLUMNS = ["target_id", "target", "pdb_id", "drug_id", "drug", "smiles", "prepared_smiles",
                "status", "score_kcal_mol", "scientific_status", "redock_rmsd_angstrom",
                "baseline_check", "redocking_check", "size_check", "check_scope", "error", "pose_file",
                "drug_cid", "drug_source", "drug_source_url", "protein_entry_url", "protein_download_url", "drug_download_url"]


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _slug(name: str) -> str:
    return (re.sub(r"[^A-Za-z0-9_-]+", "_", name).strip("_") or "item")[:55]


@dataclass
class TargetSpec:
    name: str
    pdb_id: str = ""
    file: str = ""
    chains: list[str] | None = None
    site_mode: str = "ligand"
    box: dict | None = None
    reference_project: str = ""


@dataclass
class DrugSpec:
    name: str
    smiles: str
    cid: str | int = ""
    source_url: str = ""
    source: str = ""
    retrieved: str = ""


class Study:
    def __init__(self, folder: Path, data: dict):
        self.folder = Path(folder).resolve()
        self.data = data

    def save(self) -> None:
        self.data["updated"] = _now()
        destination = self.folder / "study.json"
        temporary = self.folder / "study.json.tmp"
        temporary.write_text(json.dumps(self.data, indent=2, ensure_ascii=False), encoding="utf-8")
        temporary.replace(destination)


def create_study(root: Path, name: str, targets: list[TargetSpec | dict],
                 drugs: list[DrugSpec | dict]) -> Study:
    name = name.strip() or "Docking study"
    if not targets or not drugs:
        raise ValueError("Add at least one target and one drug before creating a study.")
    target_specs = [t if isinstance(t, TargetSpec) else TargetSpec(**t) for t in targets]
    drug_specs = [d if isinstance(d, DrugSpec) else DrugSpec(**d) for d in drugs]
    for target in target_specs:
        if not target.name.strip():
            raise ValueError("Give every target a short name.")
        if not (target.pdb_id.strip() or target.file or target.reference_project):
            raise ValueError(f"Add a PDB code or a structure file for {target.name}.")
        if target.pdb_id and not re.fullmatch(r"[A-Za-z0-9_]{4,16}", target.pdb_id.strip()):
            raise ValueError(f"The PDB code for {target.name} contains unsupported characters.")
        if target.file and not Path(target.file).is_file():
            raise ValueError(f"The structure file for {target.name} was not found.")
        _validate_site(target.site_mode, target.box)
    for drug in drug_specs:
        if not drug.name.strip() or not drug.smiles.strip():
            raise ValueError("Every drug needs a name and a SMILES structure.")
    folder = Path(root).resolve() / (f"Study {_slug(name)} {datetime.now().strftime('%Y%m%d-%H%M%S')}-{uuid4().hex[:6]}")
    folder.mkdir(parents=True)
    data = {"version": 1, "name": name, "created": _now(), "updated": _now(),
            "status": "pending", "ui_step": 0, "targets": [], "drugs": [],
            "settings": {"exhaustiveness": 8, "n_decoys": 3},
            "interpretation": "Compare poses and scores within each target's setup. Raw scores across different targets do not establish selectivity or potency."}
    study = Study(folder, data)
    for index, target in enumerate(target_specs, 1):
        entry = asdict(target)
        target_id = f"target-{index:03d}"
        project_relative = f"targets/{index:03d}-{_slug(target.name)}"
        project = Project.create(folder / project_relative, target.name.strip())
        project.close()
        entry.update(id=target_id, name=target.name.strip(), pdb_id=target.pdb_id.strip().upper(),
                     project=project_relative, status="pending", error="", site_summary="", preparation_summary="")
        if target.file:
            source = Path(target.file)
            destination = folder / project_relative / "incoming" / ("structure" + source.suffix.lower())
            shutil.copy2(source, destination)
            entry["file"] = str(destination.relative_to(folder)).replace("\\", "/")
        data["targets"].append(entry)
    for index, drug in enumerate(drug_specs, 1):
        drug_id = f"drug-{index:03d}"
        data["drugs"].append({**asdict(drug), "id": drug_id, "name": drug.name.strip(), "smiles": drug.smiles.strip(),
                              "source": drug.source or "User-supplied SMILES",
                              "internal_name": f"{drug_id}-{_slug(drug.name)}"})
    study.save()
    return study


def open_study(folder: Path | str) -> Study:
    folder = Path(folder).resolve()
    data = json.loads((folder / "study.json").read_text(encoding="utf-8"))
    if data.get("version") != 1 or not isinstance(data.get("targets"), list) or not isinstance(data.get("drugs"), list):
        raise ValueError("This folder does not contain a supported docking study.")
    study = Study(folder, data)
    for target in data["targets"]:
        _inside(study, target["project"])
    if not _lock(study).locked() and (data.get("status") in {"preparing", "docking"}
                                     or any(t.get("status") in {"preparing", "docking"} for t in data["targets"])):
        data["status"] = "interrupted"
        for target in data["targets"]:
            if target.get("status") in {"preparing", "docking"}:
                target.update(status="interrupted", error="The previous calculation was interrupted. Run this step again. Earlier files are retained.")
        study.save()
    return study


def list_studies(root: Path) -> list[dict]:
    studies = []
    for manifest in Path(root).glob("*/study.json"):
        try:
            data = json.loads(manifest.read_text(encoding="utf-8"))
            studies.append({"name": data.get("name", manifest.parent.name), "folder": str(manifest.parent),
                            "updated": data.get("updated", ""), "status": data.get("status", "pending")})
        except (OSError, ValueError):
            continue
    return sorted(studies, key=lambda item: item["updated"], reverse=True)


def _inside(study: Study, relative: str) -> Path:
    path = (study.folder / relative).resolve()
    if not path.is_relative_to(study.folder):
        raise ValueError("A study file points outside its own project folder.")
    return path


def target_project(study: Study, target_id: str) -> Project:
    target = next((row for row in study.data["targets"] if row["id"] == target_id), None)
    if target is None:
        raise ValueError("That target is not part of this study.")
    return Project.open(_inside(study, target["project"]))


def _values(project: Project, stage: str) -> dict:
    row = project.get_step(stage)
    return json.loads(row["values_json"] or "{}") if row and row["status"] == "complete" else {}


def _save_result(project: Project, stage, result: StageResult) -> StageResult:
    project.save_step(stage.stage_id, stage.title, result)
    if not result.success:
        raise ValueError(result.message or result.summary or f"{stage.title} could not finish.")
    return result


def _notify(progress: Progress, message: str) -> None:
    if progress:
        progress(message)


def _lock(study: Study) -> threading.Lock:
    with _lock_guard:
        return _locks.setdefault(str(study.folder), threading.Lock())


def _validate_site(mode: str, box: dict | None) -> None:
    if mode not in {"ligand", "detect", "manual"}:
        raise ValueError("Choose a reference-ligand, detected-pocket or manual search box.")
    if mode == "manual":
        required = ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z")
        if not box or any(key not in box for key in required):
            raise ValueError("A manual box needs three centre coordinates and three side lengths.")
        if any(not math.isfinite(float(box[key])) for key in required):
            raise ValueError("Search-box coordinates and lengths must be finite numbers.")
        if any(not 1 <= float(box[key]) <= 100 for key in required[3:]):
            raise ValueError("Search-box side lengths must be between 1 and 100 angstroms.")


def _import_reference(study: Study, target: dict, project: Project) -> bool:
    """Copy saved preparation honestly without opening a distributed DB writable."""
    source = Path(target.get("reference_project", ""))
    if not target.get("reference_project"):
        return False
    source_hash = (hashlib.sha256((source / "project.db").read_bytes()).hexdigest()
                   if (source / "project.db").is_file() else target.get("reference_hash", ""))
    prepared = _values(project, "03_receptor")
    if (target.get("reference_imported") and target.get("reference_source") == str(source)
            and target.get("reference_hash") == source_hash and prepared
            and Path(prepared.get("prepared_file", "")).is_file()):
        return True
    if not (source / "project.db").is_file():
        raise ValueError("The saved preparation source is unavailable. Add the original structure file instead.")
    from .workshop import copy_steps
    with tempfile.TemporaryDirectory(prefix="study_reference_", dir=study.folder) as temporary:
        temporary = Path(temporary)
        shutil.copy2(source / "project.db", temporary / "project.db")
        shutil.copy2(source / "project.json", temporary / "project.json")
        for stage_id in ("01_target", "02_structure", "03_receptor", "04_site"):
            origin = source / "steps" / stage_id
            if origin.exists():
                shutil.copytree(origin, temporary / "steps" / stage_id, dirs_exist_ok=True)
        snapshot = Project.open(temporary)
        try:
            copy_steps(snapshot, project, ["01_target", "02_structure", "03_receptor", "04_site"])
        finally:
            snapshot.close()
    target["reference_imported"] = True
    target["reference_source"] = str(source)
    target["reference_hash"] = source_hash
    target["source_label"] = "Saved protein preparation; drug preparation and docking are run fresh."
    return True


def _structure(study: Study, target: dict, project: Project) -> None:
    from stages.structure import StructureStage, real_ligands
    from .mmcif import cif_to_pdb
    stage = StructureStage(project)
    project.save_step("01_target", "Selected target", StageResult(
        success=True, summary=f"User-selected structure for {target['name']}",
        values={"target": {"protein": target["name"], "accession": "", "pdb_ids": [target["pdb_id"]] if target["pdb_id"] else []},
                "selection": "Direct structure selection; no target tractability claim was calculated."},
        control=ControlOutcome("Target selection", Verdict.NOT_APPLICABLE,
                               "A structure was selected directly. Target evidence was not assessed.")))
    if target.get("file"):
        source = _inside(study, target["file"])
        text = source.read_text(encoding="utf-8", errors="replace")
        if source.suffix.lower() in {".cif", ".mmcif"}:
            text = cif_to_pdb(text)
        if "ATOM  " not in text:
            raise ValueError("The structure file contains no readable protein ATOM records.")
        path = project.step_folder("02_structure") / "uploaded_structure.pdb"
        path.write_text(text, encoding="utf-8")
        chosen = {"id": target.get("pdb_id") or "uploaded", "source": "uploaded", "local_file": str(path),
                  "title": target["name"], "method": "not supplied", "resolution": None, **real_ligands(text)}
        control = ControlOutcome("Uploaded structure", Verdict.WARN,
                                 "The coordinates were supplied locally. Review the source paper and structure quality before interpreting poses.")
        project.add_provenance(ProvenanceEntry("02_structure", "User-supplied structure", "local file",
                                               parameters={"input_file": source.name, "sha256": hashlib.sha256(source.read_bytes()).hexdigest()},
                                               citation="Original coordinate source and publication must be retained by the study author."))
    else:
        entry = stage.describe_pdb(target["pdb_id"])
        if not entry:
            raise ValueError(f"The Protein Data Bank entry {target['pdb_id']} could not be retrieved. Use a downloaded PDB/mmCIF file to work offline.")
        path, problem = stage._download(entry)
        if path is None:
            raise ValueError(problem)
        chosen = {**entry, "local_file": str(path), **real_ligands(path.read_text(errors="replace"))}
        control = stage._quality_card(chosen)
        project.add_provenance(ProvenanceEntry("02_structure", "RCSB Protein Data Bank", "retrieved structure",
            parameters={"pdb_id": target["pdb_id"], "entry_url": f"https://www.rcsb.org/structure/{target['pdb_id']}",
                        "download_url": f"https://files.rcsb.org/download/{target['pdb_id']}.cif"},
            citation="Berman HM et al. The Protein Data Bank. Nucleic Acids Res. 2000;28:235-242. doi:10.1093/nar/28.1.235"))
    project.save_step("02_structure", "Protein structure", StageResult(success=True,
        summary=f"Loaded {target['name']}", values={"chosen": chosen, "alternatives": []},
        artefacts=[path], control=control))


def _site(study: Study, target: dict, project: Project) -> None:
    from stages.site import SiteStage
    stage = SiteStage(project)
    try:
        _validate_site(target.get("site_mode", "ligand"), target.get("box"))
    except ValueError as exc:
        project.save_step(stage.stage_id, stage.title, StageResult(False, message=str(exc)))
        raise
    receptor = _values(project, "03_receptor")
    if not receptor:
        raise ValueError("Prepare this target before defining its binding site.")
    mode = target.get("site_mode", "ligand")
    ligand = receptor.get("ligand", "") if mode == "ligand" else ""
    if mode == "ligand" and not ligand:
        _save_result(project, stage, StageResult(False, message="No reference ligand was retained for this target. Select a detected pocket or enter its own manual search box."))
    result = _save_result(project, stage, stage.run(pdb_path=receptor["prepared_file"],
        ligand_label=ligand, box=target.get("box") if mode == "manual" else None, padding=4.0))
    target["box"] = result.values.get("box", {})
    target["site_summary"] = result.summary
    target["site_status"] = result.control.verdict.value if result.control else "not_checked"
    target["site_error"] = ""


def run_preparation(study: Study, progress: Progress = None) -> None:
    from stages.receptor import ReceptorStage
    from stages.ligands import LigandStage
    from stages.admet import AdmetStage
    lock = _lock(study)
    if not lock.acquire(blocking=False):
        raise ValueError("This study already has a running calculation. Wait for it to finish.")
    try:
        study.data["status"] = "preparing"
        study.save()
        molecules = [{"name": drug["internal_name"], "smiles": drug["smiles"]} for drug in study.data["drugs"]]
        for target in study.data["targets"]:
            project = None
            try:
                target.update(status="preparing", error="")
                study.save()
                _notify(progress, f"{target['name']}: loading the protein structure")
                project = target_project(study, target["id"])
                input_file = _inside(study, target["file"]) if target.get("file") else None
                reference = Path(target.get("reference_project", "")) / "project.db"
                signature = json.dumps({"file": target.get("file"), "pdb_id": target.get("pdb_id"),
                                        "file_hash": hashlib.sha256(input_file.read_bytes()).hexdigest() if input_file else "",
                                        "chains": target.get("chains"), "reference_project": target.get("reference_project", ""),
                                        "reference_hash": hashlib.sha256(reference.read_bytes()).hexdigest() if target.get("reference_project") and reference.is_file() else target.get("reference_hash", "")}, sort_keys=True)
                current_prepared = _values(project, "03_receptor")
                if not (current_prepared and Path(current_prepared.get("prepared_file", "")).is_file()
                        and target.get("preparation_signature") == signature):
                    imported = _import_reference(study, target, project)
                    if not imported or target.get("chains"):
                        if not imported and not target.get("reference_imported"):
                            _structure(study, target, project)
                        raw = _values(project, "02_structure").get("chosen", {}).get("local_file", "")
                        _notify(progress, f"{target['name']}: adding hydrogens and preparing the receptor")
                        stage = ReceptorStage(project)
                        result = _save_result(project, stage, stage.run(pdb_path=raw, keep_chains=target.get("chains"),
                            remove_waters=True, repair_gaps=False, minimise=False))
                        target["preparation_summary"] = result.summary
                    else:
                        target["preparation_summary"] = project.get_step("03_receptor")["summary"]
                    target["preparation_signature"] = signature
                _notify(progress, f"{target['name']}: defining this target's own search box")
                try:
                    _site(study, target, project)
                except Exception as exc:
                    target.update(site_error=str(exc), site_summary="Binding site needs setup: " + str(exc), site_status="needs_input")
                _notify(progress, f"{target['name']}: preparing {len(molecules)} drug structures")
                project.add_provenance(ProvenanceEntry("05_ligands", "Drug input records", "study inputs",
                    parameters={"drugs": [{key: drug.get(key, "") for key in ("id", "name", "smiles", "cid", "source", "source_url", "retrieved")} for drug in study.data["drugs"]]},
                    citation="PubChem input requests are identified by their recorded CID and URL; manually supplied structures are labelled as user-supplied."))
                stage = LigandStage(project)
                _save_result(project, stage, stage.run(molecules=molecules))
                stage = AdmetStage(project)
                _save_result(project, stage, stage.run(use_online=False))
                target["warnings"] = []
                for stage_id in ("02_structure", "03_receptor", "04_site", "05_ligands", "08_admet"):
                    recorded = project.get_step(stage_id)
                    control = json.loads(recorded["control_json"] or "{}") if recorded else {}
                    if control.get("verdict") in {"warn", "fail"}:
                        target["warnings"].append({"stage": stage_id, "status": control["verdict"],
                                                   "message": control.get("plain_verdict", ""),
                                                   "what_to_do": control.get("what_to_do", "")})
                if target.get("site_error"):
                    target["warnings"].append({"stage": "04_site", "status": "needs_input",
                                               "message": target["site_error"],
                                               "what_to_do": "Define this target's binding site in the next step."})
                target["preparation_status"] = ("fail" if any(w["status"] == "fail" for w in target["warnings"])
                                                  else "warn" if target["warnings"] else "pass")
                target["status"] = "prepared"
                target["error"] = ""
            except Exception as exc:
                target.update(status="failed", error=str(exc))
            finally:
                if project:
                    project.close()
                study.save()
        states = [target["status"] for target in study.data["targets"]]
        study.data["status"] = "prepared" if all(s == "prepared" for s in states) else "partial" if "prepared" in states else "failed"
        study.save()
    finally:
        lock.release()


def define_binding_sites(study: Study, progress: Progress = None) -> None:
    lock = _lock(study)
    if not lock.acquire(blocking=False):
        raise ValueError("This study already has a running calculation. Wait for it to finish.")
    try:
        for target in study.data["targets"]:
            project = None
            try:
                _notify(progress, f"{target['name']}: updating its search box")
                project = target_project(study, target["id"])
                _site(study, target, project)
                if _values(project, "05_ligands"):
                    target.update(status="prepared", error="")
                else:
                    target.update(status="pending", error="Prepare the drug structures before docking this target.")
            except Exception as exc:
                target.update(status="failed", error=str(exc), site_error=str(exc),
                              site_status="needs_input", site_summary="Binding site needs setup: " + str(exc))
            finally:
                if project:
                    project.close()
                study.save()
        study.data["status"] = "prepared" if all(t["status"] == "prepared" for t in study.data["targets"]) else "partial"
        study.save()
    finally:
        lock.release()


def run_docking(study: Study, exhaustiveness: int = 8, n_decoys: int = 3,
                progress: Progress = None) -> None:
    from stages.dock import DockStage
    exhaustiveness, n_decoys = int(exhaustiveness), int(n_decoys)
    if not 1 <= exhaustiveness <= 128 or not 1 <= n_decoys <= 20:
        raise ValueError("Use search effort 1–128 and 1–20 comparison molecules.")
    lock = _lock(study)
    if not lock.acquire(blocking=False):
        raise ValueError("This study already has a running calculation. Wait for it to finish.")
    try:
        study.data["status"] = "docking"
        study.data["settings"] = {"exhaustiveness": exhaustiveness, "n_decoys": n_decoys}
        study.save()
        for target in study.data["targets"]:
            project = None
            try:
                _notify(progress, f"{target['name']}: docking the prepared drugs and running scientific checks")
                project = target_project(study, target["id"])
                for required in ("02_structure", "03_receptor", "04_site", "05_ligands"):
                    row = project.get_step(required)
                    if not row or row["status"] != "complete":
                        raise ValueError("This target needs current protein, binding-site and drug preparation before docking.")
                    control = json.loads(row["control_json"] or "{}")
                    if control.get("verdict") == "fail":
                        raise ValueError(control.get("plain_verdict") or "An input scientific check blocks docking for this target.")
                target.update(status="docking", error="")
                study.save()
                stage = DockStage(project)
                result = _save_result(project, stage, stage.run(exhaustiveness=exhaustiveness, n_decoys=n_decoys))
                target["docking_summary"] = result.summary
                target["status"] = "partial" if any("error" in row for row in result.values.get("results", [])) else "docked"
            except Exception as exc:
                target.update(status="failed", error=str(exc))
            finally:
                if project:
                    project.close()
                study.save()
        rows = results_rows(study)
        complete = sum(row["status"] == "completed" for row in rows)
        study.data["status"] = "completed" if complete == len(rows) else "partial" if complete else "failed"
        study.save()
    finally:
        lock.release()


def results_rows(study: Study) -> list[dict]:
    rows = []
    for target in study.data["targets"]:
        project = None
        try:
            project = target_project(study, target["id"])
            dock_row = project.get_step("06_dock")
            dock = _values(project, "06_dock")
            preparation = _values(project, "05_ligands")
            control = json.loads(dock_row["control_json"] or "{}") if dock_row and dock_row["status"] == "complete" else {}
            by_name = {entry["name"]: entry for entry in dock.get("results", [])}
            prepared_by_name = {entry["name"]: entry for entry in preparation.get("prepared", [])}
            failed = {entry["name"]: entry.get("error", "Drug preparation failed.") for entry in preparation.get("failed", [])}
            for drug in study.data["drugs"]:
                result = by_name.get(drug["internal_name"], {})
                stale = bool(dock_row and dock_row["status"] == "stale")
                failed_attempt = target.get("status") == "failed" and bool(target.get("error"))
                active = target.get("status") in {"preparing", "docking"}
                interrupted = target.get("status") == "interrupted"
                not_current = stale or failed_attempt or active or interrupted
                pair_error = result.get("error") or failed.get(drug["internal_name"], "") or target.get("error", "")
                status = ("needs_updating" if stale else "running" if active else "interrupted" if interrupted else "failed" if failed_attempt else "completed" if "score" in result else "failed" if pair_error else "pending")
                checks = dock.get("checks", {}) if not not_current else {}
                check_default = "not_current" if not_current and dock_row else "not_checked"
                pose = result.get("file", "")
                if not_current:
                    pose = ""
                if pose:
                    try:
                        pose = str(Path(pose).resolve().relative_to(study.folder)).replace("\\", "/")
                    except ValueError:
                        pose = ""
                rows.append({"target_id": target["id"], "target": target["name"], "pdb_id": target.get("pdb_id", ""),
                    "drug_id": drug["id"], "drug": drug["name"], "smiles": drug["smiles"], "status": status,
                    "prepared_smiles": prepared_by_name.get(drug["internal_name"], {}).get("smiles", ""),
                    "score_kcal_mol": result.get("score", "") if not not_current else "",
                    "scientific_status": control.get("verdict", "not_checked") if not not_current else check_default,
                    "redock_rmsd_angstrom": dock.get("redock_rmsd") if not not_current and dock.get("redock_rmsd") is not None else "",
                    "baseline_check": checks.get("baseline", {}).get("verdict", check_default),
                    "redocking_check": checks.get("re-docking", {}).get("verdict", check_default),
                    "size_check": checks.get("size", {}).get("verdict", check_default),
                    "check_scope": "Target run; baseline concerns its best-scoring drug, not every pair.",
                    "error": pair_error, "pose_file": pose, **_source_links(target, drug)})
        except Exception as exc:
            for drug in study.data["drugs"]:
                rows.append({"target_id": target["id"], "target": target["name"], "pdb_id": target.get("pdb_id", ""),
                             "drug_id": drug["id"], "drug": drug["name"], "smiles": drug["smiles"],
                             "status": "failed", "score_kcal_mol": "", "scientific_status": "not_checked", "error": str(exc), **_source_links(target, drug)})
        finally:
            if project:
                project.close()
    return rows


def _source_links(target: dict, drug: dict) -> dict:
    from .references import drug_links, structure_links
    protein_links = structure_links(target.get("pdb_id", ""))
    molecule_links = drug_links(drug["name"], drug["smiles"], drug.get("cid", ""))
    return {"drug_cid": drug.get("cid", ""), "drug_source": drug.get("source", "User-supplied SMILES"),
            "drug_source_url": drug.get("source_url", ""),
            "protein_entry_url": protein_links[0]["url"] if protein_links else "",
            "protein_download_url": protein_links[1]["url"] if len(protein_links) > 1 else "",
            "drug_download_url": molecule_links[1]["url"] if len(molecule_links) > 1 else ""}


def export_csv(study: Study) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=PAIR_COLUMNS, extrasaction="ignore")
    writer.writeheader()
    for row in results_rows(study):
        safe = dict(row)
        for key in ("target", "drug", "smiles", "error"):
            if str(safe.get(key, "")).lstrip().startswith(("=", "+", "-", "@")):
                safe[key] = "'" + str(safe[key])
        writer.writerow(safe)
    return stream.getvalue().encode("utf-8-sig")


def methods_text(study: Study) -> str:
    lines = ["Docking study methods and scientific checks", study.data["name"],
             "Generated: " + _now(), "", study.data["interpretation"], "",
             "Actual run settings: " + json.dumps(study.data.get("settings", {})), ""]
    lines += ["Original drug inputs and their sources:", json.dumps(study.data["drugs"], ensure_ascii=False), ""]
    for target in study.data["targets"]:
        lines += [f"Target {target['id']}: {target['name']}", "PDB: " + target.get("pdb_id", ""),
                  "Status: " + target.get("status", "pending"), "Error: " + target.get("error", ""),
                  "Search box: " + json.dumps(target.get("box")), ""]
        project = target_project(study, target["id"])
        try:
            for entry in project.provenance():
                lines.append(json.dumps(entry, ensure_ascii=False, default=str))
            for stage_id in ("02_structure", "03_receptor", "04_site", "05_ligands", "06_dock", "08_admet"):
                row = project.get_step(stage_id)
                if row:
                    lines += ["", f"{stage_id}: {row['status']} — {row['summary']}", row["control_json"] or "No scientific check recorded."]
                    if row["status"] == "complete" and stage_id == "06_dock":
                        lines.append("Individual docking checks: " + json.dumps(json.loads(row["values_json"] or "{}").get("checks", {}), ensure_ascii=False))
                    if row["status"] == "complete" and stage_id in {"05_ligands", "08_admet"}:
                        lines.append("Current calculated structures/descriptors: " + (row["values_json"] or "{}"))
        finally:
            project.close()
        lines.append("")
    return "\n".join(lines) + "\n"


def export_package(study: Study) -> Path:
    from .references import render_references
    study.save()
    projects = [target_project(study, target["id"]) for target in study.data["targets"]]
    try:
        references = render_references(projects)
    finally:
        for project in projects:
            project.close()
    destination = study.folder / "exports" / "docking_study.zip"
    destination.parent.mkdir(exist_ok=True)
    with zipfile.ZipFile(destination, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
        bundle.writestr("results.csv", export_csv(study))
        bundle.writestr("methods_and_checks.txt", methods_text(study))
        bundle.writestr("references.txt", references)
        bundle.write(study.folder / "study.json", "study.json")
        for path in sorted((study.folder / "targets").rglob("*")):
            if path.is_file():
                bundle.write(path, str(path.relative_to(study.folder)).replace("\\", "/"))
    return destination


def create_example_study(root: Path, name: str = "EGFR example") -> Study:
    from .workshop import EXAMPLES, MOLECULES
    return create_study(root, name, [TargetSpec("EGFR", pdb_id="1M17", reference_project=str(EXAMPLES / "egfr"))],
                        [DrugSpec(molecule["name"], molecule["smiles"]) for molecule in MOLECULES])

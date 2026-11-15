"""
Re-run protein preparation on the committed worked examples.

    python regenerate_examples.py            # report what would change
    python regenerate_examples.py --write    # do it

Why this exists. The examples under assets/workshop are pre-baked
project folders: a student pressing "Use the offline EGFR example"
gets a copy of them, so whatever is committed there is what they see.
They were generated before preparation was fixed to actually remove
water, so their prepared structures still contain it while their saved
step reports it as removed. A code fix cannot reach them. They have to
be made again.

This drives the real ReceptorStage.run() against each template project
and saves through project.save_step(), so values, summary, control and
provenance are all written by the same code path the application uses.
Nothing here edits the database by hand.

It re-runs step 3 only. Everything downstream -- the search box, the
docked poses, the ADMET results -- is left alone, and is unaffected,
because the receptor the docking actually used is built separately by
core/receptor_pdbqt.py, which has always dropped water.

Afterwards, run:   python verify_example_assets.py
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
WORKSHOP = HERE / "assets" / "workshop"
WATER = {"HOH", "WAT", "DOD", "H2O"}


def water_atoms(path: Path) -> int:
    return sum(1 for line in path.read_text(errors="replace").splitlines()
               if line.startswith(("ATOM", "HETATM"))
               and line[17:20].strip().upper() in WATER)


def say(text: str = "") -> None:
    print(text, flush=True)


def portable_location(case: Path) -> str:
    """
    The placeholder path the committed example stores, e.g.
    "/bundled-examples/egfr".

    The examples do NOT record the machine they were built on. They
    record a neutral prefix, and Project.relocate_paths() rewrites it
    to wherever the student's copy actually sits. Opening a project
    rebases it to the real folder, so regenerating here would bake a
    developer's own path into the shipped asset and leave every
    student's copy pointing at a directory that does not exist. This
    captures the placeholder before anything opens the project, so it
    can be put back afterwards.
    """
    import json
    try:
        return json.loads((case / "project.json").read_text(
            encoding="utf-8")).get("location", "")
    except (OSError, ValueError):
        return ""


def restore_location(case: Path, placeholder: str) -> None:
    """Put the neutral prefix back into project.json and the database."""
    import json
    import sqlite3
    if not placeholder:
        return
    real = str(case.resolve())

    summary = json.loads((case / "project.json").read_text(encoding="utf-8"))
    summary["location"] = placeholder
    (case / "project.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8")

    connection = sqlite3.connect(case / "project.db")
    try:
        for table, columns in (("steps", ("values_json", "control_json")),
                               ("provenance", ("parameters",))):
            existing = {row[1] for row in
                        connection.execute(f"PRAGMA table_info({table})")}
            for column in columns:
                if column in existing:
                    connection.execute(
                        f"UPDATE {table} SET {column} = "
                        f"REPLACE({column}, ?, ?)", (real, placeholder))
        connection.commit()
    finally:
        connection.close()


def _rerun_downstream(project, receptor_values: dict, say) -> None:
    """Re-run the binding site and the docking so no step is left stale."""
    import json
    import shutil
    from stages.site import SiteStage
    from stages.dock import DockStage, find_vina

    receptor = json.loads(
        project.get_step("03_receptor")["values_json"] or "{}")

    site = SiteStage(project)
    old_box = json.loads(project.get_step("04_site")["values_json"] or "{}"
                         ).get("box", {})
    result = site.run(pdb_path=receptor["prepared_file"],
                      ligand_label=receptor.get("ligand", ""), padding=4.0)
    if not result.success:
        say(f"      site step FAILED: {result.message}")
        return
    project.save_step(site.stage_id, site.title, result)
    for entry in site.provenance:
        project.add_provenance(entry)
    box = result.values.get("box", {})
    moved = max((abs(box.get(k, 0) - old_box.get(k, 0))
                 for k in ("center_x", "center_y", "center_z")), default=0)
    say(f"      search box re-derived, centre moved {moved:.2f} A")

    if not find_vina():
        say("      AutoDock Vina not found, so docking was NOT re-run.")
        say("      Step 6 is left stale. Install Vina and run this again,")
        say("      or the worked example will show 'needs updating'.")
        return

    dock = DockStage(project)
    before = {r["name"]: r.get("score")
              for r in json.loads(project.get_step("06_dock")["values_json"]
                                  or "{}").get("results", [])}
    result = dock.run()
    if not result.success:
        say(f"      docking FAILED: {result.message}")
        return
    project.save_step(dock.stage_id, dock.title, result)
    for entry in dock.provenance:
        project.add_provenance(entry)
    say("      docking re-run:")
    for row in result.values.get("results", []):
        was = before.get(row["name"])
        now = row.get("score")
        if isinstance(now, (int, float)):
            shift = f" (was {was})" if isinstance(was, (int, float)) else " (new)"
            say(f"         {row['name']:<28} {now:>8.3f}{shift}")
    rmsd = result.values.get("redock_rmsd")
    if rmsd is not None:
        say(f"      re-docking RMSD: {rmsd} A")


def regenerate(case: Path, write: bool) -> bool:
    """True if this example needed changing."""
    from core.project import Project
    from stages.receptor import ReceptorStage

    placeholder = portable_location(case)
    project = Project.open(case)
    try:
        step = project.get_step("03_receptor")
        if not step or step["status"] != "complete":
            say(f"  {case.name}: no completed preparation step, skipped")
            return False

        import json
        values = json.loads(step["values_json"] or "{}")
        prepared = Path(values.get("prepared_file", ""))
        before = water_atoms(prepared) if prepared.is_file() else 0
        if not before:
            say(f"  {case.name}: already dry, nothing to do")
            return False

        structure = json.loads(
            project.get_step("02_structure")["values_json"] or "{}")
        source = structure.get("chosen", {}).get("local_file", "")
        if not Path(source).is_file():
            say(f"  {case.name}: original structure {source} is missing, "
                "skipped")
            return False

        # Exactly the settings the saved step records, so that nothing
        # but the water changes.
        chains = list(values.get("chains", {})) or None
        kept = [h["label"] for h in values.get("kept", [])]
        gaps = bool(values.get("gaps_repaired", 0)) or bool(
            values.get("gaps_left_open"))

        say(f"  {case.name}: {before} water atoms in "
            f"{prepared.name}; re-running preparation")
        say(f"      chains={chains} keep={kept} repair_gaps={gaps}")
        if not write:
            say("      (dry run -- pass --write to apply)")
            return True

        stage = ReceptorStage(project)
        result = stage.run(pdb_path=source, keep_chains=chains,
                           keep_hets=kept, remove_waters=True,
                           repair_gaps=gaps, max_gap_length=3,
                           minimise=False)
        if not result.success:
            say(f"      FAILED: {result.message or result.summary}")
            return True
        project.save_step(stage.stage_id, stage.title, result)
        for entry in stage.provenance:
            project.add_provenance(entry)

        after_path = Path(result.values.get("prepared_file", ""))
        say(f"      {result.summary}")
        say(f"      water atoms now: {water_atoms(after_path)}")

        # Preparation feeds the search box, and the box feeds the
        # docking, so saving a new step 3 correctly marks steps 4 and 6
        # stale. An example with a stale step shows the student "needs
        # updating" instead of a result, so the chain has to be
        # finished, not just started.
        _rerun_downstream(project, values, say)
        return True
    finally:
        project.close()
        if write and placeholder:
            restore_location(case, placeholder)
            say(f"      paths reset to {placeholder}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true",
                        help="apply the changes instead of reporting them")
    args = parser.parse_args()

    sys.path.insert(0, str(HERE))
    cases = sorted(p for p in WORKSHOP.glob("*") if (p / "project.db").is_file())
    if not cases:
        say(f"  No example projects under {WORKSHOP}")
        return 1

    say()
    say("  Regenerating worked-example preparation")
    say("  " + "=" * 56)
    changed = [case for case in cases if regenerate(case, args.write)]
    say()
    if not changed:
        say("  Every example is already consistent with its own report.")
    elif args.write:
        say(f"  {len(changed)} example(s) rewritten. Now run:")
        say("      python verify_example_assets.py")
        say("  and commit the new prepared .pdb files with their project.db.")
    else:
        say(f"  {len(changed)} example(s) would change. Re-run with --write.")
    say()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""
Check the committed worked examples against what the app claims about them.

    python verify_example_assets.py

The prepared structures under assets/workshop are not test fixtures.
They are what a student opens when they press "Use the offline EGFR
example", and what the "Protein after preparation" panel draws. If the
step says the waters were removed and the file still has them, the
example teaches the opposite of the lesson.

This is a separate script rather than a unit test because a failure
here is fixed by regenerating data, not by changing code, and a test
that stays red until somebody finds time to re-run an example trains
people to ignore a red suite.

To regenerate after fixing preparation: open each example in the app,
re-run step 3 on the structure the project already has, and commit the
new assets/workshop/<case>/steps/03_receptor/*.pdb together with the
updated project.db.
"""

from __future__ import annotations

from pathlib import Path
import sqlite3
import sys

HERE = Path(__file__).resolve().parent
WORKSHOP = HERE / "assets" / "workshop"
WATER = {"HOH", "WAT", "DOD", "H2O"}


def atoms(path: Path) -> list[str]:
    return [line[17:20].strip().upper()
            for line in path.read_text(errors="replace").splitlines()
            if line.startswith(("ATOM", "HETATM"))]


def say(text: str = "") -> None:
    print(text, flush=True)


def reported_waters_removed(case: Path) -> int | None:
    """What the saved step claims, so claim and file can be compared."""
    db = case / "project.db"
    if not db.is_file():
        return None
    try:
        con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
        row = con.execute(
            "select values_json from steps where stage_id='03_receptor'"
        ).fetchone()
        con.close()
    except sqlite3.Error:
        return None
    if not row or not row[0]:
        return None
    import json
    try:
        return json.loads(row[0]).get("waters_removed")
    except ValueError:
        return None


def main() -> int:
    if not WORKSHOP.is_dir():
        say(f"  No workshop examples at {WORKSHOP}")
        return 1

    prepared = sorted(WORKSHOP.glob("*/steps/03_receptor/*_prepared.pdb"))
    if not prepared:
        say("  No prepared example structures found.")
        return 1

    say()
    say("  Checking committed worked examples")
    say("  " + "=" * 56)

    problems = []
    for path in prepared:
        case = path.parents[2]
        names = atoms(path)
        wet = sum(1 for n in names if n in WATER)
        claimed = reported_waters_removed(case)
        claim = ("the step reports no water setting" if claimed is None
                 else f"the step reports {claimed} water(s) removed")
        if wet:
            problems.append((path, wet, claim))
            say(f"  FAIL  {case.name}/{path.name}")
            say(f"        {wet} water atoms are still in the file, but "
                f"{claim}.")
        else:
            say(f"  ok    {case.name}/{path.name}  ({len(names)} atoms, "
                "no water)")

    say()
    if problems:
        say(f"  {len(problems)} example structure(s) disagree with their own "
            "report.")
        say()
        say("  Regenerate them: open each example in the app and re-run")
        say("  step 3 (Prepare protein) on the structure already chosen,")
        say("  then commit the new prepared .pdb and project.db together.")
        say()
        return 1

    say("  Every committed example matches what the app says about it.")
    say()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

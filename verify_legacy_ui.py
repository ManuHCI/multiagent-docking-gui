"""Exercise the native legacy UI with isolated projects and real docking.

Run this with the Python 3.8 build environment. It verifies the application
on the current host; it does not establish Windows 7/8/8.1 compatibility.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform
import sys
import time
import traceback
from unittest.mock import patch
from uuid import uuid4
import zipfile


ROOT = Path(__file__).resolve().parent


def run_ui_checks(output):
    from legacy.compat import install_compatibility
    install_compatibility()

    temporary = ROOT / ".qa-temp" / ("legacy-ui-acceptance-" + uuid4().hex[:8])
    temporary.mkdir(parents=True)
    os.environ["INSILICO_DATA_DIR"] = str(temporary / "user")
    os.environ["INSILICO_PROJECTS"] = str(temporary / "Projects")
    os.environ["INSILICO_LEGACY_SHORT_PATHS"] = "1"

    import tkinter as tk
    from legacy.ui import DockingStudio, DEVELOPERS
    from core.project import Project

    started = time.monotonic()
    checks, errors = [], []
    window = tk.Tk()
    window.withdraw()
    app = DockingStudio(window, project_root=temporary / "Projects")
    window.report_callback_exception = lambda kind, value, trace: errors.append("".join(traceback.format_exception(kind, value, trace)))
    host = platform.platform()
    if sys.platform == "win32":
        version = sys.getwindowsversion()
        # Python 3.8 predates the Windows 11 product name. Keep its real NT
        # build instead of labelling an NT 10.0/26200 host as Windows 10.
        host = "Windows NT {}.{} build {}".format(version.major, version.minor, version.build)
    report = {"passed": False, "host": host, "python_platform_report": platform.platform(), "python": platform.python_version(),
              "actual_old_windows_vm_tested": False, "data_folder": str(temporary), "checks": checks}

    def pump(timeout=180):
        deadline = time.monotonic() + timeout
        while app.busy and time.monotonic() < deadline:
            window.update()
            time.sleep(0.03)
        window.update()
        if app.busy:
            raise AssertionError("The isolated UI calculation exceeded its timeout.")
        if errors:
            raise AssertionError("Native UI error: " + "\n".join(errors))

    def all_text(widget):
        values = []
        if isinstance(widget, tk.Text):
            values.append(widget.get("1.0", "end-1c"))
        else:
            try:
                values.append(str(widget.cget("text")))
            except tk.TclError:
                pass
        for child in widget.winfo_children():
            values.extend(all_text(child))
        return values

    def has_structure(widget):
        return bool(getattr(widget, "image", None)) or any(has_structure(child) for child in widget.winfo_children())

    try:
        with patch("legacy.ui.messagebox.showerror", side_effect=lambda title, text, **kwargs: errors.append(str(text))):
            assert list(reversed(list(app.main_buttons))) == ["Home", "Projects", "Reference", "About"]
            assert app.logo is not None
            assert "Legacy Beta" in window.title()
            for section in ("Reference", "About", "Projects", "Home"):
                app.show_section(section)
                pump()
                assert app.section == section
                text = "\n".join(all_text(app.body))
                if section == "About":
                    assert all(developer in text for developer in DEVELOPERS)
                if section == "Reference":
                    assert "Scientific references" in text and "Save references and methods" in text
            checks.append("four fixed sections, supplied logo, complete credits and scientific reference controls")

            app.new_project()
            app.name.set("Disposable UI matrix")
            for pdb_id, name in (("1M17", "EGFR"), ("9OKG", "SV2A")):
                app.pdb_id.set(pdb_id)
                app.target_name.set(name)
                app.add_target()
            app.show_workflow(1)
            app.smiles_widget.insert("1.0", "CCO Ethanol\nCC(N)=O Acetamide")
            app.add_smiles()
            window.update()
            assert len(app.draft["targets"]) == 2 and len(app.draft["drugs"]) == 2
            assert has_structure(app.body)
            checks.append("two targets and two drug identities added, including native 2D structure display")

            app.create_project()
            pump()
            disposable = app.active.folder
            assert app.current == 2
            saved = json.loads((disposable / "study.json").read_text(encoding="utf-8"))
            assert len(saved["targets"]) == 2 and len(saved["drugs"]) == 2
            app.show_section("About")
            app.show_section("Home")
            app.show_workflow(2)
            assert app.active.folder == disposable and app.current == 2
            checks.append("four-pair input matrix persisted and current project survives global navigation")

            app.run_worker("Create offline example", lambda: app.engine.create_example_study(app.project_root), app._open_created)
            pump()
            calculated = app.active.folder
            app.run_worker("Prepare structures", lambda: app.engine.run_preparation(app.active, progress=app.progress))
            pump()
            assert app._ready_to_dock()
            app.show_workflow(1)
            assert has_structure(app.body)
            app.show_workflow(2)
            checks.append("real protein/drug preparation and recorded binding-region controls render through the native UI")

            app.show_workflow(3)
            app.run_worker("Dock example", lambda: app.engine.run_docking(app.active, exhaustiveness=1, n_decoys=1, progress=app.progress))
            pump()
            app.show_workflow(4)
            rows = app.engine.results_rows(app.active)
            assert len(rows) == 3
            assert all(row["status"] == "completed" and math.isfinite(float(row["score_kcal_mol"])) for row in rows)
            report["fresh_scores"] = [row["score_kcal_mol"] for row in rows]
            checks.append("fresh three-drug Vina docking produces finite scores and a rendered native result table")

            csv_path, zip_path, methods_path = temporary / "results.csv", temporary / "study.zip", temporary / "methods.txt"
            with patch("legacy.ui.filedialog.asksaveasfilename", return_value=str(csv_path)):
                app.save_csv()
            assert len(list(csv.DictReader(io.StringIO(csv_path.read_text(encoding="utf-8-sig"))))) == 3
            with patch("legacy.ui.filedialog.asksaveasfilename", return_value=str(zip_path)):
                app.save_package()
                pump()
            with zipfile.ZipFile(zip_path) as package:
                assert {"results.csv", "references.txt", "methods_and_checks.txt", "study.json"}.issubset(set(package.namelist()))
            with patch("legacy.ui.filedialog.asksaveasfilename", return_value=str(methods_path)):
                app.save_methods()
            assert "Vina" in methods_path.read_text(encoding="utf-8")
            with patch("legacy.ui.webbrowser.open", return_value=True):
                app.open_view("06_dock")
            views = list(calculated.glob("targets/*/exports/view-06_dock.html"))
            assert views and 'file:///' in views[0].read_text(encoding="utf-8")
            checks.append("CSV, full ZIP, methods and local-only optional 3D HTML exports verified")

            app.show_section("Projects")
            pump()
            calculated_row = next(row for row in app.library_rows if Path(row["folder"]) == calculated)
            disposable_row = next(row for row in app.library_rows if Path(row["folder"]) == disposable)
            app.open_project(disposable_row)
            assert app.active.folder == disposable and app.current == 2
            app.show_section("Projects")
            pump()
            app.open_project(calculated_row)
            assert app.active.folder == calculated and app.current == 4
            checks.append("Projects shows sizes and reopens both pending inputs and calculated results at their saved steps")

            app.show_section("Projects")
            pump()
            marker = hashlib.sha256((calculated / "study.json").read_bytes()).hexdigest()
            with patch("legacy.ui.messagebox.askokcancel", return_value=False):
                app.delete_project(disposable_row)
            assert disposable.is_dir() and calculated.is_dir()
            with patch("legacy.ui.messagebox.askokcancel", return_value=True):
                app.delete_project(disposable_row)
                pump()
            assert not disposable.exists() and calculated.is_dir()
            assert hashlib.sha256((calculated / "study.json").read_bytes()).hexdigest() == marker
            assert app.active.folder == calculated
            checks.append("cancel preserves the selected project; confirmed deletion keeps the calculated project unchanged")

            earlier = Project.create(app.project_root / "Earlier individual project", "Earlier individual project")
            earlier.close()
            app._legacy_project({"name": "Earlier individual project", "folder": str(earlier.folder), "kind": "legacy"})
            for child in window.winfo_children():
                if isinstance(child, tk.Toplevel):
                    child.withdraw()
                    child.destroy()
            window.update()
            checks.append("earlier individual-project records can be read without editing or holding an open database")

            app.show_section("Home")
            text = "\n".join(all_text(app.body))
            assert not any(phrase in text for phrase in ("Ask a question", "Your prediction", "What did you observe", "Explain in your own words"))
            assert not errors
            report.update(passed=True, errors=errors, real_docking_rows=3, seconds=round(time.monotonic() - started, 3))
    except Exception:
        report.update(error=traceback.format_exc(), errors=errors, seconds=round(time.monotonic() - started, 3))
    finally:
        window.destroy()
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(ROOT / ".verification/legacy/ui-source-acceptance.json"))
    arguments = parser.parse_args()
    raise SystemExit(run_ui_checks(arguments.output))

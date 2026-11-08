"""Isolated native legacy desktop entry point for Python 3.8 Windows builds."""
from __future__ import annotations

import argparse
import importlib
import json
import multiprocessing
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback

VERSION = "1.3.0-legacy.1"
TITLE = "In-Silico Docking Studio — Legacy Beta"


def bootstrap():
    default = Path(os.environ.get("LOCALAPPDATA", str(Path.home()))) / "InSilicoDockingStudioLegacy"
    os.environ.setdefault("INSILICO_DATA_DIR", str(default))
    os.environ["INSILICO_LEGACY_SHORT_PATHS"] = "1"
    from legacy.compat import install_compatibility
    install_compatibility()
    from core.runtime import configure, data_root
    configure()
    logs = data_root() / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    if getattr(sys, "frozen", False):
        stream = (logs / "legacy-runtime.log").open("a", encoding="utf-8", buffering=1)
        sys.stdout = sys.stderr = stream


def diagnostic(output=""):
    findings = []
    for name in ("rdkit", "numpy", "scipy", "openmm", "openmm.app", "pdbfixer", "meeko", "PIL", "tkinter", "certifi"):
        try:
            module = importlib.import_module(name)
            findings.append({"name": name, "ready": True, "version": str(getattr(module, "__version__", "bundled"))})
        except Exception as error:
            findings.append({"name": name, "ready": False, "error": str(error)})
    from stages.dock import find_vina
    vina = find_vina()
    try:
        result = subprocess.run([vina, "--version"], capture_output=True, text=True, timeout=20,
                                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        findings.append({"name": "AutoDock Vina", "ready": result.returncode == 0,
                         "version": (result.stdout or result.stderr).strip(), "path": vina})
    except Exception as error:
        findings.append({"name": "AutoDock Vina", "ready": False, "error": str(error)})
    report = {"version": VERSION, "legacy_beta": True,
              "can_start": all(item["ready"] for item in findings),
              "python": sys.version, "host": platform.platform(), "machine": platform.machine(),
              "target": "Windows 7 SP1, Windows 8 and Windows 8.1; 64-bit",
              "target_os_validated": False,
              "validation_note": "Host checks do not certify execution on the original Windows 7/8/8.1 systems.",
              "interface": "Native Tk; no WebView2 dependency", "findings": findings}
    if output:
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        Path(output).write_text(json.dumps(report, indent=2), encoding="utf-8")
    if sys.stdout is not None:
        print(json.dumps(report, indent=2))
    return 0 if report["can_start"] else 1


def window_check(output):
    import tempfile
    import tkinter as tk
    from legacy.ui import DockingStudio, DEVELOPERS
    from core.runtime import data_root
    destination = Path(output).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    report = {"passed": False, "native_tk": True, "legacy_beta": True, "target_os_validated": False}
    with tempfile.TemporaryDirectory(prefix="legacy-window-", dir=str(data_root())) as temporary:
        root = tk.Tk()
        app = DockingStudio(root, project_root=Path(temporary) / "Projects")

        def text_values(widget):
            values = []
            try:
                text = widget.cget("text")
                if text:
                    values.append(str(text))
            except tk.TclError:
                pass
            for child in widget.winfo_children():
                values.extend(text_values(child))
            return values

        def check():
            try:
                root.update_idletasks()
                nav = [name for name, button in sorted(app.main_buttons.items(), key=lambda item: item[1].winfo_x())]
                assert nav == ["Home", "Projects", "Reference", "About"], nav
                assert app.logo is not None
                home = text_values(app.body)
                assert "Start new project" in home
                app.show_section("About")
                root.update_idletasks()
                about = text_values(app.body)
                assert all(developer in about for developer in DEVELOPERS)
                app.show_section("Reference")
                root.update_idletasks()
                assert any("PubChem" in text for text in text_values(app.body))
                app.show_section("Home")
                root.update_idletasks()
                screenshot = destination.with_suffix(".png")
                try:
                    from PIL import ImageGrab
                    bounds = (root.winfo_rootx(), root.winfo_rooty(), root.winfo_rootx() + root.winfo_width(), root.winfo_rooty() + root.winfo_height())
                    ImageGrab.grab(bbox=bounds).save(str(screenshot))
                    report["screenshot"] = str(screenshot)
                except Exception as error:
                    report["screenshot_error"] = str(error)
                report.update(passed=True, navigation=nav, logo=True, credits=True, references=True,
                              no_webview_required=True)
            except Exception as error:
                report.update(error=str(error), traceback=traceback.format_exc())
            finally:
                destination.write_text(json.dumps(report, indent=2), encoding="utf-8")
                root.destroy()

        root.after(450, check)
        root.mainloop()
    return 0 if report["passed"] else 1


def main():
    bootstrap()
    parser = argparse.ArgumentParser(description=TITLE)
    parser.add_argument("--diagnose", action="store_true")
    parser.add_argument("--verify-output", default="")
    parser.add_argument("--verify-science", default="")
    parser.add_argument("--verify-study", default="")
    parser.add_argument("--window-test-output", default="", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.verify_output or args.diagnose:
        return diagnostic(args.verify_output)
    if args.verify_science:
        from core.desktop_checks import run_scientific_checks
        return run_scientific_checks(args.verify_science)
    if args.verify_study:
        from core.study_checks import run_study_checks
        return run_study_checks(args.verify_study)
    if args.window_test_output:
        return window_check(args.window_test_output)
    from legacy.ui import launch
    launch()
    return 0


if __name__ == "__main__":
    multiprocessing.freeze_support()
    try:
        raise SystemExit(main())
    except Exception:
        from core.runtime import data_root
        logs = data_root() / "logs"
        logs.mkdir(parents=True, exist_ok=True)
        error_file = logs / "legacy-startup-error.log"
        error_file.write_text(traceback.format_exc(), encoding="utf-8")
        if not any(argument.startswith("--") for argument in sys.argv[1:]):
            import ctypes
            ctypes.windll.user32.MessageBoxW(0, "The legacy app could not start. Check its required Windows updates and the log:\n" + str(error_file), TITLE, 0x10)
        if sys.stderr is not None:
            traceback.print_exc()
        raise SystemExit(1)

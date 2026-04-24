"""
Multi-Agent Molecular Docking — Desktop GUI (Skeleton)
======================================================

Author:  Dr. Manu Kumar Shetty, MAMC & Lok Nayak Hospital, New Delhi
Purpose: Thin PySide6 wrapper around the existing MultiAgent_Docking pipeline.

What this does:
    - Shows two input fields (Ligand name, Receptor PDB IDs)
    - A "Run Docking" button
    - A live log panel that captures everything the pipeline prints
    - When finished, opens the generated report folder

What this does NOT do (yet):
    - Package into .exe (add PyInstaller step later)
    - Progress bar with real % (agents don't report % — we show log stream instead)
    - Manage API keys in a settings dialog (key is read from environment for now)

Requirements:
    pip install PySide6
    (plus whatever the notebook already needs: aisuite, rdkit, meeko, vina, etc.)

How to run:
    1. Install and start Ollama (https://ollama.com/download/windows)
           then:  ollama pull llama3.2:3b
    2. Make sure docking_pipeline.py is in the same folder as this file.
    3. py -3 docking_gui.py

Architecture:
    MainWindow (UI thread)  --->  DockingWorker (QThread)  --->  executor_agent()
         ^                              |
         |  log signals (str)           |  print() captured
         +------------------------------+
"""

import os
import sys
from pathlib import Path
from contextlib import redirect_stdout, redirect_stderr

from PySide6.QtCore import QThread, Signal
from PySide6.QtGui import QFont, QTextCursor
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QPlainTextEdit,
    QMessageBox, QGroupBox, QStatusBar,
)


# =============================================================================
# PIPELINE IMPORT
# =============================================================================
# Step 1: Convert your notebook to a .py file ONCE:
#     jupyter nbconvert --to script MultiAgent_Docking.ipynb
#
# Step 2: Rename the resulting file to  docking_pipeline.py  and place it next
#         to this GUI file.
#
# Step 3: In docking_pipeline.py, comment out or remove:
#         - The os.environ["OPENAI_API_KEY"] = "sk-..."  line  (security!)
#         - The last "run the full pipeline" cells at the bottom that call
#           executor_agent() automatically -- we want to import, not run.
#
# Step 4: The GUI below will import executor_agent from it.
# =============================================================================

try:
    # This is the single entry point from your notebook
    from docking_pipeline import executor_agent, BLACKBOARD, PROJECT_ROOT
    PIPELINE_AVAILABLE = True
    PIPELINE_ERROR = None
except Exception as e:
    import traceback
    PIPELINE_AVAILABLE = False
    # Full traceback tells us exactly which line of docking_pipeline.py failed
    PIPELINE_ERROR = traceback.format_exc()
    PROJECT_ROOT = Path.cwd()  # fallback so the GUI still opens


# =============================================================================
# WORKER THREAD -- runs the long docking pipeline off the UI thread
# =============================================================================
class DockingWorker(QThread):
    """Runs executor_agent() in a background thread so the GUI stays responsive.

    Why a thread: executor_agent takes minutes. If we ran it in the UI thread,
    Windows would grey the window out and show "Not Responding".
    """

    log_line = Signal(str)          # emits each line of captured stdout
    finished_ok = Signal(dict)      # emits BLACKBOARD dict when done
    failed = Signal(str)            # emits error message on failure

    def __init__(self, ligand_name, pdb_ids, parent=None):
        super().__init__(parent)
        self.ligand_name = ligand_name
        self.pdb_ids = pdb_ids

    def run(self):
        """Entry point for the QThread. Captures stdout and forwards to UI."""

        # StreamToSignal: file-like object that emits each write() as a signal
        class StreamToSignal:
            def __init__(self, signal):
                self._signal = signal
                self._buffer = ""

            def write(self, text):
                self._buffer += text
                # Emit line-by-line so the UI updates smoothly
                while "\n" in self._buffer:
                    line, self._buffer = self._buffer.split("\n", 1)
                    self._signal.emit(line)

            def flush(self):
                if self._buffer:
                    self._signal.emit(self._buffer)
                    self._buffer = ""

        stream = StreamToSignal(self.log_line)

        try:
            with redirect_stdout(stream), redirect_stderr(stream):
                # THE ONE CALL that runs your entire notebook pipeline
                executor_agent(self.ligand_name, self.pdb_ids)
            stream.flush()
            self.finished_ok.emit(dict(BLACKBOARD))
        except Exception as e:
            stream.flush()
            import traceback
            self.failed.emit(f"{type(e).__name__}: {e}\n\n{traceback.format_exc()}")


# =============================================================================
# MAIN WINDOW
# =============================================================================
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Multi-Agent Molecular Docking  --  Dr Manu Kumar Shetty, Dept of Pharmacology, MAMC")
        self.resize(900, 700)

        self.worker = None

        self._build_ui()
        self._check_environment()

    # ---------- UI construction ----------
    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)

        # --- Input group ---------------------------------------------------
        input_box = QGroupBox("Docking inputs")
        form = QVBoxLayout(input_box)

        # Ligand
        ligand_row = QHBoxLayout()
        ligand_row.addWidget(QLabel("Ligand name:"))
        self.ligand_edit = QLineEdit()
        self.ligand_edit.setPlaceholderText("e.g. Levetiracetam")
        ligand_row.addWidget(self.ligand_edit)
        form.addLayout(ligand_row)

        # Receptors
        receptor_row = QHBoxLayout()
        receptor_row.addWidget(QLabel("Receptor PDB IDs:"))
        self.receptor_edit = QLineEdit()
        self.receptor_edit.setPlaceholderText("comma-separated, e.g. 6CM4, 6A93")
        receptor_row.addWidget(self.receptor_edit)
        form.addLayout(receptor_row)

        # Helper text
        hint = QLabel(
            "Tip: PDB IDs are 4-character codes from rcsb.org. "
            "Ligand name should match a PubChem entry."
        )
        hint.setStyleSheet("color: gray; font-style: italic;")
        form.addWidget(hint)

        root.addWidget(input_box)

        # --- Buttons -------------------------------------------------------
        button_row = QHBoxLayout()
        self.run_btn = QPushButton("\u25B6  Run Docking")
        self.run_btn.setMinimumHeight(36)
        self.run_btn.clicked.connect(self._on_run_clicked)
        button_row.addWidget(self.run_btn)

        self.open_reports_btn = QPushButton("Open Reports Folder")
        self.open_reports_btn.clicked.connect(self._open_reports_folder)
        self.open_reports_btn.setEnabled(False)
        button_row.addWidget(self.open_reports_btn)

        self.clear_log_btn = QPushButton("Clear Log")
        self.clear_log_btn.clicked.connect(lambda: self.log.clear())
        button_row.addWidget(self.clear_log_btn)

        root.addLayout(button_row)

        # --- Log panel -----------------------------------------------------
        log_box = QGroupBox("Pipeline log  (live)")
        log_layout = QVBoxLayout(log_box)
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        # Monospace font -- so ASCII art / tables from the pipeline align
        mono = QFont("Consolas" if sys.platform == "win32" else "Menlo")
        mono.setStyleHint(QFont.Monospace)
        mono.setPointSize(10)
        self.log.setFont(mono)
        self.log.setStyleSheet(
            "QPlainTextEdit { background: #1e1e1e; color: #e0e0e0; "
            "border: 1px solid #444; }"
        )
        log_layout.addWidget(self.log)
        root.addWidget(log_box, stretch=1)

        # --- Status bar ----------------------------------------------------
        self.setStatusBar(QStatusBar())
        self.statusBar().showMessage("Ready.")

    # ---------- Environment / pipeline checks ----------
    def _check_environment(self):
        """Warn the user if anything is wrong before they click Run."""
        if not PIPELINE_AVAILABLE:
            msg = (
                f"docking_pipeline.py could not be imported.\n\n"
                f"Error:\n{PIPELINE_ERROR}\n\n"
                f"Checklist:\n"
                f"  1. Did you run:  jupyter nbconvert --to script MultiAgent_Docking.ipynb\n"
                f"  2. Did you rename the output to docking_pipeline.py?\n"
                f"  3. Is it in the same folder as docking_gui.py?\n"
                f"  4. Did you comment out get_ipython() / display() lines?\n"
                f"  5. Did you comment out the cells that CALL executor_agent(...)?"
            )
            # Log it
            self._log_warn("[!] " + msg)
            # And show a popup so it is impossible to miss
            QMessageBox.critical(self, "Pipeline import failed", msg)
            self.run_btn.setEnabled(False)
            self.statusBar().showMessage("Pipeline not found -- see log.")
            return

        # Check that Ollama is reachable (the pipeline needs the local LLM server)
        try:
            import urllib.request
            urllib.request.urlopen("http://localhost:11434/api/tags", timeout=2).read()
            self._log_info("[ok] Ollama server reachable at http://localhost:11434")
        except Exception:
            self._log_warn(
                "[!] Could not reach Ollama at http://localhost:11434\n"
                "    The pipeline's planner/executor will fail when it calls the LLM.\n"
                "    Make sure Ollama is installed and running:\n"
                "        1. Download from https://ollama.com/download/windows\n"
                "        2. Pull the default model:  ollama pull llama3.2:3b\n"
                "    Ollama usually auto-starts in the background after install.\n"
                "    If not, open a command prompt and run:  ollama serve"
            )

        self._log_info(f"Project root: {PROJECT_ROOT}")
        self._log_info("Ready. Enter ligand + PDB IDs and click Run.\n")

    # ---------- Run button handler ----------
    def _on_run_clicked(self):
        ligand = self.ligand_edit.text().strip()
        pdb_text = self.receptor_edit.text().strip()

        if not ligand:
            QMessageBox.warning(self, "Missing input", "Please enter a ligand name.")
            return
        if not pdb_text:
            QMessageBox.warning(self, "Missing input", "Please enter at least one PDB ID.")
            return

        # Parse comma/space-separated PDB IDs and normalise to uppercase
        pdb_ids = [x.strip().upper() for x in pdb_text.replace(",", " ").split() if x.strip()]

        # Light validation -- PDB IDs are always 4 alphanumeric characters
        bad = [p for p in pdb_ids if len(p) != 4 or not p.isalnum()]
        if bad:
            QMessageBox.warning(
                self, "Invalid PDB ID",
                f"These don't look like valid PDB IDs (must be 4 alphanumeric chars):\n{bad}"
            )
            return

        # Disable UI controls while running
        self.run_btn.setEnabled(False)
        self.ligand_edit.setEnabled(False)
        self.receptor_edit.setEnabled(False)
        self.open_reports_btn.setEnabled(False)
        self.statusBar().showMessage("Running pipeline... this may take several minutes.")

        self._log_info(f"\n{'='*60}")
        self._log_info(f"Starting docking run")
        self._log_info(f"  Ligand:    {ligand}")
        self._log_info(f"  Receptors: {pdb_ids}")
        self._log_info(f"{'='*60}\n")

        # Kick off the worker thread
        self.worker = DockingWorker(ligand, pdb_ids)
        self.worker.log_line.connect(self._append_log)
        self.worker.finished_ok.connect(self._on_finished)
        self.worker.failed.connect(self._on_failed)
        self.worker.start()

    # ---------- Worker callbacks ----------
    def _append_log(self, line):
        self.log.appendPlainText(line)
        self.log.moveCursor(QTextCursor.End)

    def _on_finished(self, blackboard):
        self._log_info("\n[OK] Pipeline finished successfully.")
        report_path = blackboard.get("report")
        if report_path:
            self._log_info(f"   Report: {report_path}")

        self.statusBar().showMessage("Done.")
        self.run_btn.setEnabled(True)
        self.ligand_edit.setEnabled(True)
        self.receptor_edit.setEnabled(True)
        self.open_reports_btn.setEnabled(True)

        QMessageBox.information(
            self, "Docking complete",
            "Pipeline finished successfully.\nClick 'Open Reports Folder' to view the output."
        )

    def _on_failed(self, err):
        self._log_warn(f"\n[FAIL] Pipeline failed:\n{err}")
        self.statusBar().showMessage("Failed -- see log.")
        self.run_btn.setEnabled(True)
        self.ligand_edit.setEnabled(True)
        self.receptor_edit.setEnabled(True)
        QMessageBox.critical(self, "Docking failed", err.splitlines()[0])

    # ---------- Helpers ----------
    def _open_reports_folder(self):
        """Open the reports folder in the OS file manager."""
        reports = Path(PROJECT_ROOT) / "reports"
        if not reports.exists():
            QMessageBox.information(self, "No reports yet", f"{reports} does not exist.")
            return

        if sys.platform == "win32":
            os.startfile(reports)  # type: ignore[attr-defined]
        elif sys.platform == "darwin":
            os.system(f'open "{reports}"')
        else:
            os.system(f'xdg-open "{reports}"')

    def _log_info(self, msg):
        self.log.appendPlainText(msg)
        self.log.moveCursor(QTextCursor.End)

    def _log_warn(self, msg):
        # Could style this differently later; for now just append
        self.log.appendPlainText(msg)
        self.log.moveCursor(QTextCursor.End)


# =============================================================================
# ENTRY POINT
# =============================================================================
def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()

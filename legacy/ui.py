"""Native Tk desktop for the legacy edition; no browser runtime is required.

The scientific backend is shared with the modern edition. All long operations
run on worker threads and return to Tk through a queue. The optional 3D view
opens an exported local HTML file in the user's browser.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import queue
import shutil
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import webbrowser


TITLE = "In-Silico Docking Studio"
DEVELOPERS = (
    "Dr Manu Kumar Shetty, Professor, Department of Pharmacology, Maulana Azad Medical College, New Delhi",
    "Dr Pratibha Nadig, Founder, Akshaisarathi Healthcare",
)
STEPS = ("Protein targets", "Drug molecules", "Prepare and inspect", "Run docking", "Results and export")
STEP_NOTES = (
    "Add one or more protein structures. Each protein will keep its own binding region.",
    "Add one or more drug molecules. Every drug will be docked against every target.",
    "Prepare the structures, review the recorded changes and inspect each binding region.",
    "Run the complete target–drug matrix. Calculation progress appears at the bottom of the window.",
    "Read each target's results and scientific checks, then save a CSV or complete study ZIP.",
)
INK, MUTED, TEAL, BG, LINE = "#17344b", "#5d7282", "#087c80", "#f3f7f9", "#dce7ec"
SITE_LABELS = {"Reference ligand": "ligand", "Detect a pocket": "detect", "Manual box": "manual"}
BOX_KEYS = ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z")


def status_label(value):
    text = str(value or "")
    return {"pass": "Checked", "warn": "Warning", "fail": "Check failed",
            "not_applicable": "Not applicable", "not_checked": "Not checked",
            "stale": "Needs updating"}.get(text, text.replace("_", " ").capitalize())


class DockingStudio:
    """One native application window; public navigation methods support smoke checks."""

    def __init__(self, root, project_root=None):
        from core import study, references, project_library
        from core.runtime import data_root, resource_root
        from core.settings import projects_root
        self.root, self.engine, self.references, self.library = root, study, references, project_library
        self.resource_root, self.data_root = resource_root(), data_root()
        self.project_root = Path(project_root) if project_root else projects_root()
        self.project_root.mkdir(parents=True, exist_ok=True)
        self.active = None
        self.current, self.section, self.workflow, self.busy = 0, "Home", False, False
        self.draft = {"name": "Docking project", "targets": [], "drugs": []}
        self.library_rows, self.result_rows, self.events = [], [], queue.Queue()
        self.draft_started = False
        self.status = tk.StringVar(value="Ready. Start a new project or open saved work.")
        self.name = tk.StringVar(value=self.draft["name"])
        self.target_name, self.pdb_id, self.chains = tk.StringVar(), tk.StringVar(), tk.StringVar()
        self.site_mode = tk.StringVar(value="Reference ligand")
        self.box_vars = {key: tk.StringVar(value="20" if key.startswith("size") else "0") for key in BOX_KEYS}
        self.smiles_input, self.lookup_input = "", ""
        self.target_choice, self.drug_choice = tk.StringVar(), tk.StringVar()
        self.effort, self.comparisons = tk.StringVar(value="8"), tk.StringVar(value="3")
        self.logo = None
        root.title(TITLE + " — Legacy Beta")
        # Older classroom PCs often have a 1024x768 desktop. Keep the first
        # window and its progress footer inside the available screen.
        width = min(1180, max(780, root.winfo_screenwidth() - 50))
        height = min(820, max(540, root.winfo_screenheight() - 100))
        root.geometry("{}x{}".format(width, height))
        root.minsize(min(940, width), min(620, height))
        root.configure(background=BG)
        root.protocol("WM_DELETE_WINDOW", self.close)
        self._style()
        self._shell()
        self.render()
        root.after(100, self.poll_queue)

    def _style(self):
        style = ttk.Style(self.root)
        style.theme_use("clam")
        style.configure("TFrame", background=BG)
        style.configure("Card.TFrame", background="white")
        style.configure("TLabel", background=BG, foreground=INK, font=("Segoe UI", 10))
        style.configure("Card.TLabel", background="white", foreground=INK, font=("Segoe UI", 10))
        style.configure("Title.TLabel", background=BG, foreground=INK, font=("Segoe UI", 21, "bold"))
        style.configure("Heading.TLabel", background="white", foreground=INK, font=("Segoe UI", 13, "bold"))
        style.configure("Muted.TLabel", background=BG, foreground=MUTED, font=("Segoe UI", 10))
        style.configure("CardMuted.TLabel", background="white", foreground=MUTED, font=("Segoe UI", 10))
        style.configure("TButton", font=("Segoe UI", 10), padding=(12, 8), background="white", foreground=TEAL)
        style.configure("Primary.TButton", background=TEAL, foreground="white", padding=(17, 10))
        style.map("Primary.TButton", background=[("active", "#096a70"), ("disabled", "#dce7ec")])
        style.configure("Danger.TButton", foreground="#b63d48")
        style.configure("Nav.TButton", padding=(12, 9), background="white")
        style.configure("Selected.Nav.TButton", background="#d5e9e9", foreground=TEAL)
        style.configure("TEntry", padding=7, fieldbackground="white")
        style.configure("TCombobox", padding=7, fieldbackground="white")
        style.configure("Treeview", font=("Segoe UI", 10), rowheight=29, background="white", fieldbackground="white")
        style.configure("Treeview.Heading", font=("Segoe UI", 10, "bold"), background="#edf5f6", foreground=INK)
        style.map("Treeview", background=[("selected", "#d5e9e9")], foreground=[("selected", INK)])

    def _shell(self):
        header = ttk.Frame(self.root, style="Card.TFrame", padding=(18, 12))
        header.pack(fill="x")
        try:
            from PIL import Image, ImageTk
            image = Image.open(self.resource_root / "assets/logo.png").convert("RGBA")
            image.thumbnail((175, 48))
            self.logo = ImageTk.PhotoImage(image)
            ttk.Label(header, image=self.logo, style="Card.TLabel").pack(side="left", padx=(0, 18))
        except Exception:
            ttk.Label(header, text="AKSHAISARATHI", style="Heading.TLabel").pack(side="left", padx=(0, 18))
        brand = ttk.Frame(header, style="Card.TFrame")
        brand.pack(side="left")
        ttk.Label(brand, text=TITLE, style="Heading.TLabel").pack(anchor="w")
        ttk.Label(brand, text="MEDICAL EDUCATION · LEGACY BETA", style="CardMuted.TLabel").pack(anchor="w")
        self.main_buttons = {}
        for section in ("About", "Reference", "Projects", "Home"):
            button = ttk.Button(header, text=section, style="Nav.TButton", command=lambda value=section: self.show_section(value))
            button.pack(side="right", padx=2)
            self.main_buttons[section] = button
        ttk.Separator(self.root).pack(fill="x")
        self.main = ttk.Frame(self.root)
        self.main.pack(fill="both", expand=True)
        self.sidebar = ttk.Frame(self.main, padding=14, width=215)
        self.canvas = tk.Canvas(self.main, background=BG, highlightthickness=0)
        self.scrollbar = ttk.Scrollbar(self.main, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.scrollbar.set)
        self.scrollbar.pack(side="right", fill="y")
        self.canvas.pack(side="right", fill="both", expand=True)
        self.body = ttk.Frame(self.canvas, padding=(26, 25))
        self.body_id = self.canvas.create_window((0, 0), window=self.body, anchor="nw")
        self.body.bind("<Configure>", lambda event: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self.canvas.bind("<Configure>", lambda event: self.canvas.itemconfigure(self.body_id, width=event.width))
        self.root.bind("<MouseWheel>", self._wheel, add="+")
        footer = ttk.Frame(self.root, padding=(18, 9))
        footer.pack(fill="x")
        self.progress_bar = ttk.Progressbar(footer, mode="indeterminate", length=115)
        self.progress_bar.pack(side="right", padx=(12, 0))
        ttk.Label(footer, textvariable=self.status, style="Muted.TLabel", wraplength=850).pack(side="left", fill="x", expand=True)

    def _wheel(self, event):
        if event.widget.winfo_class() not in {"Text", "Treeview", "TCombobox"}:
            self.canvas.yview_scroll(-int(event.delta / 120), "units")

    def _text(self, parent, text, style="Card.TLabel", wrap=790):
        label = ttk.Label(parent, text=text, style=style, wraplength=wrap, justify="left")
        label.pack(anchor="w", fill="x", pady=(0, 10))
        return label

    def _card(self, title=""):
        card = ttk.Frame(self.body, style="Card.TFrame", padding=20)
        card.pack(fill="x", pady=(0, 16))
        if title:
            self._text(card, title, "Heading.TLabel")
        return card

    def _button(self, parent, label, command, primary=False, enabled=True, side="left"):
        button = ttk.Button(parent, text=label, command=self._guard(command), style="Primary.TButton" if primary else "TButton")
        button.pack(side=side, padx=(0, 9), pady=4)
        if not enabled or self.busy:
            button.state(["disabled"])
        return button

    def _guard(self, function):
        def invoke():
            try:
                return function()
            except Exception as error:
                self.status.set(str(error))
                messagebox.showerror(TITLE, str(error), parent=self.root)
        return invoke

    def _entry(self, parent, label, variable, width=None):
        frame = ttk.Frame(parent, style="Card.TFrame")
        frame.pack(side="left" if width else "top", fill="x" if not width else None, expand=not width, padx=(0, 12), pady=5)
        ttk.Label(frame, text=label, style="CardMuted.TLabel").pack(anchor="w", pady=(0, 4))
        entry = ttk.Entry(frame, textvariable=variable, width=width)
        entry.pack(fill="x")
        return entry

    def _tree(self, parent, headings, rows, height=5, widths=None):
        wrap = ttk.Frame(parent, style="Card.TFrame")
        wrap.pack(fill="both", expand=True, pady=(6, 9))
        names = [str(index) for index in range(len(headings))]
        tree = ttk.Treeview(wrap, columns=names, show="headings", height=height, selectmode="browse")
        for index, heading in enumerate(headings):
            tree.heading(str(index), text=heading)
            tree.column(str(index), width=(widths[index] if widths else 170), minwidth=70, stretch=True)
        for index, values in enumerate(rows):
            tree.insert("", "end", iid=str(index), values=values)
        vertical = ttk.Scrollbar(wrap, orient="vertical", command=tree.yview)
        horizontal = ttk.Scrollbar(wrap, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=vertical.set, xscrollcommand=horizontal.set)
        tree.grid(row=0, column=0, sticky="nsew")
        vertical.grid(row=0, column=1, sticky="ns")
        horizontal.grid(row=1, column=0, sticky="ew")
        wrap.rowconfigure(0, weight=1)
        wrap.columnconfigure(0, weight=1)
        return tree

    def _save_text_fields(self):
        if getattr(self, "smiles_widget", None) is not None and self.smiles_widget.winfo_exists():
            self.smiles_input = self.smiles_widget.get("1.0", "end-1c")
        if getattr(self, "lookup_widget", None) is not None and self.lookup_widget.winfo_exists():
            self.lookup_input = self.lookup_widget.get("1.0", "end-1c")
        self.draft["name"] = self.name.get().strip() or "Docking project"

    def show_section(self, section):
        self._save_text_fields()
        self.section, self.workflow = section, False
        self.render()
        if section == "Projects" and not self.busy:
            self.run_worker("Reading project sizes…", lambda: self.library.list_projects(self.project_root), self._library_loaded)

    def show_workflow(self, index):
        if self.busy or (index > 1 and self.active is None):
            return
        self._save_text_fields()
        self.section, self.workflow, self.current = "Home", True, index
        if self.active is not None:
            self.active.data["ui_step"] = index
            self.active.save()
        self.render()

    def render(self):
        self._save_text_fields()
        for child in self.body.winfo_children():
            child.destroy()
        for child in self.sidebar.winfo_children():
            child.destroy()
        self.smiles_widget = self.lookup_widget = None
        for section, button in self.main_buttons.items():
            button.configure(style="Selected.Nav.TButton" if section == self.section else "Nav.TButton")
        if self.workflow:
            self.sidebar.pack(side="left", fill="y", before=self.canvas)
            ttk.Label(self.sidebar, text="DOCKING STEPS", foreground=MUTED).pack(anchor="w", pady=(4, 12))
            for index, step in enumerate(STEPS):
                button = ttk.Button(self.sidebar, text="{}. {}".format(index + 1, step), command=lambda value=index: self.show_workflow(value), style="Selected.Nav.TButton" if index == self.current else "Nav.TButton")
                button.pack(fill="x", pady=4)
                if self.busy or (index > 1 and self.active is None):
                    button.state(["disabled"])
            ttk.Separator(self.sidebar).pack(fill="x", pady=17)
            ttk.Label(self.sidebar, text=self.active.data.get("name", "Docking project") if self.active else self.draft["name"], wraplength=180).pack(anchor="w")
            ttk.Label(self.sidebar, text="Saved locally" if self.active else "New project inputs", foreground=MUTED).pack(anchor="w", pady=(5, 12))
            ttk.Button(self.sidebar, text="Back to Home", command=lambda: self.show_section("Home")).pack(fill="x")
            self._text(self.body, "{} · Step {} of 5".format(self.active.data.get("name", "Docking project") if self.active else self.draft["name"], self.current + 1), "Muted.TLabel")
            self._text(self.body, STEPS[self.current], "Title.TLabel")
            self._text(self.body, STEP_NOTES[self.current], "Muted.TLabel")
            (self.render_targets, self.render_drugs, self.render_prepare, self.render_dock, self.render_results)[self.current]()
        else:
            self.sidebar.pack_forget()
            self._text(self.body, self.section, "Title.TLabel")
            {"Home": self.render_home, "Projects": self.render_projects, "Reference": self.render_reference, "About": self.render_about}[self.section]()
        self.canvas.yview_moveto(0)

    def run_worker(self, label, function, on_done=None):
        if self.busy:
            return
        self.busy = True
        self.status.set(label)
        self.progress_bar.start(12)
        self.render()
        def work():
            try:
                self.events.put(("done", function(), on_done))
            except Exception as error:
                self.events.put(("error", str(error), None))
        threading.Thread(target=work, daemon=True).start()

    def poll_queue(self):
        try:
            while True:
                kind, value, callback = self.events.get_nowait()
                if kind == "progress":
                    self.status.set(str(value))
                else:
                    self.busy = False
                    self.progress_bar.stop()
                    if kind == "error":
                        self.status.set(str(value))
                        messagebox.showerror(TITLE, str(value), parent=self.root)
                    else:
                        self.status.set("Completed. Review the recorded results and checks.")
                        if callback:
                            self._guard(lambda: callback(value))()
                    self.render()
        except queue.Empty:
            pass
        if self.root.winfo_exists():
            self.root.after(100, self.poll_queue)

    def progress(self, text):
        self.events.put(("progress", str(text), None))

    def close(self):
        if self.busy:
            messagebox.showinfo(TITLE, "A calculation is still running. Let it finish before closing this window.", parent=self.root)
            return
        self.root.destroy()

    def render_home(self):
        self._text(self.body, "Start a new docking project, or continue saved work in Projects.", "Muted.TLabel")
        card = self._card("New drug and target docking")
        self._text(card, "Choose proteins and drug molecules, prepare the structures, run local docking and review the results.")
        self._button(card, "Start new project", self.new_project, True)
        flow = self._card()
        self._text(flow, "1. Targets  →  2. Drugs  →  3. Prepare  →  4. Dock  →  5. Results", "CardMuted.TLabel")
        if self.active is not None or self.draft_started:
            card = self._card(self.active.data.get("name", "Docking project") if self.active else self.draft["name"])
            self._text(card, "Step {} of 5 · {}".format(self.current + 1, STEPS[self.current]), "CardMuted.TLabel")
            self._button(card, "Continue project", lambda: self.show_workflow(self.current))
        card = self._card("Try a ready example")
        self._text(card, "Use the bundled EGFR protein and example drug inputs. New docking results are produced when you run docking.")
        self._button(card, "Use the offline EGFR example", lambda: self.run_worker("Creating the example project…", lambda: self.engine.create_example_study(self.project_root), self._open_created))
        self._text(self.body, "Legacy beta: native controls do not require WebView2. Optional interactive 3D opens in your installed browser.", "Muted.TLabel")

    def new_project(self):
        if self.busy:
            return
        for widget in (getattr(self, "smiles_widget", None), getattr(self, "lookup_widget", None)):
            if widget is not None and widget.winfo_exists():
                widget.delete("1.0", "end")
        self.active = None
        self.draft = {"name": "Docking project", "targets": [], "drugs": []}
        self.name.set("Docking project")
        self.smiles_input = self.lookup_input = ""
        self.current, self.section, self.workflow, self.draft_started = 0, "Home", True, True
        self.render()

    def _open_created(self, study):
        self.active = study
        self.current, self.section, self.workflow = 2, "Home", True
        study.data["ui_step"] = 2
        study.save()

    def _library_loaded(self, rows):
        self.library_rows = rows
        self.status.set("{} saved projects · {} on disk".format(len(rows), self.library.human_size(sum(row["size_bytes"] for row in rows))))

    def render_projects(self):
        self._text(self.body, "Open saved work or delete a project you no longer need. Saved files use disk space; structures are loaded when opened.", "Muted.TLabel")
        top = self._card()
        self._text(top, "Saved in: " + str(self.project_root), "CardMuted.TLabel")
        self._text(top, "{} saved projects · {} on disk".format(len(self.library_rows), self.library.human_size(sum(row["size_bytes"] for row in self.library_rows))))
        tree = self._tree(top, ("Project", "Status", "Disk size", "Targets", "Drugs"), [(row["name"], row.get("status", "saved"), self.library.human_size(row["size_bytes"]), row.get("target_count", 0), row.get("drug_count", 0)) for row in self.library_rows], 10, (280, 160, 100, 75, 75))
        actions = ttk.Frame(top, style="Card.TFrame")
        actions.pack(fill="x")
        def selected():
            if not tree.selection():
                raise ValueError("Select a saved project in the list.")
            return self.library_rows[int(tree.selection()[0])]
        self._button(actions, "Open project", lambda: self.open_project(selected()), True)
        self._button(actions, "Delete project", lambda: self.delete_project(selected()))
        self._button(actions, "Refresh projects", lambda: self.run_worker("Reading project sizes…", lambda: self.library.list_projects(self.project_root), self._library_loaded))
        self._button(actions, "Start new project", self.new_project)
        if self.busy:
            self._text(self.body, "Project management is available after the current calculation finishes.", "Muted.TLabel")

    def open_project(self, item):
        if self.busy:
            return
        if item["kind"] != "study":
            self._legacy_project(item)
            return
        study = self.engine.open_study(item["folder"])
        self.active = study
        self.current = max(0, min(4, int(study.data.get("ui_step", 2))))
        self.section, self.workflow = "Home", True
        self.effort.set(str(study.data.get("settings", {}).get("exhaustiveness", 8)))
        self.comparisons.set(str(study.data.get("settings", {}).get("n_decoys", 3)))
        self.render()

    def _legacy_project(self, item):
        from core.project import Project
        from core.references import render_references
        project = Project.open(item["folder"])
        try:
            lines = [project.name, "", "Previously saved individual project", ""]
            for row in project.conn.execute("SELECT stage_id,title,status,summary,control_json FROM steps ORDER BY stage_id"):
                lines.extend(["{} · {} · {}".format(row[0], row[1], row[2]), row[3] or ""])
                if row[4]:
                    control = json.loads(row[4])
                    lines.append("Scientific check: {} · {}".format(control.get("verdict", "not recorded"), control.get("plain_verdict", "")))
                lines.append("")
            lines.extend(["", render_references([project])])
        finally:
            project.close()
        window = tk.Toplevel(self.root)
        window.title(TITLE + " — " + item["name"])
        window.geometry("860x600")
        self._readonly_text(window, "\n".join(lines), height=26)
        ttk.Button(window, text="Close", command=window.destroy).pack(pady=12)

    def delete_project(self, item):
        if self.busy:
            return
        description = "Permanently delete '{}' and its inputs, structures and results ({}).\n\nExport any results you need before deleting. Other projects will be kept.".format(item["name"], self.library.human_size(item["size_bytes"]))
        if not messagebox.askokcancel("Delete project", description, icon="warning", parent=self.root):
            return
        def remove():
            self.library.delete_project(self.project_root, item["folder"])
            return self.library.list_projects(self.project_root)
        def done(rows):
            if self.active and self.active.folder.resolve() == Path(item["folder"]).resolve():
                self.active = None
                self.current, self.draft_started = 0, False
            self._library_loaded(rows)
            self.status.set("Project deleted. Other projects were kept.")
        self.run_worker("Deleting the selected project…", remove, done)

    def render_targets(self):
        if self.active:
            self._text(self.body, "These saved inputs belong to this project. Start a new project from Home to use different targets or drugs.", "Muted.TLabel")
            self._target_table(self.active.data["targets"])
            self._button(self.body, "Continue to drug molecules", lambda: self.show_workflow(1), True)
            return
        card = self._card("Project name")
        self._entry(card, "Name", self.name)
        card = self._card("Add protein structures")
        row = ttk.Frame(card, style="Card.TFrame")
        row.pack(fill="x")
        self._entry(row, "PDB identifier", self.pdb_id, 16)
        self._entry(row, "Target label", self.target_name, 27)
        self._entry(row, "Chains, optional", self.chains, 20)
        ttk.Label(card, text="Binding region", style="CardMuted.TLabel").pack(anchor="w", pady=(10, 4))
        selector = ttk.Combobox(card, textvariable=self.site_mode, values=list(SITE_LABELS), state="readonly", width=35)
        selector.pack(anchor="w")
        self._text(card, "Reference ligand uses a bound molecule in the structure. Detect a pocket finds a candidate region. Manual box uses coordinates you provide.", "CardMuted.TLabel")
        box_fields = self._box_fields(card, self.box_vars)
        def toggle_box(event=None):
            if self.site_mode.get() == "Manual box":
                box_fields.pack(fill="x", pady=5)
            else:
                box_fields.pack_forget()
        selector.bind("<<ComboboxSelected>>", toggle_box)
        toggle_box()
        actions = ttk.Frame(card, style="Card.TFrame")
        actions.pack(fill="x", pady=8)
        self._button(actions, "Add protein target", self.add_target, True)
        self._button(actions, "Add local PDB or CIF files", self.add_protein_files)
        self._button(actions, "Protein downloads · RCSB PDB", lambda: webbrowser.open("https://www.rcsb.org/"))
        card = self._card("Targets in this project")
        tree = self._target_table(self.draft["targets"], card)
        actions = ttk.Frame(card, style="Card.TFrame")
        actions.pack(fill="x")
        def remove():
            if tree.selection():
                self.draft["targets"].pop(int(tree.selection()[0]))
                self.render()
        self._button(actions, "Remove selected target", remove)
        self._button(actions, "Continue to drug molecules", lambda: self.show_workflow(1), True, bool(self.draft["targets"]))

    def _box_fields(self, parent, variables):
        outer = ttk.Frame(parent, style="Card.TFrame")
        outer.pack(fill="x", pady=5)
        ttk.Label(outer, text="Manual box coordinates and side lengths · angstroms", style="CardMuted.TLabel").pack(anchor="w")
        row = ttk.Frame(outer, style="Card.TFrame")
        row.pack(fill="x", pady=5)
        for key in BOX_KEYS:
            self._entry(row, key.replace("_", " ").capitalize(), variables[key], 9)
        return outer

    def _box(self, variables):
        result = {key: float(variable.get()) for key, variable in variables.items()}
        if any(not math.isfinite(value) for value in result.values()):
            raise ValueError("Box coordinates and lengths must be finite numbers.")
        if any(not 1 <= result[key] <= 100 for key in BOX_KEYS[3:]):
            raise ValueError("Each box side must be from 1 to 100 angstroms.")
        return result

    def add_target(self, file=""):
        mode = SITE_LABELS[self.site_mode.get()]
        pdb = self.pdb_id.get().strip().upper() if not file else ""
        name = self.target_name.get().strip() or pdb or Path(file).stem
        if not name or not (pdb or file):
            raise ValueError("Enter a protein PDB identifier or add a local protein file.")
        self.draft["targets"].append({"name": name, "pdb_id": pdb, "file": file,
                                      "chains": [value.strip() for value in self.chains.get().replace(";", ",").split(",") if value.strip()] or None,
                                      "site_mode": mode, "box": self._box(self.box_vars) if mode == "manual" else None})
        self.pdb_id.set("")
        self.target_name.set("")
        self.render()

    def add_protein_files(self):
        paths = filedialog.askopenfilenames(parent=self.root, title="Choose protein structures", filetypes=[("Protein structures", "*.pdb *.cif *.mmcif"), ("All files", "*.*")])
        for filename in paths:
            self.add_target(str(Path(filename).resolve()))

    def _target_table(self, targets, parent=None):
        return self._tree(parent or self.body, ("Protein target", "Source", "Status", "Preparation / binding region"), [(row["name"], row.get("pdb_id") or Path(row.get("file", "")).name, row.get("status", "Input added"), row.get("error") or row.get("site_error") or row.get("site_summary") or row.get("preparation_summary") or "Not prepared") for row in targets], max(2, min(6, len(targets))), (160, 120, 100, 360))

    def render_drugs(self):
        if self.active:
            self._text(self.body, "Saved drug identities are kept with their source records. Start a new project to use different inputs.", "Muted.TLabel")
            self._drug_table(self.active.data["drugs"], self._card("Drug molecules"))
            self._button(self.body, "Continue to preparation", lambda: self.show_workflow(2), True)
            return
        card = self._card("Find drugs on PubChem")
        self._text(card, "Enter one drug name or PubChem CID per line. The original chemical identity and source links are saved.", "CardMuted.TLabel")
        self.lookup_widget = tk.Text(card, height=3, font=("Segoe UI", 10), relief="solid", borderwidth=1)
        self.lookup_widget.pack(fill="x", pady=(0, 9))
        self.lookup_widget.insert("1.0", self.lookup_input)
        actions = ttk.Frame(card, style="Card.TFrame")
        actions.pack(fill="x")
        self._button(actions, "Find and add drugs", self.lookup_drugs, True)
        self._button(actions, "Drug downloads · PubChem", lambda: webbrowser.open("https://pubchem.ncbi.nlm.nih.gov/"))
        card = self._card("Add molecules from SMILES or an SDF file")
        self._text(card, "Enter one SMILES followed by a short drug name per line. Example: CCO Ethanol", "CardMuted.TLabel")
        self.smiles_widget = tk.Text(card, height=3, font=("Consolas", 10), relief="solid", borderwidth=1)
        self.smiles_widget.pack(fill="x", pady=(0, 9))
        self.smiles_widget.insert("1.0", self.smiles_input)
        actions = ttk.Frame(card, style="Card.TFrame")
        actions.pack(fill="x")
        self._button(actions, "Add SMILES molecules", self.add_smiles)
        self._button(actions, "Import SDF", self.import_sdf)
        card = self._card("Drugs in this project")
        tree = self._drug_table(self.draft["drugs"], card)
        actions = ttk.Frame(card, style="Card.TFrame")
        actions.pack(fill="x")
        def remove():
            if tree.selection():
                self.draft["drugs"].pop(int(tree.selection()[0]))
                self.render()
        self._button(actions, "Remove selected drug", remove)
        self._button(actions, "Save project and prepare structures", self.create_project, True, bool(self.draft["targets"] and self.draft["drugs"]))

    def _drug_table(self, drugs, parent):
        tree = self._tree(parent, ("Drug", "SMILES", "Source"), [(row["name"], row["smiles"], row.get("source", "User-supplied structure")) for row in drugs], max(2, min(6, len(drugs))), (160, 370, 190))
        caption = ttk.Label(parent, text="Select a drug for its 2D input structure.", style="CardMuted.TLabel")
        caption.pack(anchor="w", pady=(2, 5))
        picture = ttk.Label(parent, style="Card.TLabel")
        picture.pack(anchor="w", pady=(0, 8))
        def show_structure(event=None):
            if not tree.selection():
                return
            row = drugs[int(tree.selection()[0])]
            try:
                from rdkit import Chem
                from rdkit.Chem import Draw
                from PIL import ImageTk
                molecule = Chem.MolFromSmiles(row["smiles"])
                image = Draw.MolToImage(molecule, size=(320, 145))
                photo = ImageTk.PhotoImage(image, master=self.root)
                picture.configure(image=photo)
                picture.image = photo
                caption.configure(text="{} · 2D input structure".format(row["name"]))
            except Exception:
                caption.configure(text="The molecule's name and SMILES identity are shown in the table.")
        tree.bind("<<TreeviewSelect>>", show_structure)
        if drugs:
            tree.selection_set("0")
            show_structure()
        return tree

    def _append_drugs(self, proposed):
        from rdkit import Chem
        names = {row["name"].casefold() for row in self.draft["drugs"]}
        for row in proposed:
            if not row.get("name") or Chem.MolFromSmiles(row.get("smiles", "")) is None:
                raise ValueError("Each drug needs a name and a valid SMILES structure.")
            if row["name"].casefold() in names:
                raise ValueError("A drug named '{}' is already listed. Use a distinct name for a different molecule.".format(row["name"]))
            names.add(row["name"].casefold())
        self.draft["drugs"].extend(proposed)

    def add_smiles(self):
        self._save_text_fields()
        proposed = []
        for index, line in enumerate(self.smiles_input.splitlines(), 1):
            parts = line.strip().split(None, 1)
            if parts:
                proposed.append({"smiles": parts[0], "name": parts[1].strip() if len(parts) > 1 else "Drug {}".format(len(self.draft["drugs"]) + index), "source": "User-supplied SMILES"})
        self._append_drugs(proposed)
        self.smiles_input = ""
        if self.smiles_widget:
            self.smiles_widget.delete("1.0", "end")
        self.render()

    def lookup_drugs(self):
        self._save_text_fields()
        queries = [line.strip() for line in self.lookup_input.splitlines() if line.strip()]
        if not queries:
            raise ValueError("Enter a drug name or PubChem CID.")
        cache = self.data_root / "PubChem"
        def resolve():
            rows, errors = [], []
            for query_text in queries:
                self.progress("Finding {} on PubChem…".format(query_text))
                try:
                    rows.append(self.references.resolve_drug(query_text, cache))
                except Exception as error:
                    errors.append("{}: {}".format(query_text, error))
            return rows, errors
        def done(result):
            rows, errors = result
            self._append_drugs(rows)
            self.lookup_input = ""
            if self.lookup_widget is not None and self.lookup_widget.winfo_exists():
                self.lookup_widget.delete("1.0", "end")
            if errors:
                messagebox.showwarning("Drug lookup", "\n".join(errors), parent=self.root)
            self.status.set("{} drug records added.".format(len(rows)))
        self.run_worker("Finding drug structures…", resolve, done)

    def import_sdf(self):
        from rdkit import Chem
        filename = filedialog.askopenfilename(parent=self.root, title="Choose drug structures", filetypes=[("SDF molecules", "*.sdf"), ("All files", "*.*")])
        if not filename:
            return
        proposed = []
        failed = 0
        for index, molecule in enumerate(Chem.SDMolSupplier(filename, removeHs=False), 1):
            if molecule is None:
                failed += 1
                continue
            proposed.append({"name": molecule.GetProp("_Name").strip() if molecule.HasProp("_Name") and molecule.GetProp("_Name").strip() else "{} {}".format(Path(filename).stem, index),
                             "smiles": Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True), "source": "Imported SDF"})
        self._append_drugs(proposed)
        self.render()
        self.status.set("{} molecules imported; {} unreadable records.".format(len(proposed), failed))

    def create_project(self):
        if not self.draft["targets"] or not self.draft["drugs"]:
            raise ValueError("Add at least one target and one drug.")
        self._save_text_fields()
        targets = [self.engine.TargetSpec(**dict(row)) for row in self.draft["targets"]]
        drugs = [self.engine.DrugSpec(**{key: row.get(key, "") for key in ("name", "smiles", "cid", "source_url", "source", "retrieved")}) for row in self.draft["drugs"]]
        name = self.draft["name"]
        self.run_worker("Saving the project inputs…", lambda: self.engine.create_study(self.project_root, name, targets, drugs), self._open_created)

    def _summary(self):
        card = self._card()
        targets, drugs = self.active.data["targets"], self.active.data["drugs"]
        self._text(card, "{} targets × {} drugs = {} docking combinations".format(len(targets), len(drugs), len(targets) * len(drugs)), "Heading.TLabel")
        self._text(card, "Each target has a separate structure, search box, result set and scientific checks.", "CardMuted.TLabel")

    def _ready_to_dock(self):
        return any(target.get("status") in {"prepared", "docked", "partial"} and target.get("site_status") not in {"needs_input", "fail"} for target in self.active.data["targets"])

    def render_prepare(self):
        self._summary()
        card = self._card("Prepare proteins and drugs")
        self._text(card, "Protein preparation repairs the model and adds hydrogens. Drug preparation produces 3D structures and docking inputs. Review warnings before proceeding.")
        self._button(card, "Prepare proteins and drugs", lambda: self.run_worker("Preparing proteins and drugs…", lambda: self.engine.run_preparation(self.active, progress=self.progress)), True)
        self._target_table(self.active.data["targets"], card)
        self._inspection("04_site")
        self._button(self.body, "Continue to docking", lambda: self.show_workflow(3), True, self._ready_to_dock())

    def _selected_target(self):
        targets = self.active.data["targets"]
        names = [row["name"] + " · " + row["id"] for row in targets]
        if self.target_choice.get() not in names:
            self.target_choice.set(names[0] if names else "")
        return targets[names.index(self.target_choice.get())] if names else None

    def _inspection(self, stage):
        card = self._card("Inspect a target" if stage == "04_site" else "Inspect a docked pose")
        self._text(card, "Choose a target to read its preparation and scientific checks. Interactive 3D is optional and uses a local file in your installed browser.", "CardMuted.TLabel")
        target = self._selected_target()
        choices = [row["name"] + " · " + row["id"] for row in self.active.data["targets"]]
        selector = ttk.Combobox(card, textvariable=self.target_choice, values=choices, state="readonly", width=48)
        selector.pack(anchor="w", pady=8)
        selector.bind("<<ComboboxSelected>>", lambda event: self.render())
        if not target:
            return
        project = self.engine.target_project(self.active, target["id"])
        try:
            from core import viewer
            if stage == "06_dock":
                docked = viewer.values(project, "06_dock").get("results", [])
                names = {row.get("internal_name", row["name"]): row["name"] for row in self.active.data["drugs"]}
                options = [(row["name"], names.get(row["name"], row["name"])) for row in docked if row.get("file")]
                display = [label + " · " + key for key, label in options]
                if self.drug_choice.get() not in display:
                    self.drug_choice.set(display[0] if display else "")
                if display:
                    ttk.Combobox(card, textvariable=self.drug_choice, values=display, state="readonly", width=48).pack(anchor="w", pady=8)
            self._button(card, "Open local 3D view", lambda: self.open_view(stage), enabled=not self.busy)
            texts = []
            for sid in (["03_receptor", "04_site"] if stage == "04_site" else ["06_dock"]):
                saved = project.get_step(sid)
                if saved:
                    texts.extend([saved["title"] + " · " + saved["status"], saved["summary"] or ""])
                    control = json.loads(saved["control_json"] or "{}")
                    if control:
                        texts.append("Scientific check: {} · {}".format(control.get("verdict", "not recorded"), control.get("plain_verdict", "")))
                        texts.append(control.get("action", ""))
                    data = json.loads(saved["values_json"] or "{}")
                    for key, check in (data.get("checks") or {}).items():
                        if isinstance(check, dict):
                            texts.append("{}: {} · {}".format(key.capitalize(), status_label(check.get("verdict", "not recorded")), check.get("plain_verdict", "")))
                    texts.append("")
            self._readonly_text(card, "\n".join(texts) or "The preparation and check record appears after this step runs.", height=7)
        finally:
            project.close()
        if stage == "04_site":
            self._site_editor(card, target)

    def _site_editor(self, card, target):
        self._text(card, "Binding region for this target", "Heading.TLabel")
        current_mode = next((label for label, value in SITE_LABELS.items() if value == target.get("site_mode", "ligand")), "Reference ligand")
        mode = tk.StringVar(value=current_mode)
        selector = ttk.Combobox(card, textvariable=mode, values=list(SITE_LABELS), state="readonly", width=35)
        selector.pack(anchor="w", pady=6)
        box = target.get("box") or {}
        variables = {key: tk.StringVar(value=str(box.get(key, 20 if key.startswith("size") else 0))) for key in BOX_KEYS}
        box_fields = self._box_fields(card, variables)
        def toggle_box(event=None):
            if mode.get() == "Manual box":
                box_fields.pack(fill="x", pady=5)
            else:
                box_fields.pack_forget()
        selector.bind("<<ComboboxSelected>>", toggle_box)
        toggle_box()
        self._text(card, "Coordinates and side lengths are in angstroms. Changing a region makes previous docking results outdated.", "CardMuted.TLabel")
        def apply():
            target["site_mode"] = SITE_LABELS[mode.get()]
            target["box"] = self._box(variables) if target["site_mode"] == "manual" else None
            self.active.save()
            self.run_worker("Updating the target binding regions…", lambda: self.engine.define_binding_sites(self.active, progress=self.progress))
        self._button(card, "Apply binding region", apply)

    def open_view(self, stage):
        from core import viewer
        target = self._selected_target()
        project = self.engine.target_project(self.active, target["id"])
        try:
            molecule = self.drug_choice.get().split(" · ")[-1] if stage == "06_dock" else ""
            data = viewer.payload(project, stage, molecule)
            if not data.get("protein") and not data.get("ligand"):
                raise ValueError("Run preparation or docking before opening this view.")
            script = (self.resource_root / "assets/vendor/3Dmol-min.js").resolve()
            document = viewer.document(data).replace('src="/assets/vendor/3Dmol-min.js"', 'src="' + script.as_uri() + '"')
            destination = project.folder / "exports" / ("view-" + stage + ".html")
            destination.parent.mkdir(exist_ok=True)
            destination.write_text(document, encoding="utf-8")
        finally:
            project.close()
        webbrowser.open(destination.resolve().as_uri())
        self.status.set("Local 3D opened in your browser. The data remains on this computer.")

    def render_dock(self):
        self._summary()
        card = self._card("Run every drug against every target")
        self._target_table(self.active.data["targets"], card)
        row = ttk.Frame(card, style="Card.TFrame")
        row.pack(fill="x")
        self._entry(row, "Search effort", self.effort, 20)
        self._entry(row, "Comparison molecules per target", self.comparisons, 30)
        self._text(card, "Start with a search effort of 8. Higher effort takes longer. Comparison molecules, reference-pose checks and size checks are recorded with each target.", "CardMuted.TLabel")
        def dock():
            effort, number = int(self.effort.get()), int(self.comparisons.get())
            if not 1 <= effort <= 64 or not 1 <= number <= 20:
                raise ValueError("Use a whole-number search effort from 1 to 64 and comparison count from 1 to 20.")
            def done(result):
                self.current, self.section, self.workflow = 4, "Home", True
                self.active.data["ui_step"] = 4
                self.active.save()
            self.run_worker("Running the docking study…", lambda: self.engine.run_docking(self.active, exhaustiveness=effort, n_decoys=number, progress=self.progress), done)
        self._button(card, "Run docking study", dock, True, self._ready_to_dock())
        if not self._ready_to_dock():
            self._text(card, "Prepare the structures and apply a valid binding region in step 3 before docking.", "CardMuted.TLabel")
        self._text(self.body, "Docking produces modelled poses and scores. Compare drugs within each target's setup. These outputs do not establish clinical benefit, potency or safety.", "Muted.TLabel")

    def render_results(self):
        self._summary()
        self.result_rows = self.engine.results_rows(self.active)
        card = self._card("Results for every target–drug combination")
        tree = self._tree(card, ("Target", "Drug", "Status", "Score · kcal/mol", "Scientific status"), [(row["target"], row["drug"], status_label(row.get("status", "")), row.get("score_kcal_mol", ""), status_label(row.get("scientific_status", "not checked"))) for row in self.result_rows], min(12, max(3, len(self.result_rows))), (160, 200, 100, 120, 150))
        details = self._readonly_text(card, "Select a result to read its checks and source links. Failed or incomplete rows remain in the table and exports.", height=8)
        def select(event):
            if not tree.selection():
                return
            row = self.result_rows[int(tree.selection()[0])]
            keys = [("target", "Target"), ("drug", "Drug"), ("status", "Calculation status"), ("scientific_status", "Scientific status"), ("baseline_check", "Baseline check"), ("redocking_check", "Reference-pose check"), ("size_check", "Size check"), ("redock_rmsd_angstrom", "Reference RMSD · angstroms"), ("check_scope", "Check scope"), ("error", "Recorded error"), ("drug_source", "Drug source"), ("drug_source_url", "Drug record"), ("protein_entry_url", "Protein record")]
            text = "\n".join("{}: {}".format(label, row.get(key, "")) for key, label in keys if row.get(key) not in {None, ""})
            details.configure(state="normal")
            details.delete("1.0", "end")
            details.insert("1.0", text)
            details.configure(state="disabled")
        tree.bind("<<TreeviewSelect>>", select)
        actions = ttk.Frame(card, style="Card.TFrame")
        actions.pack(fill="x", pady=8)
        self._button(actions, "Save results CSV", self.save_csv, True)
        self._button(actions, "Save complete study ZIP", self.save_package)
        self._button(actions, "Save methods and checks", self.save_methods)
        self._inspection("06_dock")

    def _readonly_text(self, parent, text, height=10):
        frame = ttk.Frame(parent, style="Card.TFrame")
        frame.pack(fill="both", expand=True, pady=(4, 9))
        box = tk.Text(frame, height=height, wrap="word", font=("Segoe UI", 10), background="white", foreground=INK, relief="solid", borderwidth=1, padx=10, pady=9)
        scrollbar = ttk.Scrollbar(frame, orient="vertical", command=box.yview)
        box.configure(yscrollcommand=scrollbar.set)
        box.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        box.insert("1.0", text)
        box.configure(state="disabled")
        return box

    def save_csv(self):
        filename = filedialog.asksaveasfilename(parent=self.root, title="Save results CSV", defaultextension=".csv", initialfile="docking-results.csv", filetypes=[("CSV table", "*.csv")])
        if filename:
            Path(filename).write_bytes(self.engine.export_csv(self.active))
            self.status.set("Results CSV saved.")

    def save_methods(self):
        filename = filedialog.asksaveasfilename(parent=self.root, title="Save methods and checks", defaultextension=".txt", initialfile="methods-and-checks.txt", filetypes=[("Text document", "*.txt")])
        if filename:
            Path(filename).write_text(self.engine.methods_text(self.active), encoding="utf-8")
            self.status.set("Methods and checks saved.")

    def save_package(self):
        filename = filedialog.asksaveasfilename(parent=self.root, title="Save complete study", defaultextension=".zip", initialfile="docking-study.zip", filetypes=[("Study ZIP archive", "*.zip")])
        if filename:
            def export():
                package = self.engine.export_package(self.active)
                if Path(filename).resolve() != package.resolve():
                    shutil.copy2(package, filename)
                return filename
            self.run_worker("Exporting structures, results, methods and references…", export, lambda result: self.status.set("Complete study saved: " + str(result)))

    def render_reference(self):
        card = self._card("Protein and drug downloads")
        self._text(card, "Download structures from the original scientific databases. Review each protein's chain and binding region before docking.", "CardMuted.TLabel")
        for item in self.references.DOWNLOAD_LINKS:
            title = item.get("title", item.get("name", "Source"))
            self._button(card, title, lambda url=item["url"]: webbrowser.open(url), side="top")
            if item.get("description"):
                self._text(card, item["description"], "CardMuted.TLabel")
        card = self._card("Scientific references")
        references = self.references.references()
        tree = self._tree(card, ("Reference", "Citation"), [(row["title"], row["citation"]) for row in references], 12, (380, 420))
        def open_reference():
            if tree.selection():
                webbrowser.open(references[int(tree.selection()[0])]["url"])
        actions = ttk.Frame(card, style="Card.TFrame")
        actions.pack(fill="x")
        self._button(actions, "Open selected reference", open_reference)
        self._button(actions, "Save references and methods", self.save_references, True)
        self._text(card, "The catalogue covers available tools and methods. Project exports identify the tools, data and settings actually used in each run.", "CardMuted.TLabel")

    def save_references(self):
        filename = filedialog.asksaveasfilename(parent=self.root, title="Save references and methods", defaultextension=".txt", initialfile="references-and-methods.txt", filetypes=[("Text document", "*.txt")])
        if filename:
            projects = [self.active.folder / target["project"] for target in self.active.data["targets"]] if self.active else None
            Path(filename).write_text(self.references.render_references(projects), encoding="utf-8")
            self.status.set("References and methods saved.")

    def render_about(self):
        card = self._card(TITLE)
        if self.logo:
            ttk.Label(card, image=self.logo, style="Card.TLabel").pack(anchor="w", pady=(0, 18))
        self._text(card, "Step-by-step molecular docking for medical education. Native legacy beta edition.")
        self._text(card, "Developed by", "Heading.TLabel")
        for developer in DEVELOPERS:
            self._text(card, developer)
        self._button(card, "Akshaisarathi Healthcare", lambda: webbrowser.open("https://akshaisarathi.com/"))
        note = self._card("Scientific use")
        self._text(note, "Docking produces hypotheses to investigate. A score or pose does not establish a medicine's potency, benefit or safety. Read the checks and source records before interpreting results.")
        self._text(note, "Full source links, scientific references, software versions and preparation records accompany study exports. Optional browser-based 3D does not affect preparation, docking or exported results.", "CardMuted.TLabel")


def launch(project_root=None):
    root = tk.Tk()
    app = DockingStudio(root, project_root=project_root)
    root.mainloop()
    return app

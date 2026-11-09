"""Practical, question-free single and batch docking workspace."""
from __future__ import annotations

import asyncio
import io
import json
from pathlib import Path
from uuid import uuid4

from nicegui import ui

from core import viewer, research
from core.education_ui import DEVELOPERS, apply_theme
from core.settings import projects_root


STEPS = [
    ("Protein targets", "Choose structures for your study", "biotech"),
    ("Drug molecules", "Add the molecules to compare", "science"),
    ("Prepare and inspect", "Prepare structures and binding regions", "view_in_ar"),
    ("Run docking", "Dock every drug against every target", "play_arrow"),
    ("Results and export", "Review checks and save the study", "table_chart"),
]


def register_study() -> None:
    @ui.page("/study")
    def study_page(section: str = "home"):
        from core import study as engine
        from core import references

        client = ui.context.client
        active = None
        current = 0
        busy = False
        message = ""
        selected_target = ""
        selected_drug = ""
        selected_export_target = ""
        selected_export_drug = ""
        selected_export_system = ""
        selected_export_duration = 1.0
        section = {name.lower(): name for name in ("Home", "Projects", "Reference", "About")}.get(section.lower(), "Home")
        workflow = False
        draft_started = False
        library_rows = []
        library_loading = False
        library_error = ""
        draft = {"name": "Docking study", "targets": [], "drugs": [],
                 "smiles": "", "lookup": ""}
        target_draft = {"pdb": "", "name": "", "chains": "", "mode": "ligand",
                        **{f"{prefix}_{axis}": 22 if prefix == "size" else None
                           for prefix in ("center", "size") for axis in "xyz"}}
        apply_theme()
        ui.add_head_html('''<style>
        .study-logo{width:165px;max-width:25vw}.study-shell{max-width:1200px;margin:0 auto}
        .study-card{background:white;border:1px solid #dce7ec;box-shadow:none;padding:23px;
          border-radius:16px}.study-title{font-size:29px;font-weight:700;line-height:1.25;
          letter-spacing:-.5px}.study-instruction{font-size:16px;color:#5d7282;line-height:1.7}
        .study-status{font-size:12px;padding:5px 9px;border-radius:7px;background:#edf5f6}
        .study-step{justify-content:flex-start;min-height:58px;text-align:left;width:100%}
        .study-step .q-btn__content{justify-content:flex-start;flex-wrap:nowrap;gap:12px}
        .study-counter{background:#e6f3f2;border:1px solid #cde6e3;border-radius:13px;
          padding:15px 20px;min-width:120px}.study-upload{width:100%;max-width:none;
          box-shadow:none;border:1px dashed #a9c8ce;border-radius:12px}
        @media(max-width:1000px){.study-logo{width:145px}.study-title{font-size:25px}}
        </style>''')

        def opened_projects():
            if active is None or busy:
                return []
            return [active.folder / target["project"] for target in active.data.get("targets", [])]

        with ui.header().classes("edu-header items-center justify-between gap-3"):
            with ui.row().classes("items-center gap-4 no-wrap"):
                steps_toggle = ui.button(icon="menu", on_click=lambda:drawer.toggle()).props("flat round").classes("lt-md").tooltip("Show project steps")
                steps_toggle.set_visibility(False)
                ui.image("/assets/logo.png").classes("study-logo").props('fit=contain alt="Akshaisarathi Healthcare logo"')
                with ui.column().classes("gap-0 gt-sm"):
                    ui.label("In-Silico Docking Studio").classes("text-lg font-bold")
                    ui.label("MEDICAL EDUCATION").classes("edu-kicker").style("font-size:10px;letter-spacing:1px")
            with ui.row().classes("items-center gap-1 flex-wrap").props('data-testid="primary-navigation"'):
                main_buttons = {}
                for title, icon in [("Home", "home"), ("Projects", "folder"), ("Reference", "menu_book"), ("About", "info_outline")]:
                    main_buttons[title] = ui.button(title, on_click=lambda name=title: show_section(name)).props(f"flat icon={icon} no-caps").classes("edu-small-button")

        with ui.left_drawer(value=False).props("width=250 breakpoint=850").style("background:#fff;border-right:1px solid #dce7ec") as drawer:
            nav = ui.column().classes("w-full gap-2 p-4")
        with ui.column().classes("study-shell w-full p-5 md:p-8 gap-5"):
            with ui.row().classes("edu-info w-full items-center gap-3") as working:
                ui.spinner(size="sm", color="primary")
                progress_label = ui.label().classes("text-sm font-medium")
            working.set_visibility(False)
            body = ui.column().classes("w-full gap-5")

        def refresh_progress():
            working.set_visibility(busy)
            if busy:
                progress_label.set_text(message or "Processing the study…")

        ui.timer(0.25, refresh_progress)

        def progress(text):
            nonlocal message
            message = str(text)

        def notify_error(error):
            with client:
                ui.notify(str(error), type="negative", timeout=10000)

        def notify_success(text, timeout=5000):
            with client:
                ui.notify(text, type="positive", timeout=timeout)

        async def load_library():
            nonlocal library_rows, library_loading, library_error
            from core import project_library
            if library_loading:
                return
            library_loading = True
            library_error = ""
            render()
            try:
                library_rows = await asyncio.to_thread(project_library.list_projects, projects_root())
            except Exception as error:
                library_error = str(error)
            finally:
                library_loading = False
                render()

        async def show_section(name):
            nonlocal section, workflow
            section = name
            workflow = False
            render()
            if name == "Projects":
                await load_library()

        def continue_project():
            nonlocal section, workflow
            section, workflow = "Home", True
            render()

        def navigate(index):
            nonlocal current, section, workflow
            if busy:
                return
            if index > 1 and active is None:
                return
            current = index
            section, workflow = "Home", True
            if active is not None:
                active.data["ui_step"] = current
                active.save()
            render()

        async def operation(label, function, next_step=None):
            nonlocal busy, message, current
            if busy:
                return
            busy = True
            message = label
            render()
            try:
                result = await asyncio.to_thread(function)
                if next_step is not None:
                    current = next_step
                    if active is not None:
                        active.data["ui_step"] = current
                        active.save()
                return result
            except Exception as error:
                notify_error(error)
            finally:
                busy = False
                with client:
                    render()

        def new_study(copy_inputs=False):
            nonlocal active, current, selected_target, selected_drug, section, workflow, draft_started
            if busy:
                return
            if copy_inputs and active is not None:
                draft["name"] = active.data.get("name", "Docking study") + " copy"
                draft["targets"] = [{key: row.get(key) for key in
                                      ("name", "pdb_id", "file", "chains", "site_mode", "box", "reference_project")}
                                     for row in active.data.get("targets", [])]
                for target in draft["targets"]:
                    if target.get("file") and not Path(target["file"]).is_absolute():
                        target["file"] = str(active.folder / target["file"])
                draft["drugs"] = [dict(row) for row in active.data.get("drugs", [])]
            else:
                draft.update(name="Docking study", targets=[], drugs=[], smiles="", lookup="")
            target_draft.update(pdb="", name="", chains="", mode="ligand",
                                **{f"{prefix}_{axis}": 22 if prefix == "size" else None
                                   for prefix in ("center", "size") for axis in "xyz"})
            active = None
            current = 0
            section, workflow, draft_started = "Home", True, True
            selected_target = selected_drug = ""
            render()

        def resume(folder):
            nonlocal active, current, selected_target, selected_drug, section, workflow
            if busy:
                return
            try:
                active = engine.open_study(folder)
                current = min(max(int(active.data.get("ui_step", 2)), 0), len(STEPS)-1)
                section, workflow = "Home", True
                selected_target = selected_drug = ""
                render()
            except Exception as error:
                notify_error(error)

        def create():
            nonlocal active, current
            try:
                if draft["smiles"].strip():
                    add_smiles(False)
                if not draft["targets"]:
                    raise ValueError("Add at least one protein target.")
                if not draft["drugs"]:
                    raise ValueError("Add at least one drug molecule.")
                targets = [engine.TargetSpec(**{key: row.get(key) for key in
                           ("name", "pdb_id", "file", "chains", "site_mode", "box", "reference_project")})
                           for row in draft["targets"]]
                drugs = [engine.DrugSpec(**{key: row.get(key, "") for key in
                         ("name", "smiles", "cid", "source_url", "source", "retrieved")})
                         for row in draft["drugs"]]
                active = engine.create_study(projects_root(), draft["name"], targets, drugs)
                active.data["drug_sources"] = {row["name"]: {k: row[k] for k in
                    ("cid", "source_url", "source", "retrieved") if row.get(k)} for row in draft["drugs"]}
                active.data["ui_step"] = current = 2
                active.save()
                render()
                notify_success("Study saved. Prepare the structures in step 3.")
            except Exception as error:
                notify_error(error)

        def example():
            nonlocal active, current, section, workflow
            if busy:
                return
            try:
                active = engine.create_example_study(projects_root())
                active.data["ui_step"] = current = 2
                section, workflow = "Home", True
                active.save()
                render()
                notify_success("EGFR example saved. Its reference preparation is labelled; docking runs are generated when you run step 4.", timeout=8000)
            except Exception as error:
                notify_error(error)

        def add_drug(row):
            from rdkit import Chem
            name = str(row.get("name", "")).strip()
            smiles = str(row.get("smiles", "")).strip()
            if not name or Chem.MolFromSmiles(smiles) is None:
                raise ValueError("Each drug needs a name and a valid SMILES structure.")
            if any(existing["name"].casefold() == name.casefold() for existing in draft["drugs"]):
                raise ValueError(f"A drug named {name} is already listed. Use a distinct name for a different molecule.")
            draft["drugs"].append({**row, "name": name, "smiles": smiles})

        def add_smiles(refresh=True):
            lines = [line.strip() for line in draft["smiles"].splitlines() if line.strip()]
            proposed = []
            for index, line in enumerate(lines, 1):
                parts = line.split(None, 1)
                proposed.append({"smiles": parts[0], "name": parts[1].strip() if len(parts)>1 else f"Drug {len(draft['drugs'])+index}",
                                 "source": "User-supplied SMILES"})
            before = list(draft["drugs"])
            try:
                for row in proposed:
                    add_drug(row)
            except Exception:
                draft["drugs"] = before
                raise
            draft["smiles"] = ""
            if refresh:
                render()

        def add_smiles_clicked():
            try:
                add_smiles()
            except Exception as error:
                notify_error(error)

        async def lookup_drugs():
            queries = [line.strip() for line in draft["lookup"].splitlines() if line.strip()]
            if not queries:
                notify_error("Enter a drug name or PubChem CID before adding molecules.")
                return
            found, failed = [], []
            def fetch():
                for index, query in enumerate(queries, 1):
                    progress(f"PubChem: {query} · {index} of {len(queries)}")
                    try:
                        found.append(references.resolve_drug(query, projects_root()/".drug-cache"))
                    except Exception as error:
                        failed.append((query, str(error)))
            await operation("Fetching drug structures from PubChem…", fetch)
            for row in found:
                try:
                    add_drug(row)
                except ValueError as error:
                    failed.append((row.get("name", "Drug"), str(error)))
            draft["lookup"] = "\n".join(row[0] for row in failed)
            render()
            if failed:
                notify_error("Some drugs were not added: " + "; ".join(f"{name}: {error}" for name, error in failed))

        async def upload_drugs(event):
            try:
                from rdkit import Chem
                filename = Path(event.file.name).name
                if Path(filename).suffix.lower() != ".sdf":
                    raise ValueError("Upload an SDF file for drug structures.")
                data = await event.file.read()
                rows = []
                for index, molecule in enumerate(Chem.ForwardSDMolSupplier(io.BytesIO(data)), 1):
                    if index > 500:
                        raise ValueError("An SDF import supports up to 500 molecules per file.")
                    if molecule is None:
                        raise ValueError(f"Molecule {index} in the SDF file could not be read.")
                    rows.append({"name": molecule.GetProp("_Name").strip() if molecule.HasProp("_Name") and molecule.GetProp("_Name").strip() else f"{Path(filename).stem} {index}",
                                 "smiles": Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True),
                                 "source": f"Uploaded SDF: {filename}"})
                before = list(draft["drugs"])
                try:
                    for row in rows:
                        add_drug(row)
                except Exception:
                    draft["drugs"] = before
                    raise
                render()
                notify_success(f"Added {len(rows)} drug structures.")
            except Exception as error:
                notify_error(error)

        def remove_row(key, index):
            draft[key].pop(index)
            render()

        def summary():
            targets = active.data.get("targets", []) if active is not None else draft["targets"]
            drugs = active.data.get("drugs", []) if active is not None else draft["drugs"]
            with ui.row().classes("gap-3 flex-wrap"):
                for count, label in [(len(targets), "protein targets"), (len(drugs), "drug molecules"), (len(targets)*len(drugs), "target–drug pairs")]:
                    with ui.column().classes("study-counter gap-0"):
                        ui.label(str(count)).classes("text-2xl font-bold")
                        ui.label(label).classes("text-xs edu-muted")

        def render_control(project, sid):
            row = project.get_step(sid)
            if not row:
                return
            if row.get("status") != "complete":
                ui.label("Needs updating. Run this stage again before using its earlier checks.").classes("edu-info text-sm w-full")
                return
            control = json.loads(row.get("control_json") or "{}")
            if not control:
                return
            if sid in {"03_receptor", "04_site", "06_dock"}:
                from core.workshop import teaching_check
                control = teaching_check(sid, control, recorded=json.loads(row.get("values_json") or "{}"))
            verdict = control.get("verdict", "not_applicable")
            colours = {"pass":"#26754d", "warn":"#a86a10", "fail":"#b63d48", "not_applicable":"#62798a"}
            labels = {"pass":"Checked", "warn":"Needs review", "fail":"Criterion not met", "not_applicable":"Not checked"}
            def display_text(value):
                text = str(value)
                for drug in sorted(active.data.get("drugs", []), key=lambda item:len(item.get("internal_name", "")), reverse=True):
                    internal = drug.get("internal_name", "")
                    if internal:
                        text = text.replace(internal, drug["name"])
                return text
            with ui.card().classes("w-full gap-2 p-4").style(f"box-shadow:none;border:1px solid #dce7ec;border-left:4px solid {colours.get(verdict,'#62798a')}"):
                heading = {"03_receptor":"Preparation changes", "04_site":"Site geometry", "06_dock":"Docking review"}.get(sid,"Scientific check")
                ui.label(("Review recorded checks" if sid == "06_dock" else labels.get(verdict, verdict)) + " · " + heading).classes("font-semibold")
                ui.label(display_text(control.get("plain_verdict", ""))).classes("text-sm").style("line-height:1.65")
                if control.get("what_to_do"):
                    ui.label("Next action: " + display_text(control["what_to_do"])).classes("text-sm font-medium")
                if sid == "06_dock":
                    from core.workshop import docking_checks
                    checks = docking_checks(json.loads(row.get("values_json") or "{}"), control)
                    for name, check in checks.items():
                        status = check.get("verdict", "not_applicable")
                        status = getattr(status,"value",status)
                        label = {"baseline":"Score comparison", "re-docking":"Known pose recovery", "size":"Effect of molecular size"}.get(name, name.capitalize())
                        ui.label(f"{label} · {labels.get(status,status)}: {display_text(check.get('plain_verdict',''))}").classes("text-sm")
                        if name == "baseline":
                            detail = check.get("detail", {})
                            if "best decoy" in detail:
                                ui.label(f"Best comparison score: {detail['best decoy']} kcal/mol. Comparison molecules run: {detail.get('decoys run', 'Not recorded')}.").classes("text-sm")
                    ui.label("A score comparison does not overturn experimental evidence for a known drug and target. Comparison molecules are not proven nonbinders.").classes("text-sm")

        def inspect(stage):
            nonlocal selected_target, selected_drug
            targets = active.data.get("targets", [])
            if not targets:
                return
            ids = {row["id"]: row["name"] for row in targets}
            if selected_target not in ids:
                selected_target = next(iter(ids))
            def choose_target(event):
                nonlocal selected_target, selected_drug
                selected_target = event.value
                selected_drug = ""
                render()
            ui.select(ids, value=selected_target, label="Protein target to inspect", on_change=choose_target).props("outlined").classes("w-full")
            try:
                project = engine.target_project(active, selected_target)
                try:
                    if stage == "06_dock":
                        results = viewer.values(project, stage).get("results", [])
                        names = {row.get("internal_name", row["name"]): row["name"] for row in active.data.get("drugs", [])}
                        options = {row["name"]: names.get(row["name"],row["name"]) for row in results if row.get("file")}
                        if options:
                            if selected_drug not in options:
                                selected_drug = next(iter(options))
                            def choose_drug(event):
                                nonlocal selected_drug
                                selected_drug = event.value
                                render()
                            ui.select(options, value=selected_drug, label="Docked drug to view", on_change=choose_drug).props("outlined").classes("w-full")
                    if stage == "04_site":
                        with ui.row().classes("w-full items-start gap-4 flex-wrap"):
                            for original, label in [(True, "Protein before preparation"), (False, "Protein after preparation")]:
                                with ui.column().classes("gap-2").style("flex:1;min-width:300px"):
                                    ui.label(label).classes("font-semibold")
                                    viewer.render(project, "03_receptor", original=original)
                        ui.label("Use Show atoms to inspect the added atoms. Small red spheres show water molecules. Other retained components are shown as sticks.").classes("text-sm edu-muted")
                        review = research.preparation_review(project)
                        if review:
                            ui.label("What preparation changed").classes("text-lg font-semibold")
                            ui.table(columns=[{"name":"change","label":"Change","field":"change","align":"left"},
                                              {"name":"value","label":"Number","field":"value"},
                                              {"name":"meaning","label":"Simple meaning","field":"meaning","align":"left"}],
                                     rows=review, row_key="change").classes("w-full").props("flat wrap-cells")
                            receptor = viewer.values(project, "03_receptor")
                            for key, label in [("kept", "Components kept"), ("removed", "Components removed")]:
                                ui.label(label + ": " + (", ".join(x.get("label", "") + " " + x.get("what", "") for x in receptor.get(key, [])) or "None recorded")).classes("text-sm")
                        for sid, key, label in [("02_structure", "local_file", "Download original protein"), ("03_receptor", "prepared_file", "Download prepared protein")]:
                            data = viewer.values(project, sid)
                            filename = (data.get("chosen", {}) if sid == "02_structure" else data).get(key, "")
                            if viewer.local_text(project, filename):
                                ui.button(label, on_click=lambda path=filename: ui.download(path)).props("outline no-caps icon=download")
                        ui.label("Prepared drug structure").classes("text-lg font-semibold")
                        prepared = viewer.values(project, "05_ligands").get("prepared", [])
                        if prepared:
                            names = {d["internal_name"]: d["name"] for d in active.data["drugs"]}
                            options = {d["name"]: names.get(d["name"], d["name"]) for d in prepared}
                            if selected_drug not in options:
                                selected_drug = next(iter(options))
                            def choose_prepared(event):
                                nonlocal selected_drug
                                selected_drug = event.value
                                render()
                            ui.select(options, value=selected_drug, label="Prepared drug to inspect", on_change=choose_prepared).props("outlined").classes("w-full")
                            molecule = next(d for d in prepared if d["name"] == selected_drug)
                            viewer.render(project, "05_ligands", selected_drug)
                            for line in research.drug_review(molecule).splitlines():
                                ui.label(line).classes("text-sm")
                            ui.label("Prepared SMILES: " + molecule["smiles"]).classes("text-xs edu-muted").style("overflow-wrap:anywhere")
                            if viewer.local_text(project, molecule.get("file", "")):
                                ui.button("Download prepared drug SDF", on_click=lambda path=molecule["file"]: ui.download(path)).props("outline no-caps icon=download")
                        else:
                            ui.label("Prepare the drug to see its 3D structure and recorded changes.").classes("text-sm")
                        ui.label("Protein pocket and docking box").classes("text-lg font-semibold")
                        viewer.render(project, stage)
                    else:
                        viewer.render(project, stage, selected_drug)
                    for sid in (["03_receptor", "04_site"] if stage == "04_site" else ["06_dock"]):
                        render_control(project, sid)
                finally:
                    project.close()
            except Exception as error:
                ui.label(f"The view is not ready: {error}").classes("text-sm edu-muted")

        def target_rows():
            rows = []
            for target in active.data.get("targets", []):
                detail = target.get("error") or target.get("site_error") or target.get("site_summary") or target.get("preparation_summary") or "Ready to prepare"
                if target.get("source_label"):
                    detail = target["source_label"] + " " + detail
                status = {"pending":"Not prepared", "prepared":"Ready", "docked":"Docked", "partial":"Partly complete", "failed":"Failed", "interrupted":"Interrupted", "running":"Running"}.get(target.get("status"),target.get("status","Not prepared"))
                rows.append({"name": target["name"], "status": status, "summary":detail})
            ui.table(columns=[{"name":"name","label":"Protein target","field":"name","align":"left"},
                              {"name":"status","label":"Status","field":"status","align":"left"},
                              {"name":"summary","label":"Preparation and site summary","field":"summary","align":"left"}],
                     rows=rows,row_key="name").classes("w-full").props("wrap-cells flat")

        async def apply_site(target_id, mode_field, fields):
            try:
                mode = mode_field.value
                box = None
                if mode == "manual":
                    box = {key: float(field.value) for key, field in fields.items()}
                    if any(box[key] <= 0 for key in ("size_x", "size_y", "size_z")):
                        raise ValueError("Search-box dimensions must be positive.")
                target = next(row for row in active.data["targets"] if row["id"] == target_id)
                target["site_mode"] = mode
                target["box"] = box
                active.save()
                await operation("Updating the binding regions…", lambda: engine.define_binding_sites(active, progress=progress))
            except Exception as error:
                notify_error(error)

        def site_editors():
            for target in active.data.get("targets", []):
                with ui.expansion(target["name"] + " · Binding-region settings", icon="my_location").classes("w-full study-card"):
                    ui.label("Set the search region for this protein. Review its location in the 3D view before docking.").classes("text-sm edu-muted")
                    if target.get("site_error"):
                        ui.label("Binding region needs setup: " + target["site_error"]).classes("edu-info text-sm w-full")
                    mode = ui.select({"ligand":"Use the bound reference molecule", "detect":"Detect a candidate pocket", "manual":"Enter search-box coordinates"},
                                     value=target.get("site_mode") or "ligand", label="Binding-region method").props("outlined").classes("w-full")
                    box = target.get("box") or {}
                    fields = {}
                    with ui.row().classes("w-full gap-3 flex-wrap") as coordinates:
                        for prefix, label in [("center", "Centre"), ("size", "Size")]:
                            for axis in "xyz":
                                key = prefix + "_" + axis
                                fields[key] = ui.number(f"{label} {axis.upper()} (Å)", value=box.get(key), step=0.5).props("outlined dense").style("width:130px")
                    coordinates.bind_visibility_from(mode, "value", lambda value: value == "manual")
                    ui.button("Apply binding-region settings", on_click=lambda tid=target["id"], m=mode, f=fields: apply_site(tid,m,f)).props("outline no-caps icon=check")

        def render_targets():
            if active is not None:
                ui.label("This study's target inputs are saved. Copy its inputs into a new study to change the target or drug set while keeping earlier results.").classes("study-instruction")
                target_rows()
                ui.button("Copy inputs into a new study", on_click=lambda: new_study(True)).props("outline no-caps icon=content_copy")
                ui.button("Continue to drug molecules", on_click=lambda: navigate(1)).props("unelevated color=primary no-caps icon-right=arrow_forward")
                return
            with ui.card().classes("study-card w-full gap-4"):
                ui.input("Study name").bind_value(draft,"name").props("outlined").classes("w-full")
                ui.label("Add a protein target").classes("text-lg font-semibold")
                ui.label("Use a PDB entry ID for an online structure download, or upload a local PDB/mmCIF file. Add each target separately so its binding region can be reviewed.").classes("study-instruction")
                with ui.row().classes("w-full gap-3 flex-wrap"):
                    pdb = ui.input("PDB entry ID", placeholder="1M17").bind_value(target_draft,"pdb").props("outlined").style("flex:1;min-width:160px")
                    name = ui.input("Target label", placeholder="EGFR").bind_value(target_draft,"name").props("outlined").style("flex:2;min-width:180px")
                    chains = ui.input("Protein chains (optional)",placeholder="A or A,B").bind_value(target_draft,"chains").props("outlined").style("flex:1;min-width:180px")
                mode = ui.select({"ligand":"Use the bound reference molecule", "detect":"Detect a candidate pocket", "manual":"Enter search-box coordinates"},
                                 label="Binding-region method").bind_value(target_draft,"mode").props("outlined").classes("w-full")
                manual_fields = {}
                with ui.row().classes("w-full gap-3 flex-wrap") as manual_coordinates:
                    for prefix, label in [("center", "Centre"), ("size", "Size")]:
                        for axis in "xyz":
                            key = prefix + "_" + axis
                            manual_fields[key] = ui.number(f"{label} {axis.upper()} (Å)", step=0.5).bind_value(target_draft,key).props("outlined dense").style("width:130px")
                manual_coordinates.bind_visibility_from(mode,"value",lambda value:value=="manual")
                ui.label("A bound reference molecule locates a measured site. Pocket detection proposes a site for inspection; it does not establish a binding site.").classes("text-sm edu-muted")
                def append_target(local_file=""):
                    code = (pdb.value or "").strip().upper()
                    label = (name.value or "").strip() or code or Path(local_file).stem
                    if not code and not local_file:
                        raise ValueError("Enter a PDB ID or upload a local protein structure.")
                    if any(row["name"].casefold() == label.casefold() for row in draft["targets"]):
                        raise ValueError("Use a distinct target label for each structure.")
                    box = None
                    if mode.value == "manual":
                        try:
                            box = {key:float(field.value) for key,field in manual_fields.items()}
                        except (TypeError,ValueError):
                            raise ValueError("Enter all six search-box coordinates and side lengths.") from None
                    draft["targets"].append({"name":label,"pdb_id":code if not local_file else "", "file":local_file,
                         "chains":[s.strip() for s in (chains.value or "").replace(";",",").split(",") if s.strip()] or None,
                         "site_mode":mode.value,"box":box,"reference_project":""})
                    target_draft.update(pdb="", name="", chains="")
                def add_target():
                    try:
                        append_target()
                        render()
                    except Exception as error:
                        notify_error(error)
                ui.button("Add protein target", on_click=add_target).props("unelevated color=primary no-caps icon=add")
                async def upload_target(event):
                    try:
                        filename = Path(event.file.name).name
                        suffix = Path(filename).suffix.lower()
                        if suffix not in {".pdb", ".cif", ".mmcif"}:
                            raise ValueError("Upload a PDB or mmCIF protein structure.")
                        data = await event.file.read()
                        text = data.decode("utf-8", errors="replace")
                        if not (any(line.startswith("ATOM  ") for line in text.splitlines()) or "_atom_site." in text):
                            raise ValueError("The uploaded file does not contain protein atom records.")
                        folder = projects_root()/".study-uploads"/uuid4().hex
                        folder.mkdir(parents=True, exist_ok=True)
                        path = folder/filename
                        path.write_bytes(data)
                        append_target(str(path))
                        render()
                    except Exception as error:
                        notify_error(error)
                ui.upload(label="Upload local protein structures", on_upload=upload_target, auto_upload=True,
                          max_file_size=25_000_000).props('accept=".pdb,.cif,.mmcif"').classes("study-upload")
                with ui.row().classes("gap-5 flex-wrap"):
                    ui.link("Find and download proteins · RCSB PDB", "https://www.rcsb.org/", new_tab=True)
                    ui.link("Protein information · UniProt", "https://www.uniprot.org/", new_tab=True)
            if draft["targets"]:
                with ui.card().classes("study-card w-full gap-3"):
                    ui.label("Targets in this study").classes("text-lg font-semibold")
                    for index, target in enumerate(draft["targets"]):
                        with ui.row().classes("w-full items-center gap-3"):
                            ui.icon("biotech", color="primary")
                            with ui.column().classes("gap-0"):
                                ui.label(target["name"]).classes("font-semibold")
                                ui.label(target.get("pdb_id") or Path(target.get("file","")).name).classes("text-xs edu-muted")
                                if target.get("pdb_id"):
                                    for link in references.structure_links(target["pdb_id"]):
                                        ui.link(link["label"],link["url"],new_tab=True).classes("text-xs")
                            ui.space()
                            ui.button(icon="close", on_click=lambda i=index: remove_row("targets",i)).props("flat round").tooltip("Remove target")
                    ui.button("Continue to drug molecules", on_click=lambda: navigate(1)).props("unelevated color=primary no-caps icon-right=arrow_forward")
        def render_drugs():
            if active is not None:
                ui.label("Saved molecules").classes("text-lg font-semibold")
                render_drug_list(active.data.get("drugs", []), removable=False)
                ui.button("Continue to preparation", on_click=lambda: navigate(2)).props("unelevated color=primary no-caps icon-right=arrow_forward")
                return
            ui.label("Add drug names from PubChem, paste SMILES, or import an SDF file. Every listed molecule will be docked against every listed target.").classes("study-instruction")
            with ui.card().classes("study-card w-full gap-3"):
                ui.label("Add from PubChem · online").classes("text-lg font-semibold")
                ui.textarea("Drug names or PubChem CIDs — one per line", placeholder="erlotinib\ngefitinib").bind_value(draft,"lookup").props("outlined rows=3").classes("w-full")
                ui.button("Add drugs from PubChem", on_click=lookup_drugs).props("unelevated color=primary no-caps icon=download")
                ui.link("Open PubChem to find and download drug structures", "https://pubchem.ncbi.nlm.nih.gov/", new_tab=True)
            with ui.expansion("Paste molecular structures · SMILES", value=True, icon="content_paste").classes("study-card w-full"):
                ui.label("Paste a SMILES code followed by the molecule name. One molecule per line.").classes("text-sm edu-muted")
                ui.textarea("SMILES and drug name", placeholder="CC(=O)Oc1ccccc1C(=O)O Aspirin\nCCO Ethanol").bind_value(draft,"smiles").props("outlined rows=4").classes("w-full")
                ui.button("Add listed molecules", on_click=add_smiles_clicked).props("outline no-caps icon=add")
            ui.upload(label="Import local drug structures · SDF",on_upload=upload_drugs,auto_upload=True,
                      max_file_size=25_000_000).props('accept=".sdf"').classes("study-upload")
            render_drug_list(draft["drugs"], removable=True)
            summary()
            ui.button("Save study and prepare structures", on_click=create).props("unelevated color=primary no-caps icon-right=arrow_forward").set_enabled(bool(draft["targets"]))

        def render_drug_list(drugs, removable):
            if not drugs:
                return
            with ui.card().classes("study-card w-full gap-3"):
                ui.label(f"{len(drugs)} drug molecules").classes("text-lg font-semibold")
                for index, drug in enumerate(drugs):
                    with ui.row().classes("w-full items-start gap-3 no-wrap"):
                        ui.icon("science", color="primary")
                        with ui.column().classes("gap-1").style("min-width:0;flex:1"):
                            ui.label(drug["name"]).classes("font-semibold")
                            ui.label(drug["smiles"]).classes("text-xs edu-muted").style("word-break:break-all")
                            for link in references.drug_links(drug["name"], drug["smiles"], str(drug.get("cid", ""))):
                                ui.link(link["label"],link["url"],new_tab=True).classes("text-xs")
                        if removable:
                            ui.button(icon="close", on_click=lambda i=index: remove_row("drugs",i)).props("flat round").tooltip("Remove drug")

        def render_prepare():
            summary()
            ui.label("Prepare each protein and drug, then inspect each target's binding region. Preparation choices and scientific checks are recorded with the study.").classes("study-instruction")
            ui.button("Prepare proteins and drugs", on_click=lambda: operation("Preparing protein and drug structures…", lambda: engine.run_preparation(active, progress=progress))).props("unelevated color=primary no-caps icon=build").set_enabled(not busy)
            target_rows()
            with ui.card().classes("study-card w-full gap-4"):
                ui.label("Inspect the protein and search box").classes("text-lg font-semibold")
                inspect("04_site")
            site_editors()
            try:
                project = engine.target_project(active, active.data["targets"][0]["id"])
                try:
                    descriptors = viewer.values(project,"08_admet").get("molecules",[])
                    names = {row.get("internal_name",row["name"]):row["name"] for row in active.data["drugs"]}
                    if descriptors:
                        with ui.expansion("Calculated drug descriptors", icon="medication").classes("study-card w-full"):
                            ui.label("These values describe the prepared chemical structures. They are not measured pharmacokinetics or toxicity predictions.").classes("text-sm edu-muted")
                            values = [{"name":names.get(row["name"],row["name"]),**row.get("properties",{})} for row in descriptors]
                            ui.table(columns=[{"name":"name","label":"Drug","field":"name","align":"left"},
                                              {"name":"mw","label":"Molecular weight","field":"mw"},
                                              {"name":"logp","label":"Estimated logP","field":"logp"},
                                              {"name":"tpsa","label":"Polar surface area","field":"tpsa"},
                                              {"name":"hbd","label":"H-bond donors","field":"hbd"},
                                              {"name":"hba","label":"H-bond acceptors","field":"hba"}],rows=values,row_key="name").classes("w-full").props("flat wrap-cells")
                finally:
                    project.close()
            except (OSError,ValueError):
                pass
            render_swissadme()
            ui.label("Changing a binding region makes earlier docking results outdated. Run docking again after applying a change.").classes("edu-info text-sm w-full")
            ready = any(target.get("status") in {"prepared","docked","partial"} and target.get("site_status") not in {"needs_input","fail"} for target in active.data["targets"])
            ui.button("Continue to docking", on_click=lambda: navigate(3)).props("unelevated color=primary no-caps icon-right=arrow_forward").set_enabled(ready)

        def render_dock():
            summary()
            ui.label("Run the complete target–drug matrix").classes("text-lg font-semibold")
            ui.label("Each target keeps its own prepared structure and search box. Targets are processed in sequence; a failed target is recorded without discarding other targets' results.").classes("study-instruction")
            target_rows()
            with ui.card().classes("study-card w-full gap-3"):
                effort = ui.number("Search effort", value=active.data.get("settings",{}).get("exhaustiveness",8),min=1,max=64,step=1,format="%.0f").props("outlined").classes("w-full")
                ui.label("Start at 8. Higher effort takes longer. Reference-pose, comparison-molecule and size checks remain part of the workflow.").classes("text-sm edu-muted")
                comparisons = ui.number("Comparison molecules per target",value=active.data.get("settings",{}).get("n_decoys",3),min=1,max=20,step=1,format="%.0f").props("outlined").classes("w-full")
                def run():
                    try:
                        value = int(effort.value)
                        number = int(comparisons.value)
                        if value != effort.value or number != comparisons.value or not 1 <= value <= 64 or not 1 <= number <= 20:
                            raise ValueError("Enter a whole-number search effort from 1 to 64 and comparison count from 1 to 20.")
                        return operation("Running molecular docking…", lambda: engine.run_docking(active,exhaustiveness=value,n_decoys=number,progress=progress),next_step=4)
                    except (TypeError,ValueError) as error:
                        notify_error(error)
                ready = any(target.get("status") in {"prepared","docked","partial"} and target.get("site_status") not in {"needs_input","fail"} for target in active.data["targets"])
                ui.button("Run docking study",on_click=run).props("unelevated color=primary no-caps icon=play_arrow size=lg").set_enabled(not busy and ready)
                if not ready:
                    ui.label("Prepare the structures in step 3 before running docking.").classes("text-sm edu-muted")
            ui.label("Docking produces modelled poses and scores. Read all checks before interpreting a ranking; these outputs do not establish potency, clinical benefit or safety.").classes("edu-info text-sm w-full")

        async def download_package():
            package = await operation("Preparing the complete study package…", lambda: engine.export_package(active))
            if package:
                with client:
                    ui.download(str(package))

        def render_swissadme():
            with ui.expansion("Drug behaviour in the body and SwissADME", icon="medication").classes("study-card w-full"):
                ui.label("SwissADME predicts properties related to gut absorption, entry into the brain and drug metabolism. These predictions are not measured pharmacokinetics or dose advice.").classes("text-sm")
                ui.label("Open SwissADME, paste the SMILES below and press Run. You can download its report as CSV from that website.").classes("text-sm")
                ui.link("Open SwissADME", research.SWISSADME_URL, new_tab=True)
                smiles = research.swissadme_input(active)
                ui.code(smiles.decode("utf-8")).classes("w-full").style("max-height:220px;overflow:auto")
                ui.button("Download SMILES for SwissADME", on_click=lambda: ui.download(smiles, filename="SwissADME_input.txt")).props("outline no-caps icon=download")

        def render_code_exports():
            nonlocal selected_export_target, selected_export_drug, selected_export_system
            with ui.card().classes("study-card w-full gap-3"):
                ui.label("Download code for a selected protein and drug").classes("text-lg font-semibold")
                targets = {t["id"]: t["name"] for t in active.data["targets"]}
                drugs = {d["id"]: d["name"] for d in active.data["drugs"]}
                if selected_export_target not in targets:
                    selected_export_target = selected_target if selected_target in targets else next(iter(targets))
                selected_id = next((d["id"] for d in active.data["drugs"] if d["internal_name"] == selected_drug), next(iter(drugs)))
                if selected_export_drug not in drugs:
                    selected_export_drug = selected_id
                if not selected_export_system:
                    selected_export_system = "membrane" if "SV2" in targets[selected_export_target].upper() else "water"
                def select_export(event, field):
                    nonlocal selected_export_target, selected_export_drug, selected_export_system, selected_export_duration
                    if field == "target":
                        selected_export_target = event.value
                        if "SV2" in targets.get(event.value, "").upper():
                            system_field.set_value("membrane")
                    elif field == "drug":
                        selected_export_drug = event.value
                    elif field == "system":
                        selected_export_system = event.value
                    else:
                        selected_export_duration = event.value
                target_field = ui.select(targets, value=selected_export_target, label="Protein for code export", on_change=lambda e: select_export(e, "target")).props("outlined").classes("w-full")
                drug_field = ui.select(drugs, value=selected_export_drug, label="Drug for code export", on_change=lambda e: select_export(e, "drug")).props("outlined").classes("w-full")
                async def export_code(kind):
                    target_id, drug_id = target_field.value, drug_field.value
                    mode, duration = system_field.value, duration_field.value
                    package = await operation("Preparing the selected code package…", lambda: research.export_pair(active, target_id, drug_id, kind, system_type=mode, nanoseconds=duration))
                    if package:
                        with client:
                            ui.download(package, filename="Docking_Colab.zip" if kind == "docking" else "MD_Colab.zip")
                ui.label("Docking code includes the actual prepared inputs, search box, recorded settings and a notebook for Google Colab. The package reruns the selected pair. It keeps the original comparison outcomes in its report.").classes("text-sm")
                with ui.row().classes("gap-3 flex-wrap"):
                    ui.button("Download docking code and inputs", on_click=lambda: export_code("docking")).props("unelevated color=primary no-caps icon=download")
                    ui.button("Download docking notebook", on_click=lambda: ui.download(research.notebook("Dock your selected protein and drug", "run_docking.py"), filename="Run_Docking.ipynb")).props("outline no-caps")
                ui.label("MD simulation code").classes("font-semibold")
                system_field = ui.select({"water":"Soluble protein in water", "membrane":"Membrane protein such as SV2A"}, value=selected_export_system, label="MD system", on_change=lambda e: select_export(e, "system")).props("outlined").classes("w-full")
                duration_field = ui.number("Production simulation length in ns", value=selected_export_duration, min=0.001, max=1000, on_change=lambda e: select_export(e, "duration")).props("outlined").classes("w-full")
                ui.label("The MD package contains the selected docked drug with its chemical bond information. Membrane mode needs a separately prepared lipid system. Code is exported here; the simulation runs in Colab or on your own computer.").classes("text-sm")
                ui.link("Membrane system preparation with CHARMM GUI", "https://www.charmm-gui.org/", new_tab=True)
                with ui.row().classes("gap-3 flex-wrap"):
                    ui.button("Download MD code and inputs", on_click=lambda: export_code("md")).props("unelevated color=primary no-caps icon=download")
                    ui.button("Download MD notebook", on_click=lambda: ui.download(research.notebook("MD for your selected protein and drug", "run_md.py", md=True), filename="Run_MD.ipynb")).props("outline no-caps")
                ui.link("Open Google Colab", "https://colab.research.google.com/", new_tab=True)

        def render_results():
            summary()
            rows = engine.results_rows(active)
            normalized = []
            for index, row in enumerate(rows):
                normalized.append({**row,"pair_id":str(index),
                    "target":row.get("target",row.get("target_name","")),
                    "drug":row.get("drug",row.get("drug_name",row.get("molecule",""))),
                    "score":row.get("score",row.get("score_kcal_mol","")),
                    "status":{"completed":"Completed","failed":"Failed","pending":"Not run","needs_updating":"Needs updating","running":"Running","interrupted":"Interrupted"}.get(row.get("status"),row.get("status","")),
                    "scientific_status":{"pass":"Recorded criteria met","warn":"Needs review","fail":"Needs review","not_applicable":"Not checked","not_checked":"Not checked","not_current":"Needs updating"}.get(row.get("scientific_status"),row.get("scientific_status","Not checked")),
                    "detail":row.get("error") or row.get("check_scope") or row.get("check") or row.get("validation") or ""})
            ui.label("Compare scores within the same target and setup. Scores from different targets do not establish selectivity.").classes("edu-info text-sm w-full")
            ui.table(columns=[{"name":"target","label":"Protein target","field":"target","align":"left","sortable":True},
                              {"name":"drug","label":"Drug molecule","field":"drug","align":"left","sortable":True},
                              {"name":"score","label":"Vina score (kcal/mol)","field":"score","sortable":True},
                              {"name":"status","label":"Status","field":"status","align":"left"},
                              {"name":"scientific_status","label":"Scientific checks","field":"scientific_status","align":"left"},
                              {"name":"detail","label":"Check or error","field":"detail","align":"left"}],
                     rows=normalized,row_key="pair_id",pagination={"rowsPerPage":20}).classes("w-full").props("wrap-cells flat")
            with ui.row().classes("gap-3 flex-wrap"):
                ui.button("Download results CSV",on_click=lambda:ui.download(engine.export_csv(active),filename="docking_results.csv")).props("unelevated color=primary no-caps icon=download")
                ui.button("Download complete study",on_click=download_package).props("outline no-caps icon=folder_zip")
                ui.button("Download simple docking report", on_click=lambda: ui.download(research.report_text(active).encode("utf-8"), filename="Docking_report.txt")).props("outline no-caps icon=description")
                ui.button("Download printable report", on_click=lambda: ui.download(research.report_html(active), filename="Docking_report.html")).props("outline no-caps icon=print")
                ui.button("References and methods",on_click=lambda:show_section("Reference")).props("outline no-caps icon=menu_book")
            with ui.card().classes("study-card w-full gap-4"):
                ui.label("Inspect docked poses and scientific checks").classes("text-lg font-semibold")
                inspect("06_dock")
            with ui.expansion("Read the docking report in simple words", icon="description").classes("study-card w-full"):
                ui.label(research.report_text(active)).classes("text-sm").style("white-space:pre-wrap;overflow-wrap:anywhere")
            render_code_exports()
            render_swissadme()
            ui.button("Return to docking settings",on_click=lambda:navigate(3)).props("flat no-caps icon=arrow_back")

        def render_home():
            ui.label("Home").classes("study-title")
            ui.label("Start a new docking project, or continue work saved in Projects.").classes("study-instruction")
            with ui.card().classes("study-card w-full gap-4"):
                ui.icon("biotech", color="primary", size="36px")
                ui.label("New drug and target docking").classes("text-xl font-semibold")
                ui.label("Add one or more protein targets and drug molecules. Follow five steps to prepare structures, run docking and review results.").classes("study-instruction")
                ui.button("Start new project", on_click=new_study).props("unelevated color=primary no-caps icon=add size=lg").set_enabled(not busy)
                ui.label("1. Targets  →  2. Drugs  →  3. Prepare  →  4. Dock  →  5. Results").classes("text-sm edu-muted")
            if active is not None or draft_started:
                with ui.card().classes("study-card w-full gap-3"):
                    ui.label(active.data.get("name", "Docking project") if active is not None else draft["name"]).classes("text-lg font-semibold")
                    ui.label(f"Step {current+1} of {len(STEPS)} · {STEPS[current][0]}" + (" · Saved locally" if active is not None else " · New project inputs")).classes("text-sm edu-muted")
                    ui.button("Continue project", on_click=continue_project).props("unelevated color=primary no-caps icon=arrow_forward")
            with ui.card().classes("study-card w-full gap-3"):
                ui.label("Try a ready example").classes("text-lg font-semibold")
                ui.label("Use the offline EGFR protein and example drugs. Preparation and docking run when you select those actions.").classes("text-sm edu-muted")
                ui.button("Use the offline EGFR example", on_click=example).props("outline color=primary no-caps icon=science").set_enabled(not busy)

        async def delete_saved_project(item, dialog):
            nonlocal active, current, selected_target, selected_drug, busy, message, draft_started
            if busy:
                return
            from core import project_library
            dialog.close()
            busy, message = True, "Deleting the selected project…"
            render()
            try:
                await asyncio.to_thread(project_library.delete_project, projects_root(), item["folder"])
                if active is not None and active.folder.resolve() == Path(item["folder"]).resolve():
                    active = None
                    current = 0
                    selected_target = selected_drug = ""
                    draft_started = False
                    draft.update(name="Docking study", targets=[], drugs=[], smiles="", lookup="")
                notify_success("Project deleted. Its storage space has been released.")
            except Exception as error:
                notify_error(error)
            finally:
                busy = False
                await load_library()
                render()
                if not dialog.is_deleted:
                    dialog.delete()

        def confirm_delete(item):
            from core import project_library
            if busy:
                return
            with client, ui.dialog() as dialog, ui.card().classes("gap-4 p-6").style("width:520px;max-width:94vw"):
                ui.label("Delete project").classes("text-xl font-semibold")
                ui.label(item["name"]).classes("font-semibold")
                ui.label(f"This permanently removes this project's inputs, structures and results ({project_library.human_size(item.get('size_bytes', 0))}). Export any results you need before deleting.").classes("text-sm edu-muted")
                ui.label("Other projects and the application will be kept.").classes("text-sm")
                with ui.row().classes("w-full justify-end gap-2"):
                    ui.button("Cancel", on_click=dialog.close).props("flat no-caps")
                    ui.button("Delete project", on_click=lambda:delete_saved_project(item, dialog)).props("unelevated color=negative no-caps icon=delete")
            dialog.on("hide", lambda:dialog.delete() if not dialog.is_deleted else None)
            dialog.open()

        def render_projects():
            from core import project_library
            with ui.row().classes("w-full items-center gap-3"):
                ui.label("Projects").classes("study-title")
                ui.space()
                ui.button("Refresh projects", on_click=load_library).props("flat no-caps icon=refresh").set_enabled(not library_loading)
                ui.button("Start new project", on_click=new_study).props("unelevated color=primary no-caps icon=add").set_enabled(not busy)
            ui.label("Open saved work or delete projects you no longer need.").classes("study-instruction")
            ui.label("Saved projects use disk space. Their structures and results are loaded only when opened. Deleting a project frees its files from disk.").classes("edu-info text-sm w-full")
            if library_loading:
                with ui.row().classes("items-center gap-3"):
                    ui.spinner(size="sm", color="primary")
                    ui.label("Loading saved projects and storage sizes…").classes("text-sm edu-muted")
                return
            if library_error:
                ui.label("Projects could not be loaded: " + library_error).classes("edu-info text-sm w-full")
                return
            if not library_rows:
                ui.label("No saved projects yet. Start a new project from Home.").classes("study-card w-full")
                return
            total = sum(row.get("size_bytes", 0) for row in library_rows)
            ui.label(f"{len(library_rows)} saved projects · {project_library.human_size(total)} on disk").classes("text-sm font-semibold")
            ui.label("Saved in: " + str(projects_root())).classes("text-xs edu-muted").style("overflow-wrap:anywhere")
            for item in library_rows:
                with ui.card().classes("study-card w-full gap-3").props('data-testid="project-card"') as card:
                    card._props["data-project-folder"] = Path(item["folder"]).name
                    with ui.row().classes("w-full items-start gap-3"):
                        ui.icon("folder", color="primary", size="26px")
                        with ui.column().classes("gap-1").style("flex:1;min-width:180px"):
                            ui.label(item["name"]).classes("text-lg font-semibold")
                            status = str(item.get("status", "saved")).replace("_", " ").capitalize()
                            ui.label(f"{status} · {project_library.human_size(item.get('size_bytes', 0))}").classes("text-sm edu-muted")
                            if item.get("target_count") or item.get("drug_count"):
                                ui.label(f"{item.get('target_count', 0)} targets · {item.get('drug_count', 0)} drugs").classes("text-xs edu-muted")
                            if item.get("updated"):
                                ui.label("Last saved: " + str(item["updated"]).replace("T", " ")[:19]).classes("text-xs edu-muted")
                        if item.get("kind") == "legacy":
                            from urllib.parse import quote
                            ui.button("Open project", on_click=lambda folder=Path(item["folder"]).name:ui.navigate.to("/legacy?project="+quote(folder))).props("outline color=primary no-caps icon=folder_open").set_enabled(not busy)
                        else:
                            ui.button("Open project", on_click=lambda folder=item["folder"]:resume(folder)).props("outline color=primary no-caps icon=folder_open").set_enabled(not busy)
                        ui.button("Delete project", on_click=lambda row=item:confirm_delete(row)).props("flat color=negative no-caps icon=delete_outline").set_enabled(not busy)

        def download_references():
            try:
                ui.download(references.render_references(opened_projects()).encode("utf-8"), filename="References_and_methods.txt")
            except Exception as error:
                notify_error(error)

        def render_reference():
            ui.label("Reference").classes("study-title")
            ui.label("Protein and drug downloads, scientific references and methods.").classes("study-instruction")
            if active is not None:
                ui.label("Current project: " + active.data.get("name", "Docking project")).classes("text-sm font-semibold")
                ui.button("Continue project", on_click=continue_project).props("outline no-caps icon=arrow_back")
            with ui.card().classes("study-card w-full gap-3"):
                ui.label("Protein and drug downloads").classes("text-xl font-semibold")
                for item in references.DOWNLOAD_LINKS:
                    ui.link(item.get("title", item.get("name", "Download source")), item["url"], new_tab=True).classes("font-medium")
                    if item.get("description"):
                        ui.label(item["description"]).classes("text-sm edu-muted -mt-2")
            with ui.card().classes("study-card w-full gap-3"):
                ui.label("Scientific references").classes("text-xl font-semibold")
                ui.label("The catalog covers the tools, databases and interpretation methods available in this application. The methods exported with each project identify the tools and data used in that run.").classes("text-sm edu-muted")
                ui.button("Download references and methods", on_click=download_references).props("unelevated color=primary no-caps icon=download")
                for reference in references.references():
                    with ui.column().classes("gap-1 py-2"):
                        ui.link(reference["title"], reference["url"], new_tab=True).classes("font-medium")
                        ui.label(reference["citation"]).classes("text-sm edu-muted")

        def render_about():
            ui.label("About").classes("study-title")
            with ui.card().classes("study-card w-full gap-4"):
                ui.image("/assets/logo.png").props('alt="Akshaisarathi Healthcare logo"').style("width:360px;max-width:100%;height:auto")
                ui.label("In-Silico Docking Studio").classes("text-xl font-bold")
                ui.label("Step-by-step molecular docking for medical education and in-silico studies.").classes("study-instruction")
                ui.separator()
                ui.label("Developed by").classes("edu-kicker")
                for developer in DEVELOPERS:
                    ui.label(developer).classes("text-base font-medium").style("line-height:1.7")
                ui.link("Akshaisarathi Healthcare", "https://akshaisarathi.com/", new_tab=True)
                ui.label("Scientific tools, methods and source links are listed in Reference and included with project exports.").classes("edu-info text-sm")

        def render():
            with client:
                nav.clear()
                body.clear()
                in_project = section == "Home" and workflow
                drawer.set_value(in_project)
                steps_toggle.set_visibility(in_project)
                for name, button in main_buttons.items():
                    button.style("background:#e6f2f1;color:#176b70" if section == name else "background:transparent;color:#4c6074")
                    button._props["aria-current"] = "page" if section == name else "false"
                    button.update()
                if in_project:
                  with nav:
                    ui.label("PROJECT STEPS").classes("edu-kicker mt-2 mb-2")
                    for index,(title,detail,icon) in enumerate(STEPS):
                        with ui.button(on_click=lambda i=index:navigate(i)).props("flat no-caps").classes("study-step").style("background:#e6f2f1;color:#176b70" if current==index else "color:#4c6074") as button:
                            ui.icon(icon,size="22px")
                            with ui.column().classes("gap-0"):
                                ui.label(f"{index+1}. {title}").classes("font-semibold")
                                ui.label(detail).classes("text-xs").style("line-height:1.4")
                        button.set_enabled(not busy and (index < 2 or active is not None))
                    ui.separator().classes("my-3")
                    if active is not None:
                        ui.label(active.data.get("name","Docking study")).classes("font-semibold text-sm")
                        ui.label("Saved locally · Reopen from Projects").classes("edu-status-note")
                        ui.button("Copy inputs to a new study",on_click=lambda:new_study(True)).props("flat no-caps icon=content_copy align=left").classes("w-full text-xs").set_enabled(not busy)
                    ui.button("Back to Home",on_click=lambda:show_section("Home")).props("outline no-caps icon=home").classes("w-full")
                    ui.label("Follow the steps to prepare, dock and review the study.").classes("edu-status-note mt-3")
                with body:
                    if not in_project:
                        {"Home":render_home, "Projects":render_projects, "Reference":render_reference, "About":render_about}[section]()
                    else:
                        with ui.row().classes("w-full items-center gap-3"):
                            ui.label(active.data.get("name", "Docking project") if active is not None else draft["name"]).classes("text-base font-semibold")
                            ui.space()
                            ui.label(f"Step {current+1} of {len(STEPS)}").classes("study-status").props('data-testid="project-step"')
                        ui.label(STEPS[current][0]).classes("study-title")
                        if busy:
                            ui.label("Processing the current action. The saved project and previous results are retained.").classes("edu-info text-sm w-full")
                            summary()
                        else:
                            [render_targets,render_drugs,render_prepare,render_dock,render_results][current]()
                refresh_progress()

        render()
        if section == "Projects":
            ui.timer(0.01, load_library, once=True)

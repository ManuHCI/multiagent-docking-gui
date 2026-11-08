"""Shared, local-only visual language and help for the teaching interface."""
from __future__ import annotations

from nicegui import ui

DEVELOPERS = (
    "Dr Manu Kumar Shetty, Professor, Department of Pharmacology, Maulana Azad Medical College, New Delhi",
    "Dr Pratibha Nadig, Founder, Akshaisarathi Healthcare",
)


def about_dialog():
    """Shared credits in every working interface."""
    with ui.dialog() as dialog, ui.card().classes("gap-4 p-7").style("width:700px;max-width:94vw"):
        ui.image("/assets/logo.png").props('alt="Akshaisarathi Healthcare logo"').style("width:360px;max-width:100%;height:auto")
        ui.label("In-Silico Docking Studio").classes("text-2xl font-bold")
        ui.label("Step-by-step molecular docking for medical education and in-silico studies.").classes("edu-muted")
        ui.separator()
        ui.label("Developed by").classes("edu-kicker")
        for developer in DEVELOPERS:
            ui.label(developer).classes("text-base font-medium").style("line-height:1.7")
        ui.link("Akshaisarathi Healthcare", "https://akshaisarathi.com/", new_tab=True)
        ui.label("Scientific tools, methods and source links are listed under Reference and included with study exports.").classes("edu-info text-sm")
        ui.button("Close", on_click=dialog.close).props("unelevated no-caps color=primary").classes("self-end")
    return dialog


def resources_dialog(projects=None):
    """Official download links, references and current-study provenance."""
    from core.references import DOWNLOAD_LINKS, references, render_references
    with ui.dialog() as dialog, ui.card().classes("gap-4 p-6").style("width:900px;max-width:96vw;max-height:90vh;overflow:auto"):
        with ui.row().classes("w-full items-center justify-between"):
            ui.label("References & downloads").classes("text-2xl font-bold")
            ui.button(icon="close", on_click=dialog.close).props("flat round").tooltip("Close references")
        ui.label("Download protein structures and drug records from their original databases. Review the structure, chain and binding site before docking.").classes("edu-info text-sm")
        ui.label("Protein and drug downloads").classes("text-lg font-semibold")
        for item in DOWNLOAD_LINKS:
            if isinstance(item, dict):
                ui.link(item.get("title", item.get("name", "Download source")), item["url"], new_tab=True)
                if item.get("description"):
                    ui.label(item["description"]).classes("text-sm edu-muted -mt-2")
        ui.separator()
        ui.label("Scientific references").classes("text-lg font-semibold")
        ui.label("The catalog covers the tools, databases and interpretation methods available in this application. The exported methods record identifies the tools and data actually used in each run.").classes("text-sm edu-muted")
        for reference in references():
            with ui.column().classes("w-full gap-1 py-2"):
                ui.link(reference["title"], reference["url"], new_tab=True).classes("font-medium")
                ui.label(reference["citation"]).classes("text-sm edu-muted")
        def download():
            current = projects() if callable(projects) else projects
            try:
                text = render_references(current)
                ui.download(text.encode("utf-8"), filename="References_and_methods.txt")
            except Exception as exc:
                ui.notify(f"References could not be exported: {exc}", type="negative")
        ui.button("Download references and methods", on_click=download).props("unelevated icon=download no-caps color=primary")
    return dialog


def apply_theme() -> None:
    """Keep both learning routes readable and consistent on laptop screens."""
    ui.colors(primary="#087c80", secondary="#173c59", accent="#d99935",
              positive="#26754d", warning="#a86a10", negative="#b63d48")
    ui.add_head_html('''<style>
    :root{--edu-ink:#17344b;--edu-muted:#5d7282;--edu-teal:#087c80;
      --edu-line:#dce7ec;--edu-bg:#f3f7f9}
    body{font-family:'Segoe UI',system-ui,sans-serif;background:var(--edu-bg);
      color:var(--edu-ink);font-size:15px;-webkit-font-smoothing:antialiased}
    .nicegui-content{padding:0}.q-header{box-shadow:none!important}
    .q-card{border-radius:16px}.q-btn{border-radius:9px;letter-spacing:0}
    .q-btn:focus-visible,.q-field:focus-within{outline:3px solid #55b4bf;
      outline-offset:3px}.q-field--outlined .q-field__control{border-radius:10px}
    .q-field--outlined .q-field__control:before{border-color:#bfd0db}
    .q-table__container{border-radius:12px;border:1px solid var(--edu-line);
      box-shadow:none}.q-table thead{background:#edf5f6}
    .q-table th{font-weight:650;color:#264d62}.q-table td{font-size:14px}
    .edu-header{min-height:78px;background:#fff!important;color:var(--edu-ink)!important;
      border-bottom:1px solid var(--edu-line);padding:12px 22px}
    .edu-brand{background:#e4f3f2;border:1px solid #cbe7e5;border-radius:13px;
      width:45px;height:45px;display:flex;align-items:center;justify-content:center}
    .edu-kicker{font-size:11px;font-weight:750;letter-spacing:1.5px;color:#54717f;
      text-transform:uppercase}.edu-muted{color:var(--edu-muted)}
    .edu-hero{background:linear-gradient(125deg,#123c53 0%,#15556a 65%,#107a7e 100%);
      padding:34px;border-radius:22px;color:#fff;position:relative;overflow:hidden}
    .edu-hero h1{font-size:36px;font-weight:700;line-height:1.16;letter-spacing:-1px;
      margin:0;max-width:720px}.edu-hero .edu-subtitle{color:#d9eaf0;font-size:17px;
      line-height:1.65;max-width:660px}.edu-hero .q-badge{background:#ffffff1c!important;
      border:1px solid #ffffff35;color:#e3f4f4!important;padding:7px 11px}
    .edu-hero .q-btn{background:#fff!important;color:#135f68!important;
      min-height:48px;font-weight:650;padding:0 20px}
    .edu-feature{padding:22px;min-width:200px;flex:1;gap:12px;
      border:1px solid var(--edu-line);box-shadow:none!important;background:#fff}
    .edu-feature-icon{background:#eaf5f5;color:#087c80;width:44px;height:44px;
      border-radius:12px;display:flex;align-items:center;justify-content:center}
    .edu-info{background:#ecf5f7;border:1px solid #d8e9ee;border-radius:12px;
      padding:14px 16px}.edu-phase{font-size:11px;font-weight:750;letter-spacing:1px;
      color:#6d8390;text-transform:uppercase;margin:18px 0 5px}
    .edu-small-button{font-size:13px!important}.edu-home-card{border:1px solid #dce7ec;
      box-shadow:none!important;padding:22px}.edu-help .q-item{padding:12px 0}
    .edu-help p{line-height:1.65}.edu-resume{border:1px solid #dce7ec;border-radius:12px;
      background:#fff;min-height:54px;padding:9px 14px}
    .ws-card{border:1px solid #dce7ec;box-shadow:none!important;border-radius:16px!important}
    .ws-action{background:#087c80!important;color:#fff!important;min-height:44px}
    .ws-note{font-size:16px;line-height:1.7}.ws-nav .q-btn__content{flex-wrap:nowrap;
      justify-content:flex-start;text-align:left;gap:10px}
    .ws-nav{font-size:13px!important;line-height:1.4!important;min-height:45px;
      padding:8px 10px!important}.ws-nav .q-icon{font-size:19px!important;flex-shrink:0}
    .ws-nav:disabled{opacity:.6!important}.ws-title{letter-spacing:-.4px}
    .edu-status-note{font-size:12px;line-height:1.5;color:#5c7482}
    @media(max-width:1100px){.edu-header{padding:10px 14px}.edu-header .edu-tagline{display:none}}
    @media(max-width:700px){.edu-hero{padding:24px}.edu-hero h1{font-size:28px}
      .ws-title{font-size:25px!important}.ws-columns{flex-direction:column!important}
      .edu-header{min-height:68px}.edu-header .edu-brand{width:36px;height:36px}
      .edu-header .edu-brand-text{font-size:15px!important}}
    </style>''')


def learning_help():
    """An always-available glossary with an explanation of result labels."""
    from core.workshop import GLOSSARY
    with ui.dialog() as dialog, ui.card().classes("edu-help gap-4 p-6").style("width:680px;max-width:94vw"):
        with ui.row().classes("w-full items-center justify-between"):
            ui.label("A quick guide to molecular docking").classes("text-xl font-bold")
            ui.button(icon="close", on_click=dialog.close).props("flat round").tooltip("Close guide")
        ui.label("Add protein targets and drug molecules, prepare each target's binding site, run docking and review the checks. No quizzes or written answers are required.").classes("ws-note edu-muted")
        with ui.expansion("Result labels", value=True).classes("w-full"):
            for icon, colour, title, text in [
                ("check_circle", "#26754d", "Checked", "The stated check passed. Read its scope and limitations."),
                ("warning", "#a86a10", "Warning", "The calculation needs careful interpretation. Read the suggested next action."),
                ("cancel", "#b63d48", "Check failed", "A scientific check failed, even if the calculation finished."),
                ("update", "#62798a", "Needs updating", "Inputs changed. Run the step again before using the old result."),
            ]:
                with ui.row().classes("no-wrap items-start gap-3 py-2"):
                    ui.icon(icon).style(f"color:{colour};font-size:22px")
                    with ui.column().classes("gap-0"):
                        ui.label(title).classes("font-semibold")
                        ui.label(text).classes("text-sm edu-muted")
        with ui.expansion("Explain the key terms").classes("w-full"):
            for term, text in GLOSSARY.items():
                with ui.column().classes("gap-0 py-2"):
                    ui.label(term).classes("font-semibold")
                    ui.label(text).classes("text-sm edu-muted")
        ui.label("For medical education. A docking pose or score is a hypothesis to investigate; it does not establish a medicine's potency, benefit or safety.").classes("edu-info text-sm")
        ui.button("Continue learning", on_click=dialog.close).props("unelevated no-caps color=primary").classes("self-end")
    return dialog

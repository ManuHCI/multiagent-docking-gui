"""
The interface.

A desktop window, written in Python, built around one idea: the user is
a pharmacologist who has never done this before, and the thing that
will decide whether this software helps them is not whether it can
dock a molecule. It is whether they understand what came out.

So every stage screen has the same four parts, in the same places:

    WHAT THIS STEP DOES     one sentence, always visible
    EXPLAIN THIS STEP       the full teaching text, on demand
    DO IT                   the controls and the button
    WHAT CAME OUT           the result, and beneath it the check that
                            says whether the result means anything

The fourth part is never optional and never collapsed by default. A
result without its control is the thing this software exists to
prevent.

Run with:  python app.py
"""

from __future__ import annotations

import asyncio
import json
from pathlib import Path
from typing import Any, Callable

from nicegui import app as nicegui_app
from nicegui import ui

from core.contract import Verdict
from core.education_ui import apply_theme, about_dialog as make_about_dialog, resources_dialog as make_resources_dialog
from core.project import Project
from forms import StageForm

# Where projects are saved. Read through core.settings rather than
# fixed here, because on a departmental Windows machine the home
# directory is on a small C: drive and the user needs to be able to
# point this at a bigger one without editing code.
from core.settings import (                                   # noqa: E402
    default_root, problem_with_root, projects_root, set_projects_root,
)


def list_projects() -> list[dict[str, Any]]:
    """Every project on this computer, most recently used first."""
    import datetime

    out: list[dict[str, Any]] = []
    root = projects_root()
    if not root.exists():
        return out
    for folder in root.iterdir():
        db = folder / "project.db"
        if not db.exists():
            continue
        try:
            when = datetime.datetime.fromtimestamp(db.stat().st_mtime)
        except OSError:
            continue
        done = 0
        try:
            import sqlite3
            conn = sqlite3.connect(db)
            done = conn.execute(
                "SELECT COUNT(*) FROM steps WHERE status='complete'"
            ).fetchone()[0]
            conn.close()
        except Exception:  # noqa: BLE001
            pass
        out.append({"name": folder.name, "folder": folder,
                    "when": when, "steps": done})
    out.sort(key=lambda p: p["when"], reverse=True)
    return out


ASSETS = Path(__file__).parent / "assets"
nicegui_app.add_static_files("/assets", str(ASSETS))

# Branding and attribution. Kept here as constants rather than scattered
# through the layout, so they can be changed in one place -- including
# when this is eventually served from akshaisarathi.com rather than run
# on a laptop.
SITE_URL = "https://akshaisarathi.com/"
AUTHOR = ("Dr Manu Kumar Shetty, Professor, Department of Pharmacology, "
          "Maulana Azad Medical College, New Delhi")
GUIDE = "Dr Pratibha Nadig, Founder, Akshaisarathi Healthcare"

# ----------------------------------------------------------------------
# Look. Restrained on purpose: this is a teaching tool used by tired
# people in a departmental computer room, not a dashboard.
# ----------------------------------------------------------------------

VERDICT_STYLE: dict[str, tuple[str, str, str]] = {
    # colour, icon, heading
    "pass": ("#2e7d32", "check_circle", "Checked"),
    "warn": ("#b26a00", "warning", "Read this before going on"),
    "fail": ("#b3261e", "cancel", "Stop"),
    "not_applicable": ("#5f6368", "help", "Not checked"),
}

STAGES: list[tuple[str, str, str, str, str]] = [
    # id, title, one-line blurb, icon, colour
    #
    # The colour is not decoration. Seven identical grey steps give a
    # novice nothing to navigate by; a colour and a picture per step
    # mean "I am on the green one, molecules" without reading. The
    # family runs blue -> teal -> green -> amber -> purple, which is
    # the order of the work, so position is visible at a glance.
    ("01_target", "Target",
     "Pick the protein and review existing target evidence.",
     "biotech", "#1a73e8"),
    ("02_structure", "Structure",
     "Get its 3D shape from the world database.",
     "view_in_ar", "#00838f"),
    ("03_receptor", "Prepare protein",
     "Clean the structure without throwing away what matters.",
     "cleaning_services", "#00897b"),
    ("04_site", "Binding site",
     "Tell the computer which pocket to aim at.",
     "my_location", "#2e7d32"),
    ("05_ligands", "Prepare molecules",
     "Turn your drugs into proper 3D molecules.",
     "science", "#f9a825"),
    ("06_dock", "Dock",
     "Predict poses and read the setup, baseline and size checks.",
     "join_inner", "#e8710a"),
    ("08_admet", "Pharmacokinetics",
     "Calculate descriptors and explore ADME questions.",
     "medication", "#7b1fa2"),
    ("09_md", "Simulation (MD) · optional",
     "Explore motion under a model. Prepare a cloud job.",
     "waves", "#c2185b"),
]


class Session:
    """What the window is currently showing."""

    def __init__(self) -> None:
        self.project: Project | None = None
        self.current: str = "01_target"
        self.busy = False
        self.render: Callable | None = None


class ClientSession:
    """Keep projects and navigation separate for each browser client."""
    def _get(self) -> Session:
        storage = nicegui_app.storage.client
        if "teaching_session" not in storage:
            storage["teaching_session"] = Session()
        return storage["teaching_session"]

    def __getattr__(self, name):
        return getattr(self._get(), name)

    def __setattr__(self, name, value):
        setattr(self._get(), name, value)


S = ClientSession()


# ----------------------------------------------------------------------
# Reusable pieces
# ----------------------------------------------------------------------

def section_heading(text: str, icon: str, colour: str) -> None:
    """
    A coloured strip that separates the parts of a stage screen.

    Grey capitals all look alike. Tying each section to the step's own
    colour makes the page scannable, which matters most to the reader
    who has scrolled past the explanation and is looking for the
    button.
    """
    with ui.row().classes("items-center no-wrap q-mt-md q-mb-xs"):
        ui.icon(icon).style(f"color:{colour}; font-size:1.15rem")
        ui.label(text).style(
            f"color:{colour}; font-weight:700; font-size:0.78rem; "
            "letter-spacing:0.09em")


def control_panel(control: Any) -> None:
    """
    The check that says whether the result means anything.

    Always rendered when a control exists, never behind a toggle, and
    placed directly under the result rather than at the bottom of the
    page. Specification Section 9: the controls cannot be switched off.
    """
    if control is None:
        return
    verdict = (control.verdict.value if hasattr(control.verdict, "value")
               else str(control.verdict))
    colour, icon, heading = VERDICT_STYLE.get(verdict,
                                              VERDICT_STYLE["not_applicable"])

    with ui.card().classes("w-full q-mt-sm").style(
        f"border-left: 6px solid {colour}; background: {colour}0d"
    ):
        with ui.row().classes("items-center no-wrap"):
            ui.icon(icon).style(f"color:{colour}; font-size:1.5rem")
            ui.label(heading).style(
                f"color:{colour}; font-weight:600; font-size:1.05rem")
            ui.space()
            ui.label("Scientific check" if "?" in control.name else control.name).style("color:#5f6368; font-size:0.85rem")
        ui.label(control.plain_verdict).classes("q-mt-xs").style(
            "line-height:1.55")

        # The action. A finding with no next step is only a worry, and a
        # student who does not know what to do about it will either
        # ignore the check or stop altogether.
        if getattr(control, "what_to_do", ""):
            with ui.row().classes("no-wrap items-start q-mt-sm q-pa-sm").style(
                "background:#ffffff; border-radius:6px; "
                f"border:1px solid {colour}44; width:100%"
            ):
                ui.icon("arrow_forward").style(
                    f"color:{colour}; font-size:1.2rem; margin-top:2px")
                with ui.column().classes("q-gutter-none").style("flex:1"):
                    ui.label("WHAT TO DO").style(
                        f"color:{colour}; font-weight:700; "
                        "font-size:0.72rem; letter-spacing:0.08em")
                    ui.label(control.what_to_do).style(
                        "line-height:1.55; color:#202124")

        if control.detail:
            with ui.expansion("The numbers behind this").classes("w-full"):
                with ui.column().classes("q-gutter-none"):
                    for key, value in control.detail.items():
                        with ui.row().classes("no-wrap"):
                            ui.label(f"{key}:").style(
                                "color:#5f6368; min-width:230px")
                            ui.label(str(value))
        # Only shown when there is no explicit action, because the
        # WHAT TO DO block above already says what to do and a vague
        # "read the above first" underneath it contradicts a Stop.
        if control.blocks_progress and not getattr(control, "what_to_do", ""):
            ui.label(
                "Read the above before going on - it changes what your "
                "result is worth."
            ).classes("q-mt-sm").style(
                f"color:{colour}; font-style:italic; font-size:0.9rem")


def guard(fn: Callable[..., Any], what: str) -> Callable[..., Any]:
    """
    Wrap a click handler so it can never fail in silence.

    NiceGUI swallows an exception raised inside a handler: the click
    does nothing, no message appears, nothing reaches the console. That
    happened here - a sidebar item called a function that had been
    deleted, and the only symptom was that clicking did nothing at all.

    For this audience that is the worst possible failure. A resident
    who clicks and gets no response has no way to tell whether the
    software is broken, busy, or waiting for something they have not
    done, and will usually assume they did something wrong. An ugly
    error is far better than none.
    """
    def wrapped(*args: Any, **kwargs: Any) -> Any:
        try:
            return fn(*args, **kwargs)
        except Exception as exc:  # noqa: BLE001
            ui.notify(
                f"Something went wrong while {what}. This is a fault in "
                f"the software, not in anything you did. "
                f"({type(exc).__name__}: {str(exc)[:120]})",
                type="negative", timeout=12000, multi_line=True,
            )
            import traceback
            traceback.print_exc()
    return wrapped


def teaching_panel(teaching: Any, accent: str = "#2e5496") -> None:
    """
    Plain words first, detail only if asked for.

    The order matters. An earlier version led with a paragraph written
    for someone who already knew what a crystal structure was, and a
    pharmacology professor reading the protein-preparation step said
    even MD students would not follow it. So the plainest version is
    now what you see; everything else is a click away for the reader
    who wants it.
    """
    if teaching is None:
        return

    # The one-line purpose is already on the coloured banner above, so
    # it is not repeated here.
    simple = getattr(teaching, "simple", "") or ""
    if simple:
        with ui.card().classes("w-full").style(
            f"background:{accent}0f; border-left:6px solid {accent}; "
            "box-shadow:none"
        ):
            with ui.row().classes("items-center no-wrap q-mb-xs"):
                ui.icon("school").style(f"color:{accent}")
                ui.label("READ THIS FIRST - what this step is, in plain "
                         "words").style(
                    f"font-weight:700; color:{accent}; font-size:0.8rem; "
                    "letter-spacing:0.06em")
            for para in simple.split("\n\n"):
                ui.label(para).classes("q-mb-sm").style(
                    "line-height:1.7; font-size:1.02rem; color:#20303f")

    with ui.expansion(
        "Tell me more (the proper explanation)"
    ).classes("w-full q-mt-sm").props("dense"):
        for para in teaching.panel.split("\n\n"):
            ui.label(para).classes("q-mb-sm").style("line-height:1.65")

    if teaching.not_mean:
        with ui.card().classes("w-full q-mt-sm").style(
            "background:#fff8e1; border-left:6px solid #c47f00; "
            "box-shadow:none"
        ):
            with ui.row().classes("items-center no-wrap"):
                ui.icon("report_problem").style("color:#c47f00")
                ui.label("Interpretation limits").style(
                    "font-weight:600; color:#7a3b00")
            for item in teaching.not_mean:
                with ui.row().classes("no-wrap items-start q-mt-xs"):
                    ui.label("—").style("color:#c47f00; min-width:14px")
                    ui.label(item).style("line-height:1.55; color:#5a4a2a")

    if teaching.learn_more:
        with ui.row().classes("q-mt-sm q-gutter-sm items-center"):
            ui.label("Read more elsewhere:").style(
                "color:#5f6368; font-size:0.9rem")
            for label, url in teaching.learn_more:
                ui.link(label, url, new_tab=True).style("font-size:0.9rem")


def stage_status(stage_id: str) -> tuple[str, str]:
    """Icon and colour for the sidebar."""
    if S.project is None:
        return "radio_button_unchecked", "#9aa0a6"
    icon, colour, _ = S.project.display_status(stage_id)
    return icon, colour


# ----------------------------------------------------------------------
# The page
# ----------------------------------------------------------------------

WORKFLOW_HELP = [
    ("Purpose of the software",
     "It walks you through designing a drug on a computer, one step at "
     "a time, and - just as importantly - it tells you how much each "
     "result is actually worth. Computers are good at producing "
     "confident-looking numbers in this field. Most of those numbers "
     "mean far less than they appear to, and this software is built to "
     "say so."),
    ("Start with simple instructions",
     "Every step explains itself "
     "in plain words before you run it. The guided docking workspace provides "
     "a prepared EGFR example, 3D views and practical controls. In the "
     "regular workflow, review the inputs and checks at each step before "
     "drawing a conclusion."),
    ("Result labels and colours",
     "BLUE box - what this step is, in plain words. Read it first.\n\n"
     "AMBER box, 'What this does NOT mean' - the wrong conclusions "
     "people commonly draw from this kind of result. This is the most "
     "useful part of the whole software.\n\n"
     "After you run a step, a box appears under the result with a "
     "green tick, an amber warning, or a red stop sign. That is the "
     "check. Inside it there is always a line headed WHAT TO DO, which "
     "tells you the next action - continue, go back and change "
     "something, or stop."),
    ("Scientific checks",
     "Each step verifies its own work. The docking step does two "
     "things: it takes a drug whose real position is already known and "
     "asks the software to find it again, and it docks a small comparison "
     "set matched for size and lipophilicity. Failure to separate from "
     "that baseline weakens a score-based claim; it does not establish "
     "why the scores are similar. These checks are attempted automatically "
     "and report when a required reference or comparison is unavailable."),
    ("The eight steps, in order",
     "1. Target - pick the protein, and find out whether anyone has "
     "ever made a drug against it.\n"
     "2. Structure - get its 3D shape from the world database.\n"
     "3. Prepare protein - clean the structure without throwing away "
     "a metal or cofactor that the protein needs.\n"
     "4. Binding site - draw the box that says where to look.\n"
     "5. Prepare molecules - turn your drugs into 3D molecules.\n"
     "6. Dock - inspect the pose, re-docking, baseline and size checks.\n"
     "7. Pharmacokinetics - calculate molecular descriptors and review ADME properties.\n"
     "8. Simulation (optional) - prepare a job to explore motion under "
     "a force-field model. It runs separately, after this workshop."),
    ("Saved project files",
     "Into one folder on this computer. You can copy it, email it, or "
     "hand it to a supervisor, and it contains every setting used so "
     "the work can be repeated or checked."),
]


@ui.page("/station/swissadme")
def swissadme_launcher(project: str = "") -> None:
    """
    Serve the prepared handoff page in a tab of its own.

    Downloading an HTML file and asking a novice to find it in their
    downloads folder is a step too many. Serving it means one click to
    get there, and another - theirs, on that page - to submit. The
    submission stays user-initiated, which is what SwissADME's terms
    require.
    """
    from stations.swissadme import SwissADMEStation

    folder = (projects_root() / project).resolve()
    if not project or folder.parent != projects_root().resolve() or not (folder / "project.db").exists():
        ui.label("No project is open.")
        return
    selected = Project.open(folder)
    station = SwissADMEStation(selected)
    try:
        page = station.build_launcher()
        ui.add_body_html(page.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001
        ui.label(f"The handoff could not be prepared: {exc}")
    finally:
        selected.close()


@ui.page("/")
def home_page() -> None:
    """All launchers enter the same four-section application."""
    ui.navigate.to("/study")


@ui.page("/legacy")
def main_page(project: str = "") -> None:
    """Keep existing individual projects accessible from the Projects library."""
    root = projects_root().resolve()
    folder = (root / project).resolve()
    if not project or folder.parent != root or not (folder / "project.db").is_file():
        ui.navigate.to("/study?section=projects")
        return
    from core.project_library import project_lock
    legacy_lock = project_lock(folder)
    if not legacy_lock.acquire(blocking=False):
        ui.label("This project is already open. Close its other window before opening it here.")
        ui.link("Back to Projects", "/study?section=projects")
        return
    try:
        S.project = Project.open(folder)
    except Exception:
        legacy_lock.release()
        raise
    selected_project = S.project
    released = False

    def release_project():
        nonlocal released
        if not released:
            selected_project.close()
            legacy_lock.release()
            released = True

    def leave_project(destination):
        release_project()
        ui.navigate.to("/study?section=" + destination)

    S.current = "01_target"
    ui.context.client.on_disconnect(release_project)
    apply_theme()
    with ui.dialog() as help_dialog, ui.card().style("max-width:760px"):
        ui.label("How this works").style(
            "font-size:1.35rem; font-weight:600; color:#1f3864")
        ui.label(
            "Written for someone who has never done this before. "
            "Nothing here assumes you know any chemistry software."
        ).style("color:#5f6368")
        ui.separator().classes("q-my-sm")
        for question, answer in WORKFLOW_HELP:
            with ui.expansion(question).classes("w-full"):
                for para in answer.split("\n"):
                    if para.strip():
                        ui.label(para).classes("q-mb-xs").style(
                            "line-height:1.7")
        ui.button("Close", on_click=lambda: help_dialog.close()).props(
            "flat").classes("q-mt-md")

    about_dialog = make_about_dialog()
    references_dialog = make_resources_dialog(lambda: [S.project] if S.project else [])

    with ui.header().classes("edu-header"):
        with ui.row().classes("items-center no-wrap w-full"):
            ui.button(icon="menu", on_click=lambda: drawer.toggle()).props(
                "flat round").classes("lt-lg").tooltip("Show workflow steps")
            ui.image("/assets/logo.png").props('alt="Akshaisarathi Healthcare logo"').style("width:175px;height:auto").classes("gt-sm")
            with ui.column().classes("q-gutter-none"):
                ui.label("In-Silico Docking Studio").classes(
                    "edu-brand-text text-lg font-bold")
                ui.label("MEDICAL EDUCATION").classes(
                    "edu-tagline edu-kicker").style("font-size:10px;letter-spacing:1px")
            ui.space()
            def close_project() -> None:
                if S.busy:
                    ui.notify("Wait for the current calculation to finish.", type="info")
                    return
                if S.project is not None:
                    S.project.close()
                S.project = None
                project_label.text = "no project open"
                S.current = "01_target"
                refresh_sidebar()
                render()

            project_label = ui.label("no project open").classes("gt-sm").style(
                "font-size:0.8rem; color:#5d7282;max-width:160px;overflow:hidden;text-overflow:ellipsis")
            for label, section, icon_name in (("Home", "home", "home"), ("Projects", "projects", "folder"),
                                              ("Reference", "reference", "menu_book"), ("About", "about", "info_outline")):
                ui.button(label, on_click=guard(lambda destination=section: leave_project(
                    destination), "leaving the project")).props(
                        f"flat icon={icon_name} no-caps").classes("edu-small-button")

    with ui.left_drawer(fixed=True).props("width=284 breakpoint=1024").style(
        "background:#ffffff;border-right:1px solid #dce7ec;padding:18px 12px"
    ) as drawer:
        ui.label("YOUR INVESTIGATION").classes("edu-kicker mb-2")
        ui.label("Follow the steps. Read each check.").classes("text-xs edu-muted mb-4")
        stage_rows: dict[str, Any] = {}
        for n, (stage_id, title, blurb, icon_name, accent) in enumerate(
                STAGES, start=1):
            with ui.row().classes(
                "items-center no-wrap cursor-pointer q-pa-sm q-mb-xs"
            ).style(
                f"border-radius:8px; background:#ffffff; "
                f"border-left:5px solid {accent}"
            ) as row:
                # The step's own colour and picture, so the sidebar can
                # be navigated by sight rather than by reading seven
                # similar lines of text.
                with ui.element("div").style(
                    f"background:{accent}1a; border-radius:50%; "
                    "width:34px; height:34px; display:flex; "
                    "align-items:center; justify-content:center; "
                    "flex:0 0 34px"
                ):
                    ui.icon(icon_name).style(
                        f"color:{accent}; font-size:1.2rem")
                with ui.column().classes("q-gutter-none").style("flex:1"):
                    ui.label(f"{n}. {title}").style(
                        "font-weight:600; font-size:0.95rem; color:#202124")
                    ui.label(blurb).style(
                        "font-size:0.75rem; color:#5f6368; line-height:1.3")
                    status_text = ui.label("Not run").style(
                        "font-size:0.7rem; line-height:1.3; color:#5f6368")
                icon, colour = stage_status(stage_id)
                ic = ui.icon(icon).style(f"color:{colour}; font-size:1.1rem")
                with ic:
                    status_tooltip = ui.tooltip("Not run")
                row.on("click", guard(lambda s=stage_id: select(s),
                                      "opening that step"))
                stage_rows[stage_id] = (row, ic, status_text, status_tooltip)

    body = ui.column().classes("w-full p-5 md:p-8").style("max-width:1040px;margin:0 auto;gap:0")

    from types import SimpleNamespace
    progress = SimpleNamespace(text="")
    with ui.row().classes("q-px-lg items-center q-gutter-sm"):
        spinner = ui.spinner(size="1.4rem")
        spinner.set_visibility(False)

    def select(stage_id: str) -> None:
        """Jump straight to a step from the sidebar."""
        if S.busy:
            ui.notify("Wait for the current calculation to finish.", type="info")
            return
        S.current = stage_id
        render()

    def go(delta: int) -> None:
        if S.busy:
            ui.notify("Wait for the current calculation to finish.", type="info")
            return
        ids = [st[0] for st in STAGES]
        i = ids.index(S.current) + delta
        if 0 <= i < len(ids):
            S.current = ids[i]
            render()

    with ui.footer().style(
        "background:#ffffff; border-top:1px solid #d8dee4; padding:10px 18px"
    ):
        with ui.row().classes("items-center no-wrap w-full q-gutter-sm"):
            back_btn = ui.button("Back", on_click=guard(lambda: go(-1), "going back")).props(
                "flat icon=arrow_back")
            ui.space()
            step_label = ui.label().style("color:#5f6368; font-size:0.9rem")
            ui.space()
            next_btn = ui.button("Next step",
                                 on_click=guard(lambda: go(1),
                                                "moving on")).props(
                "unelevated color=primary icon-right=arrow_forward")

    async def do_run(stage_id: str, form: StageForm, button: Any) -> None:
        """Run one stage without freezing the window."""
        if S.busy:
            return
        project = S.project
        stage = build_stage(stage_id)
        if stage is None:
            ui.notify("This step is not available yet.", type="warning")
            return

        state = stage.check()
        if state.readiness.value == "blocked":
            ui.notify(f"{state.reason} {state.fix}", type="negative",
                      timeout=8000)
            return
        if state.readiness.value == "needs_input":
            ui.notify(f"{state.reason} {state.fix}", type="warning",
                      timeout=8000)
            return

        button.disable()
        S.busy = True
        spinner.set_visibility(True)
        progress.text = "Working. Docking can take a few minutes."
        try:
            values = form.values()
            result = await asyncio.to_thread(lambda: stage.run(**values))
        except Exception as exc:  # noqa: BLE001
            ui.notify(f"{type(exc).__name__}: {exc}", type="negative",
                      timeout=10000)
            return
        finally:
            spinner.set_visibility(False)
            progress.text = ""
            button.enable()
            S.busy = False

        project.save_step(stage_id, stage.title, result)
        if result.success:
            ui.notify(result.summary or "Done", type="positive")
        else:
            ui.notify(result.message or "That did not work",
                      type="negative", timeout=10000)
        refresh_sidebar()
        render()

    def refresh_sidebar() -> None:
        accents = {st[0]: st[4] for st in STAGES}
        for stage_id, (row, icon_el, status_text, status_tooltip) in stage_rows.items():
            icon, colour = stage_status(stage_id)
            status = S.project.display_status(stage_id)[2] if S.project else "Not run"
            icon_el.name = icon
            icon_el.style(f"color:{colour}; font-size:1.1rem")
            status_text.text = status
            status_text.style(f"color:{colour}")
            status_tooltip.text = status
            # Where am I? Without this the sidebar shows seven steps and
            # no answer to that question.
            here = stage_id == S.current
            accent = accents[stage_id]
            row.style(
                "border-radius:8px; "
                + (f"background:{accent}1f; border-left:5px solid {accent}; "
                   "box-shadow:0 1px 3px rgba(0,0,0,0.12)"
                   if here else
                   f"background:#ffffff; border-left:5px solid {accent}")
            )

    def render() -> None:
        body.clear()
        # The sidebar must agree with the page about which step this is,
        # so the highlight is repainted on every render rather than only
        # when a step finishes.
        refresh_sidebar()
        ids = [st[0] for st in STAGES]
        i = ids.index(S.current)
        _, title, blurb, icon_name, accent = STAGES[i]
        if S.project is None:
            # No project open: the heading should say so rather than
            # naming a step the user cannot reach yet.
            title = "Your project workspace"
            step_label.text = ""
            back_btn.set_visibility(False)
            next_btn.set_visibility(False)
        else:
            step_label.text = f"Step {i + 1} of {len(STAGES)}: {title}"
            back_btn.set_visibility(i > 0)
            next_btn.set_visibility(i < len(STAGES) - 1)
        with body:
            if S.project is None:
                ui.label("REGULAR WORKFLOW").classes("edu-kicker mb-2")
                ui.label(title).style(
                    "font-size:1.9rem; font-weight:700; color:#17344b;letter-spacing:-.5px")
            else:
                # A coloured banner carrying the step number, its
                # picture and its one-line purpose. The number answers
                # "where am I", the colour matches the sidebar, and the
                # blurb answers "why am I here" before any explanation
                # is read.
                with ui.row().classes(
                    "items-center no-wrap w-full q-pa-md q-mb-md"
                ).style(
                    f"background:linear-gradient(90deg,{accent}26,"
                    f"{accent}08); border-left:7px solid {accent}; "
                    "border-radius:10px"
                ):
                    with ui.element("div").style(
                        f"background:{accent}; border-radius:50%; "
                        "width:52px; height:52px; display:flex; "
                        "align-items:center; justify-content:center; "
                        "flex:0 0 52px"
                    ):
                        ui.icon(icon_name).style(
                            "color:#ffffff; font-size:1.7rem")
                    with ui.column().classes("q-gutter-none").style(
                            "flex:1"):
                        ui.label(
                            f"STEP {i + 1} OF {len(STAGES)}"
                        ).style(
                            f"color:{accent}; font-weight:700; "
                            "font-size:0.72rem; letter-spacing:0.1em")
                        ui.label(title).style(
                            "font-size:1.55rem; font-weight:600; "
                            "color:#202124; line-height:1.2")
                        ui.label(blurb).style(
                            "font-size:0.95rem; color:#49505a")

            if S.project is None:
                ui.label(
                    "Investigate your own protein and molecules through eight steps. "
                    "Open a saved project or create a new one. Your structures, settings "
                    "and results stay together in its project folder."
                ).classes("q-mt-sm edu-muted").style("line-height:1.7;font-size:16px")
                with ui.row().classes("edu-info w-full items-center gap-3 mt-5 mb-3"):
                    ui.icon("school",size="28px",color="primary")
                    with ui.column().classes("gap-0").style("flex:1;min-width:200px"):
                        ui.label("Run a docking study with simple steps.").classes("font-semibold")
                        ui.label("Multiple protein targets, multiple drugs, 3D views and downloadable results.").classes("text-sm edu-muted")
                    ui.button("Open guided docking",on_click=lambda:ui.navigate.to("/study")).props("unelevated icon=arrow_forward no-caps color=primary")

                existing = list_projects()
                if existing:
                    ui.label("YOUR PROJECTS").classes("q-mt-md").style(
                        "font-size:0.78rem; letter-spacing:0.08em; "
                        "color:#5f6368")
                    for item in existing:
                        def open_it(folder=item["folder"],
                                    nm=item["name"]) -> None:
                            S.project = Project.open(folder)
                            project_label.text = nm
                            ui.notify(f"Opened {nm}", type="positive")
                            refresh_sidebar()
                            render()

                        with ui.card().classes(
                            "w-full q-mt-xs cursor-pointer"
                        ).style("box-shadow:none; border:1px solid #d8dee4"
                                ) as card:
                            with ui.row().classes("items-center no-wrap"):
                                ui.icon("folder_open").style("color:#2e5496")
                                with ui.column().classes("q-gutter-none"):
                                    ui.label(item["name"]).style(
                                        "font-weight:500")
                                    ui.label(
                                        f"{item['steps']} of {len(STAGES)} steps done "
                                        f"- last used "
                                        f"{item['when']:%d %b %Y, %H:%M}"
                                    ).style("font-size:0.82rem; "
                                            "color:#5f6368")
                                ui.space()
                                ui.icon("chevron_right").style(
                                    "color:#9aa0a6")
                            card.on("click", guard(open_it, "opening that project"))

                # Where the work is saved, and how to move it. On a
                # departmental Windows machine the home folder is on a
                # small C: drive, and the person who needs to move this
                # is not going to edit a Python file to do it.
                trouble = problem_with_root()
                if trouble:
                    with ui.card().classes("w-full q-mt-md").style(
                        "background:#fff8e1; border-left:6px solid "
                        "#c47f00; box-shadow:none"
                    ):
                        ui.label(trouble).style("line-height:1.55")

                with ui.expansion(
                    f"Saving in: {projects_root()}"
                ).classes("w-full q-mt-md").props("dense icon=folder"):
                    ui.label(
                        "Everything a project needs is inside its own "
                        "folder here - the structures, the poses, the "
                        "settings used and a copy of everything the "
                        "websites said. You can copy a folder to "
                        "another computer and open it there."
                    ).style("line-height:1.6; font-size:0.9rem")
                    ui.label(
                        "If this drive is short of space, put them "
                        "somewhere else. Type the full folder path, "
                        "for example  E:\\InSilicoProjects"
                    ).classes("q-mt-sm").style(
                        "line-height:1.6; font-size:0.9rem")
                    with ui.row().classes("w-full no-wrap items-center "
                                          "q-mt-xs q-gutter-sm"):
                        where = ui.input(
                            "Folder for new projects",
                            value=str(projects_root()),
                        ).props("outlined dense").style("flex:1")

                        def move_root() -> None:
                            ok, said = set_projects_root(where.value)
                            ui.notify(said,
                                      type="positive" if ok else "warning",
                                      timeout=12000, multi_line=True)
                            if ok:
                                render()

                        ui.button("Use this folder",
                                  on_click=guard(move_root,
                                                 "changing the folder")
                                  ).props("unelevated color=primary")
                    ui.label(
                        f"Leave it empty and press the button to go "
                        f"back to {default_root()}."
                    ).classes("q-mt-xs").style(
                        "font-size:0.82rem; color:#5f6368")

                ui.label("START SOMETHING NEW").classes("q-mt-lg").style(
                    "font-size:0.78rem; letter-spacing:0.08em; color:#5f6368")
                with ui.row().classes("q-mt-xs items-center"):
                    name = ui.input("Project name",
                                    value="My first docking").props("outlined")

                    def create() -> None:
                        folder = projects_root() / name.value
                        if folder.exists():
                            ui.notify(
                                "A project with that name already exists - "
                                "open it above, or choose another name.",
                                type="warning", timeout=7000)
                            return
                        S.project = Project.create(folder, name=name.value)
                        project_label.text = name.value
                        ui.notify(f"Created in {folder}", type="positive")
                        refresh_sidebar()
                        render()

                    ui.button("Create project",
                              on_click=guard(create,
                                             "creating the project")).props(
                        "unelevated color=primary")
                return

            teaching, control, summary = load_stage(S.current)
            old_step = S.project.get_step(S.current)
            if old_step and old_step["status"] == "stale":
                ui.label("Needs updating: inputs changed after this result was calculated. Run this step again. Older files are retained in the project folder.").classes("q-pa-md w-full").style(
                    "background:#fff4d7;border-left:5px solid #b26a00;color:#795200")
            teaching_panel(teaching, accent)

            if S.current in {"02_structure", "03_receptor", "04_site", "05_ligands", "06_dock"}:
                with ui.expansion("Explore the current structure in 3D", value=True).classes("w-full q-mt-md"):
                    from core.viewer import render as render_view
                    render_view(S.project, S.current)

            ui.separator().classes("q-my-md")

            # --- inputs and the button
            section_heading("SET IT UP AND RUN IT", "tune", accent)
            form = StageForm(S.current, S.project)
            with ui.column().classes("w-full q-gutter-sm q-mt-xs"):
                form.build()

            with ui.row().classes("items-center q-mt-md q-gutter-md"):
                run_button = ui.button(
                    f"Run step {i + 1}: {title}",
                ).props("unelevated size=lg icon=play_arrow").style(
                    f"background:{accent} !important; color:#ffffff")
                run_button.on("click",
                              lambda: do_run(S.current, form, run_button))
                ui.label().bind_text_from(progress, "text").style(
                    "color:#5f6368; font-size:0.9rem")

            # --- result and its control
            if summary:
                ui.separator().classes("q-my-md")
                section_heading("WHAT CAME OUT", "fact_check", accent)
                ui.label(summary).classes("q-mt-xs").style(
                    "font-size:1.05rem; line-height:1.5")
                control_panel(control)
                if S.current == "06_dock":
                    docking_check_panels()
                download_panel(S.current, accent)
                notebook_panel(S.current, accent)
                detail_panel(S.current)
            station_panel(S.current)

    S.render = render
    project_label.text = folder.name
    render()


def station_panel(stage_id: str) -> None:
    """
    The handoff to an outside website, and the way back.

    SwissADME has no programming interface and its terms forbid
    automated access, so this software never contacts it. It prepares
    the input, opens the site in the user's own browser with the
    molecules already loaded, and reads the file they download back
    into the project.

    The return path is the part that matters. Without it the user ends
    up with a browser tab full of numbers and no idea which of them
    bear on their question - which is what every list of useful links
    leaves them with.
    """
    if stage_id != "08_admet" or S.project is None:
        return
    if not S.project.is_complete("05_ligands"):
        return

    from stations.swissadme import SwissADMEStation

    station = SwissADMEStation(S.project)
    spec = station.spec()

    ui.separator().classes("q-my-md")
    ui.label("GET MORE, FROM SWISSADME").style(
        "font-size:0.78rem; letter-spacing:0.08em; color:#5f6368")

    with ui.card().classes("w-full q-mt-xs").style(
        "box-shadow:none; border:1px solid #d8dee4"
    ):
        ui.label(
            "The properties above are calculated here. SwissADME adds "
            "predictions this software does not make - gut absorption, "
            "whether the molecule is likely to reach the brain, which "
            "liver enzymes it may block, and a picture called the "
            "BOILED-Egg that shows both at a glance."
        ).style("line-height:1.6")
        ui.label(
            "It is a free academic service and its terms do not allow "
            "software to contact it automatically, so you open it "
            "yourself. Your molecules are loaded in for you - press the "
            "button, then Run on their page."
        ).classes("q-mt-sm").style("font-size:0.88rem; color:#5f6368")

        def open_station() -> None:
            # Prepare it first, so a failure is reported here rather
            # than as a blank tab the user cannot interpret.
            station.build_launcher()
            from urllib.parse import quote
            ui.navigate.to("/station/swissadme?project=" + quote(S.project.folder.name), new_tab=True)
            ui.notify(
                "Opened in a new tab. Press the button there, then Run "
                "on the SwissADME page.",
                type="positive", timeout=8000, multi_line=True)

        with ui.row().classes("q-mt-md items-center q-gutter-md"):
            ui.button("Open SwissADME with my molecules",
                      on_click=guard(open_station, "opening SwissADME")
                      ).props("unelevated color=primary icon=open_in_new")

        with ui.expansion("What to do on their page").classes(
            "w-full q-mt-sm"
        ):
            for i, line in enumerate(spec.instructions, start=1):
                ui.label(f"{i}. {line}").style("line-height:1.6")
            ui.label(spec.terms_note).classes("q-mt-sm").style(
                "font-size:0.85rem; color:#5f6368")

        ui.separator().classes("q-my-md")
        ui.label("Bring the results back").style("font-weight:600")
        ui.label(
            spec.expect_back + " Drop it here and it will be read in, "
            "checked, and added to your report and your record of "
            "methods."
        ).style("font-size:0.9rem; color:#5f6368; line-height:1.55")

        def ingest(event: Any) -> None:
            incoming = S.project.folder / "incoming"
            incoming.mkdir(exist_ok=True)
            path = incoming / event.name
            path.write_bytes(event.content.read())
            result = station.ingest(path)
            if not result.success:
                ui.notify(result.message, type="negative", timeout=10000,
                          multi_line=True)
                return
            S.project.save_step("08_admet_swissadme",
                                "SwissADME results", result)
            ui.notify(result.summary, type="positive")
            S.render()

        ui.upload(on_upload=ingest, auto_upload=True,
                  label="Drop swissadme.csv here").props(
            "accept=.csv flat bordered").classes("w-full q-mt-sm")

    saved = S.project.get_step("08_admet_swissadme")
    if saved and saved["status"] == "complete":
        ui.label("FROM SWISSADME").classes("q-mt-md").style(
            "font-size:0.78rem; letter-spacing:0.08em; color:#5f6368")
        ui.label(saved["summary"]).style("line-height:1.5")
        control = rehydrate_control(saved["control_json"])
        control_panel(control)
        values = json.loads(saved["values_json"] or "{}")
        for m in values.get("molecules", []):
            with ui.expansion(m.get("name") or "molecule").classes("w-full"):
                for key, value in m.items():
                    if key in ("name", "_all"):
                        continue
                    with ui.row().classes("no-wrap"):
                        ui.label(f"{key}:").style(
                            "color:#5f6368; min-width:230px")
                        ui.label(str(value))



def notebook_panel(stage_id: str, accent: str) -> None:
    """
    The same docking, as Python the student can run themselves.

    The point of handing this over is not convenience. A student who
    only ever presses buttons learns a workflow; one who can read the
    forty lines that place the box and the twenty that measure the
    re-docking RMSD learns what the workflow is for. It is also the
    only way their work outlives this application.
    """
    if S.project is None:
        return
    step = S.project.get_step(stage_id)
    if not step or step["status"] != "complete" or not step["values_json"]:
        return
    notebook = json.loads(step["values_json"]).get("colab_notebook", "")
    if not notebook or not Path(notebook).exists():
        return

    with ui.card().classes("w-full q-mt-sm").style(
        f"border:2px dashed {accent}; background:{accent}0a; "
        "box-shadow:none"
    ):
        with ui.row().classes("items-center no-wrap w-full"):
            ui.icon("code").style(f"color:{accent}; font-size:2rem")
            with ui.column().classes("q-gutter-none").style("flex:1"):
                ui.label("Run this yourself, in Python").style(
                    "font-weight:600; font-size:1.05rem")
                ui.label(
                    "A Google Colab notebook that does everything on "
                    "this screen - the same receptor conversion, the "
                    "same decoy library, the same three checks. Change "
                    "one cell to use your own protein and molecules. "
                    "Free, no installation, runs in a browser."
                ).style("font-size:0.88rem; color:#5f6368")
            ui.button(
                "Download notebook",
                on_click=lambda p=notebook: ui.download(Path(p)),
            ).props("unelevated icon=download").style(
                f"background:{accent} !important; color:#ffffff")
        with ui.expansion("How a student uses it").classes(
                "w-full").props("dense"):
            for line in (
                "1. Download it, then go to colab.research.google.com "
                "and sign in with any Google account.",
                "2. File -> Upload notebook -> choose the file.",
                "3. Run the cells from the top. The first one installs "
                "everything and takes about two minutes.",
                "4. Cell 3 is the only one to edit: it holds the PDB "
                "code and the SMILES of the molecules.",
            ):
                ui.label(line).classes("q-mb-xs").style(
                    "font-size:0.88rem; line-height:1.5")


def download_panel(stage_id: str, accent: str) -> None:
    """
    The way a job package actually reaches the user.

    An EXPORT stage writes a folder of files, and a folder inside the
    project directory is of no use to somebody who does not know where
    the project directory is. This hands them the zip.
    """
    if S.project is None:
        return
    step = S.project.get_step(stage_id)
    if not step or step["status"] != "complete" or not step["values_json"]:
        return
    values = json.loads(step["values_json"])
    archive = values.get("zip", "")
    if not archive or not Path(archive).exists():
        return

    size_mb = Path(archive).stat().st_size / 1_000_000
    with ui.card().classes("w-full q-mt-sm").style(
        f"border:2px dashed {accent}; background:{accent}0a; "
        "box-shadow:none"
    ):
        with ui.row().classes("items-center no-wrap w-full"):
            ui.icon("cloud_download").style(
                f"color:{accent}; font-size:2rem")
            with ui.column().classes("q-gutter-none").style("flex:1"):
                ui.label("Your job package is ready").style(
                    "font-weight:600; font-size:1.05rem")
                ui.label(
                    f"{Path(archive).name} - {size_mb:.1f} MB. Unzip it, "
                    "open README.txt, and follow the Google Colab route."
                ).style("font-size:0.88rem; color:#5f6368")
            ui.button(
                "Download",
                on_click=lambda a=archive: ui.download(Path(a)),
            ).props("unelevated icon=download").style(
                f"background:{accent} !important; color:#ffffff")
        ui.label(f"It is also saved at: {archive}").style(
            "font-size:0.78rem; color:#9aa0a6")


def detail_panel(stage_id: str) -> None:
    """The numbers from this stage, under the control, collapsed."""
    if S.project is None:
        return
    step = S.project.get_step(stage_id)
    if not step or step["status"] != "complete" or not step["values_json"]:
        return
    import json as _json
    values = _json.loads(step["values_json"])

    rows = _summarise(stage_id, values)
    if not rows:
        return
    with ui.expansion("The details").classes("w-full q-mt-sm"):
        for key, value in rows:
            with ui.row().classes("no-wrap items-start q-mb-xs"):
                ui.label(f"{key}:").style(
                    "color:#5f6368; min-width:210px")
                ui.label(str(value)).style("line-height:1.5")


def _summarise(stage_id: str, v: dict) -> list[tuple[str, Any]]:
    """What is worth showing from each stage, in plain words."""
    if stage_id == "01_target":
        t = v.get("target", {})
        h = v.get("history", {})
        return [
            ("Protein", t.get("protein", "")),
            ("Gene", ", ".join(t.get("genes", []))),
            ("UniProt", t.get("accession", "")),
            ("Length", f"{t.get('length', 0)} residues"),
            ("Experimental structures", len(t.get("pdb_ids", []))),
            ("Evidence it can be drugged",
             ", ".join(h.get("sm_buckets", [])) or "none found"),
        ]
    if stage_id == "02_structure":
        c = v.get("chosen", {})
        rows = [(k.replace("_", " ").capitalize(), c.get(k))
                for k in ("id", "method", "resolution", "released",
                          "n_ligands", "residues_missing", "mean_plddt")
                if c.get(k) is not None]
        if c.get("note"):
            rows.append(("Note", c["note"]))
        # The alternatives were ranked and then hidden, which leaves
        # "how do I choose the best structure?" unanswerable from
        # inside the software. They are the answer: the shortlist, the
        # numbers that decided it, and the reason the winner won.
        alternatives = v.get("alternatives") or []
        if alternatives:
            rows.append(("", ""))
            rows.append(("Also considered",
                         "ranked below the one chosen - to use one "
                         "instead, type its code in the box above and "
                         "run again"))
            for alt in alternatives[:4]:
                bits = []
                if alt.get("resolution"):
                    bits.append(f"{alt['resolution']} A")
                if alt.get("method"):
                    bits.append(str(alt["method"]).split()[0].lower())
                bits.append("drug bound" if alt.get("n_ligands")
                            else "nothing bound")
                if alt.get("residues_missing") is not None:
                    bits.append(f"{alt['residues_missing']} missing")
                rows.append((f"   {alt.get('id', '?')}", ", ".join(bits)))
        if c.get("chains_dropped"):
            rows.append(("Chains left out",
                         ", ".join(c["chains_dropped"])
                         + " - their labels are too long for the file "
                           "format; check none of them is your protein"))
        return rows
    if stage_id == "03_receptor":
        return [
            ("Ligand found", v.get("ligand", "none")),
            ("Kept", ", ".join(k["label"] for k in v.get("kept", []))),
            ("Removed", ", ".join(k["label"] for k in v.get("removed", []))
             or "nothing"),
            ("Waters removed", v.get("waters_removed", 0)),
            ("Hydrogens added", v.get("hydrogens_added", 0)),
            ("Gaps left open", len(v.get("gaps_left_open", []))),
        ]
    if stage_id == "04_site":
        b = v.get("box", {})
        sh = v.get("shape", {})
        return [
            ("Box", f"{b.get('size_x')} x {b.get('size_y')} x "
                    f"{b.get('size_z')} A"),
            ("Volume", f"{b.get('volume_A3')} A3"),
            ("Chosen from", b.get("source", "")),
            ("Space a ligand could occupy",
             f"{sh.get('usable_fraction', 0):.0%}"),
            ("Enclosed by protein", f"{sh.get('buriedness', 0):.0%}"),
        ]
    if stage_id == "05_ligands":
        rows = []
        for m in v.get("prepared", []):
            changes = "; ".join(m.get("changes", [])) or "unchanged"
            rows.append((m.get("name") or "molecule", changes))
        for m in v.get("failed", []):
            rows.append((m.get("name") or "molecule",
                         f"could not be read: {m.get('error', '')}"))
        return rows
    if stage_id == "08_admet":
        rows = []
        for m in v.get("molecules", []):
            pr = m["properties"]
            ru = m["rules"]
            rows.append((m["name"],
                         f"MW {pr['mw']}, logP {pr['logp']}, "
                         f"polar surface {pr['tpsa']}, "
                         f"{pr['rot']} rotatable bonds"))
            lip = ru["lipinski_detail"]
            rows.append((f"  {m['name']} - Lipinski",
                         "; ".join(lip) if lip
                         else "no violations (a description, not a pass mark)"))
        return rows
    if stage_id == "09_md":
        job = v.get("job", {})
        size = job.get("system", {})
        rows = [
            ("Molecule", job.get("molecule", "")),
            ("Simulation length", f"{job.get('nanoseconds', 0):g} ns"),
            ("Atoms once in water", f"{size.get('total_atoms', 0):,}"),
            ("Water box", f"{size.get('box_nm', 0)} nm per side"),
        ]
        for r in v.get("routes", []):
            cost = (f"about ${r['cost_usd']:.2f}" if r["cost_usd"]
                    else "free")
            rows.append((f"  {r['route']}", f"{r['human']} - {cost}"))
        return rows
    if stage_id == "06_dock":
        rows = []
        # The three checks, named, before the scores. The RMSD in
        # particular has to be a field of its own: it is the number
        # that says whether any of the scores below mean anything,
        # and it is the one people are told to report.
        source = v.get("box_source", "")
        if source and "detected" in source:
            rows.append(("Where the box came from",
                         f"{source} - NOT from a bound drug, so it may "
                         "not be the real binding site"))
        elif source:
            rows.append(("Where the box came from", source))
        rmsd = v.get("redock_rmsd")
        if rmsd is not None:
            rows.append(("Re-docking RMSD",
                         f"{rmsd} A against the crystal pose of "
                         f"{v.get('redock_molecule', 'the bound ligand')}"
                         + ("  - passes the 2 A reference-pose threshold" if rmsd <= 2.0
                            else "  - does not pass the 2 A reference-pose threshold")))
        elif v.get("redock_note"):
            rows.append(("Re-docking RMSD",
                         f"could not be measured: {v['redock_note']}"))
        confound = v.get("size_confound") or {}
        if confound.get("usable"):
            rows.append(("Score-size association (R-squared)",
                         f"{confound['r_squared']:.2f}  "
                         f"(spread {confound['spread']} kcal/mol; "
                         "association does not establish its cause)"))
        if rows:
            rows.append(("", ""))
        for r in v.get("results", []):
            if "score" in r:
                rows.append((r["name"] or "molecule",
                             f"{r['score']} kcal/mol"))
            else:
                rows.append((r.get("name", "molecule"),
                             f"failed: {r.get('error', '')}"))
        return rows
    return []


def docking_check_panels() -> None:
    """Present every docking check without changing the recorded result."""
    if S.project is None:
        return
    step = S.project.get_step("06_dock")
    if not step or step["status"] != "complete":
        return
    from core.workshop import docking_checks

    recorded = json.loads(step["values_json"] or "{}")
    original = json.loads(step["control_json"] or "{}")
    checks = docking_checks(recorded, original)
    if not checks:
        return
    ui.label("INDIVIDUAL DOCKING CHECKS").classes("q-mt-md").style(
        "font-size:0.8rem; font-weight:700; color:#5f6368")
    for name, outcome in checks.items():
        ui.label(name.capitalize() + " check").classes("q-mt-sm").style(
            "font-weight:600; color:#1f3864")
        control_panel(rehydrate_control(json.dumps(outcome)))


def load_stage(stage_id: str) -> tuple[Any, Any, str]:
    """Teaching content always; result and control if the step has run."""
    teaching = control = None
    summary = ""

    try:
        stage = build_stage(stage_id)
        if stage is not None:
            teaching = stage.explain()
    except Exception:  # noqa: BLE001
        teaching = None

    if S.project is not None:
        step = S.project.get_step(stage_id)
        if step and step["status"] == "complete":
            summary = step["summary"] or ""
            blob = step["control_json"]
            from core.workshop import teaching_check
            blob = json.dumps(teaching_check(
                stage_id, json.loads(blob or "{}"),
                recorded=json.loads(step["values_json"] or "{}")))
            control = rehydrate_control(blob) if blob and blob != "{}" else None
    return teaching, control, summary


def build_stage(stage_id: str) -> Any:
    """A stage instance, only for its teaching text."""
    if S.project is None:
        return None
    from stages.ligands import LigandStage
    from stages.receptor import ReceptorStage
    from stages.site import SiteStage
    from stages.structure import StructureStage
    from stages.target import TargetStage

    from stages.admet import AdmetStage
    from stages.dock import DockStage
    from stages.md import MdStage

    table = {
        "01_target": TargetStage, "02_structure": StructureStage,
        "03_receptor": ReceptorStage, "04_site": SiteStage,
        "05_ligands": LigandStage, "06_dock": DockStage,
        "09_md": MdStage, "08_admet": AdmetStage,
    }
    cls = table.get(stage_id)
    return cls(S.project) if cls else None


def rehydrate_control(blob: str | None) -> Any:
    if not blob:
        return None
    import json
    from types import SimpleNamespace

    data = json.loads(blob)
    data["verdict"] = SimpleNamespace(value=data.get("verdict", "not_applicable"))
    return SimpleNamespace(**data)


from workshop_ui import register_workshop
register_workshop()

from study_ui import register_study
register_study()


if __name__ in {"__main__", "__mp_main__"}:
    import os
    ui.run(title="In-Silico Docking Studio", reload=False,
           host=os.getenv("INSILICO_HOST", "127.0.0.1"),
           port=int(os.getenv("INSILICO_PORT", "8111")), show=False, favicon="🧬")

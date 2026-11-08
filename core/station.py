"""
Station: the pattern for a step this software cannot run itself.

Some of the best free tools in this field cannot lawfully or practically
be automated. SwissADME's terms of use prohibit accessing it with any
crawler or automated data-retrieval tool, and it offers no programming
interface; pkCSM likewise has none. A person clicking a button and
reading the result is ordinary permitted use.

So a station does four things -- specification Section 6.2:

    WHY    what this step tells you, and what it cannot tell you
    GO     one button, with the input already prepared
    THERE  numbered instructions and a picture of the real page
    BACK   a drop-zone that parses the file the user downloaded

The BACK zone is the whole difference between this software and a page
of bookmarks.

Rule enforced in code below: the launcher page NEVER submits by itself.
A human presses the button. Automatic submission would cross back into
the automated access that the terms prohibit.
"""

from __future__ import annotations

import html
import webbrowser
from abc import abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import urlencode

from .contract import Stage, StageResult, Tier


@dataclass
class StationSpec:
    """
    Everything about one external service, held as data.

    External sites are redesigned every few years and screenshots rot.
    Keeping this as data -- loaded from an editable file shipped
    separately from the program -- means a broken station can be fixed
    and redistributed without a new software release.
    Specification Section 6.3, maintenance note.
    """

    key: str                     # e.g. "swissadme"
    name: str                    # e.g. "SwissADME"
    url: str                     # where the form posts, or the page to open
    method: str = "post"         # "get" (deep link) or "post" (launcher)
    field_name: str = "smiles"   # the form field the input goes into

    why: str = ""                # WHY zone
    cannot_tell_you: str = ""    # the honest half of the WHY zone
    instructions: list[str] = field(default_factory=list)   # THERE zone
    screenshot: str = ""         # file name inside the station assets folder
    expect_back: str = ""        # what file the user should return with
    citation: str = ""
    terms_note: str = ""         # shown to the user, e.g. academic use only


# ----------------------------------------------------------------------
# The launcher page
# ----------------------------------------------------------------------

_LAUNCHER_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{name} - step ready</title>
<style>
  :root {{ color-scheme: light dark; }}
  body {{ font-family: system-ui, "Segoe UI", sans-serif; max-width: 720px;
         margin: 0 auto; padding: 40px 20px; line-height: 1.6; }}
  h1 {{ font-size: 1.4rem; margin-bottom: .2em; }}
  .why {{ background: #eef3f8; border-left: 4px solid #2e5496;
          padding: 14px 18px; margin: 20px 0; border-radius: 4px; }}
  .cannot {{ color: #7a3b00; background: #fff6e8; border-left: 4px solid #c47f00;
             padding: 14px 18px; margin: 20px 0; border-radius: 4px; }}
  textarea {{ width: 100%; font-family: ui-monospace, monospace; font-size: .95rem;
              padding: 10px; box-sizing: border-box; }}
  button {{ font-size: 1.15rem; padding: 14px 28px; cursor: pointer;
            background: #2e5496; color: #fff; border: 0; border-radius: 6px; }}
  ol {{ padding-left: 1.3em; }}
  .back {{ margin-top: 30px; padding-top: 18px; border-top: 1px solid #ccd4dd;
           color: #444; }}
  .terms {{ font-size: .85rem; color: #666; margin-top: 28px; }}
  @media (prefers-color-scheme: dark) {{
    body {{ background: #15181c; color: #e6e8ea; }}
    .why {{ background: #1d2733; }}
    .cannot {{ background: #2b2213; color: #ffd59b; }}
  }}
</style>
</head>
<body>

<h1>{name}</h1>

<div class="why"><strong>What this step tells you.</strong> {why}</div>
<div class="cannot"><strong>What it cannot tell you.</strong> {cannot}</div>

<form action="{url}" method="{method}" target="_blank">
  <label for="payload"><strong>Your input is already prepared:</strong></label>
  <textarea id="payload" name="{field_name}" rows="4">{payload}</textarea>
  <p><button type="submit">Open {name} &rarr;</button></p>
</form>

<h2>Once you are there</h2>
<ol>
{instructions}
</ol>

{screenshot_block}

<div class="back">
  <strong>Then come back.</strong> {expect_back}
  Return to the software and drop that file onto the step - it will be
  read in, checked, and added to your report and your record of methods.
</div>

<p class="terms">{terms}{citation}</p>

</body>
</html>
"""


class StationStage(Stage):
    """
    A stage that hands the user to an external website.

    Subclasses provide spec() and ingest(). run() prepares the input,
    writes the launcher page and opens it in the user's own browser.
    """

    tier = Tier.STATION

    @abstractmethod
    def spec(self) -> StationSpec:
        """The service this station hands off to."""

    @abstractmethod
    def payload(self) -> str:
        """The prepared input, e.g. a list of SMILES, one per line."""

    # -- building the handoff -------------------------------------------

    def build_launcher(self) -> Path:
        """Write the local launcher page and return its path."""
        spec = self.spec()
        payload = self.payload()

        steps = "\n".join(
            f"  <li>{html.escape(s)}</li>" for s in spec.instructions
        ) or "  <li>Follow the instructions on the page.</li>"

        shot = ""
        if spec.screenshot:
            asset = self.project.step_folder(self.stage_id) / spec.screenshot
            if asset.exists():
                shot = (
                    f'<p><img src="{html.escape(asset.name)}" '
                    f'alt="What the {html.escape(spec.name)} page looks like" '
                    'style="max-width:100%;border:1px solid #ccd4dd;border-radius:4px"></p>'
                )

        terms = f"{spec.terms_note} " if spec.terms_note else ""

        page = _LAUNCHER_TEMPLATE.format(
            name=html.escape(spec.name),
            why=html.escape(spec.why),
            cannot=html.escape(spec.cannot_tell_you),
            url=html.escape(spec.url),
            method=html.escape(spec.method),
            field_name=html.escape(spec.field_name),
            payload=html.escape(payload),
            instructions=steps,
            screenshot_block=shot,
            expect_back=html.escape(spec.expect_back),
            terms=html.escape(terms),
            citation=html.escape(spec.citation),
        )

        # Guard rail, checked on every build: the page must not submit
        # itself. See the module docstring.
        lowered = page.lower()
        for forbidden in ("submit()", "onload=", "setinterval", "settimeout"):
            if forbidden in lowered:
                raise RuntimeError(
                    "The launcher page must never submit by itself. "
                    f"Found {forbidden!r} in the generated page."
                )

        path = self.project.step_folder(self.stage_id) / f"open_{spec.key}.html"
        path.write_text(page, encoding="utf-8")
        return path

    def run(self, open_browser: bool = True, **kwargs) -> StageResult:
        spec = self.spec()
        path = self.build_launcher()
        if open_browser:
            webbrowser.open(path.as_uri())
        return StageResult(
            success=True,
            summary=f"{spec.name} opened with your input prepared.",
            values={"station": spec.key, "awaiting": spec.expect_back},
            artefacts=[path],
        )

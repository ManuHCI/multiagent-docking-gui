"""
Where this software keeps things, and how that gets changed.

One setting so far, and it matters more than one setting usually
would: which folder the projects go in.

The default is a folder in the user's home directory, which on
Windows means the C: drive. That is the right default -- it always
exists and always writable -- and it is wrong for a good number of
the people this software is for. A departmental machine often has a
small system drive and a large second one, and a docking project
carries downloaded structures, prepared receptors, every docked pose
and a cache of every answer the websites gave. It is not enormous,
but it accumulates, and the person who most needs the space is
usually the one least able to move a folder afterwards.

So: an environment variable for anyone scripting this, a small file
for the setting made through the interface, and the home directory
when neither says otherwise. The file lives beside app.py rather than
in the user's profile so that copying the whole folder to another
machine carries the setting with it.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

APP_DIR = Path(__file__).resolve().parent.parent
from .runtime import data_root

PACKAGED = bool(getattr(sys, "frozen", False))
SETTINGS_FILE = (data_root() if PACKAGED or os.environ.get("INSILICO_DATA_DIR")
                 else APP_DIR) / "settings.json"
ENV_VAR = "INSILICO_PROJECTS"
DEFAULT_NAME = "InSilicoProjects"


def _read() -> dict:
    try:
        return json.loads(SETTINGS_FILE.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def _write(data: dict) -> bool:
    try:
        SETTINGS_FILE.parent.mkdir(parents=True, exist_ok=True)
        SETTINGS_FILE.write_text(json.dumps(data, indent=2),
                                 encoding="utf-8")
        return True
    except OSError:
        # An installation folder the user cannot write to. Rare, but
        # it should not take the software down with it: the setting is
        # a convenience and the default always works.
        return False


def configured_root() -> Path | None:
    """What has been asked for, before checking whether it works."""
    from_env = os.environ.get(ENV_VAR, "").strip()
    if from_env:
        return Path(from_env)
    saved = _read().get("projects_root", "")
    return Path(saved) if saved else None


def default_root() -> Path:
    if PACKAGED or os.environ.get("INSILICO_DATA_DIR"):
        return data_root() / "Projects"
    return Path.home() / DEFAULT_NAME


def projects_root() -> Path:
    """
    The folder to use now.

    Falls back to the home directory if the configured one cannot be
    made -- a drive that is not plugged in, a path on a share that is
    not mounted today. Silently losing the setting would be wrong, so
    `problem_with_root()` says so when something is asked.
    """
    wanted = configured_root()
    if wanted is not None:
        try:
            wanted.mkdir(parents=True, exist_ok=True)
            probe = wanted / ".writable"
            probe.touch()
            probe.unlink()
            return wanted
        except OSError:
            pass
    root = default_root()
    root.mkdir(parents=True, exist_ok=True)
    return root


def problem_with_root() -> str:
    """A sentence to show when the chosen folder is not usable."""
    wanted = configured_root()
    if wanted is None:
        return ""
    if projects_root() == wanted:
        return ""
    return (
        f"Projects were set to be saved in {wanted}, but that folder "
        "cannot be written to right now - the drive may not be "
        f"connected. Using {default_root()} instead for now. The "
        "setting has not been changed."
    )


def set_projects_root(folder: str | Path) -> tuple[bool, str]:
    """
    Change where projects are saved.

    Returns whether it worked and a sentence for the user. Existing
    projects are deliberately not moved: moving a user's data without
    being asked is how data gets lost, and the old folder is named so
    they can move it themselves.
    """
    text = str(folder).strip().strip('"')
    if not text:
        return False, "No folder was given."
    folder = Path(text)

    # A drive letter that is not there. Windows raises on mkdir, but
    # on any other system "E:/Projects" is read as a relative path and
    # a folder called "E:" gets quietly created in the current
    # directory -- which is not what anybody meant.
    if len(text) > 1 and text[1] == ":" and os.name != "nt":
        return False, (
            f"{text} is a Windows drive path, and this is not "
            "Windows. Give a path for this system instead."
        )

    previous = projects_root()
    try:
        folder.mkdir(parents=True, exist_ok=True)
        probe = folder / ".writable"
        probe.touch()
        probe.unlink()
    except OSError as exc:
        return False, (
            f"That folder cannot be used: {exc.strerror or exc}. Check "
            "the drive letter and that the drive is connected."
        )

    data = _read()
    data["projects_root"] = str(folder)
    if not _write(data):
        return False, (
            "The setting could not be saved - this folder is "
            "read-only. You can set the "
            f"{ENV_VAR} environment variable instead."
        )

    if previous == folder:
        return True, f"Projects are saved in {folder}."
    return True, (
        f"New projects will be saved in {folder}. Projects already in "
        f"{previous} stay where they are - move the folders across "
        "yourself if you want them together."
    )

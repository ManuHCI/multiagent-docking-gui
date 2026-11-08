"""
First-run self-test.

Checks every bundled tool and every external service, and reports in
plain English what is working -- before the user begins, not when a
step fails halfway through. Specification Section 11.

Design rules, all of which exist because of the audience:

  * Nothing here raises. A missing tool is a finding, not a crash.
  * Every finding says what it means and what to do about it.
  * The network part is optional. The application must open, teach and
    run the worked cases with no connection at all, so an offline
    machine gets "working offline", not a wall of red.
  * The whole thing finishes in a few seconds. Network checks run with
    a short timeout and in parallel.
"""

from __future__ import annotations

import importlib
import platform
import shutil
import subprocess
import sys

from core.install import how_to_install
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path


class State(str, Enum):
    OK = "ok"
    MISSING = "missing"      # needed, not there
    OPTIONAL = "optional"    # not there, but the software still works
    OFFLINE = "offline"      # a service that could not be reached


@dataclass
class Finding:
    name: str
    state: State
    detail: str = ""          # version, or what is wrong
    meaning: str = ""         # what this means for the user
    fix: str = ""             # what to do about it

    @property
    def blocking(self) -> bool:
        return self.state is State.MISSING


@dataclass
class Report:
    findings: list[Finding] = field(default_factory=list)
    machine: dict[str, str] = field(default_factory=dict)

    @property
    def can_start(self) -> bool:
        return not any(f.blocking for f in self.findings)

    @property
    def online(self) -> bool:
        return not any(f.state is State.OFFLINE for f in self.findings)

    def plain_summary(self) -> str:
        if not self.can_start:
            broken = [f.name for f in self.findings if f.blocking]
            return (
                "Something needed is missing: "
                + ", ".join(broken)
                + ". The software cannot run until this is fixed."
            )
        if not self.online:
            return (
                "Everything on this computer is working. Some websites could "
                "not be reached, so you are working offline - the worked "
                "examples and all the teaching will still run."
            )
        return "Everything is working. You can begin."


# ----------------------------------------------------------------------
# What must be present
# ----------------------------------------------------------------------

# (import name, friendly name, what it does, required?)
PYTHON_PACKAGES = [
    ("rdkit", "RDKit", "draws and prepares molecules", True),
    ("openmm", "OpenMM", "tidies up protein structures", True),
    ("pdbfixer", "PDBFixer", "repairs gaps in protein structures", True),
    ("meeko", "Meeko", "prepares docking input files", True),
    ("nicegui", "NiceGUI", "opens the teaching interface", True),
    ("numpy", "NumPy", "calculates pocket geometry", True),
    ("scipy", "SciPy", "compares molecular geometry", True),
]

# (executable, friendly name, what it does, version flag)
EXECUTABLES = [
    ("vina", "AutoDock Vina", "does the docking", "--version"),
]

# (label, url) -- reachability only, no data fetched
SERVICES = [
    ("UniProt", "https://rest.uniprot.org/uniprotkb/P00533.json"),
    ("Protein Data Bank", "https://data.rcsb.org/rest/v1/core/entry/1M17"),
    ("AlphaFold Database", "https://alphafold.ebi.ac.uk/api/prediction/P00533"),
    ("PubChem", "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/CanonicalSMILES/JSON"),
    ("SwissADME (opened by you, never automatically)", "https://www.swissadme.ch/"),
]


# ----------------------------------------------------------------------
# Checks
# ----------------------------------------------------------------------

def _check_package(mod: str, name: str, does: str, required: bool) -> Finding:
    try:
        m = importlib.import_module(mod)
    except Exception as exc:  # noqa: BLE001
        return Finding(
            name=name,
            state=State.MISSING if required else State.OPTIONAL,
            detail=type(exc).__name__,
            meaning=(
                f"{name} is not installed. It {does}."
                if required
                else f"{name} is not installed, so the step that {does} "
                     "will be unavailable. Everything else works."
            ),
            fix="Run:  " + how_to_install(mod),
        )
    version = getattr(m, "__version__", "") or "installed"
    return Finding(name=name, state=State.OK, detail=str(version),
                   meaning=f"Ready - it {does}.")


def _check_executable(exe: str, name: str, does: str, flag: str) -> Finding:
    if exe == "vina":
        from stages.dock import find_vina
        path = find_vina()
    else:
        path = shutil.which(exe)
    if not path:
        return Finding(
            name=name,
            state=State.MISSING,
            detail="not found",
            meaning=f"{name} is not installed. It {does}, so docking cannot run.",
            fix="Run:  " + how_to_install(exe),
        )
    try:
        out = subprocess.run(
            [path, flag], capture_output=True, text=True, timeout=20
        )
        detail = (out.stdout or out.stderr).strip().splitlines()[0][:70]
    except Exception:  # noqa: BLE001
        detail = "found, version unknown"
    return Finding(name=name, state=State.OK, detail=detail,
                   meaning=f"Ready - it {does}.")


def _check_service(label: str, url: str, timeout: float = 6.0) -> Finding:
    import urllib.error
    import urllib.request

    req = urllib.request.Request(url, method="GET",
                                 headers={"User-Agent": "InSilicoTeachingSuite"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            code = resp.status
    except urllib.error.HTTPError as e:
        code = e.code          # a 4xx still means the site answered
    except Exception as exc:   # noqa: BLE001
        return Finding(
            name=label,
            state=State.OFFLINE,
            detail=type(exc).__name__,
            meaning=f"Could not reach {label}.",
            fix="Check your internet connection, or carry on offline.",
        )
    return Finding(name=label, state=State.OK, detail=f"responded ({code})",
                   meaning="Reachable.")


def run(check_network: bool = True) -> Report:
    """Run every check. Never raises."""
    report = Report(machine={
        "python": sys.version.split()[0],
        "system": f"{platform.system()} {platform.release()}",
        "machine": platform.machine(),
        "folder": str(Path(sys.prefix)),
    })

    for mod, name, does, required in PYTHON_PACKAGES:
        report.findings.append(_check_package(mod, name, does, required))

    for exe, name, does, flag in EXECUTABLES:
        report.findings.append(_check_executable(exe, name, does, flag))

    root = Path(__file__).resolve().parent.parent
    for name, relative in [("Offline 3D viewer", "assets/vendor/3Dmol-min.js"),
                           ("Prepared EGFR workshop", "assets/workshop/egfr/project.db")]:
        present = (root / relative).is_file()
        report.findings.append(Finding(
            name=name, state=State.OK if present else State.OPTIONAL,
            detail="bundled" if present else "not found",
            meaning="Ready for offline teaching." if present else
                    "The regular calculations work, but this workshop feature is unavailable.",
            fix="Restore the full assets folder supplied with this app." if not present else ""))

    if check_network:
        with ThreadPoolExecutor(max_workers=len(SERVICES)) as pool:
            futures = [pool.submit(_check_service, lbl, url)
                       for lbl, url in SERVICES]
            report.findings.extend(f.result() for f in futures)

    return report


# ----------------------------------------------------------------------
# Printing -- used by the installer's final step and the Help menu
# ----------------------------------------------------------------------

_MARK = {
    State.OK: "  ok  ",
    State.MISSING: " MISSING ",
    State.OPTIONAL: " absent ",
    State.OFFLINE: " offline ",
}


def print_report(report: Report) -> None:
    print()
    print("In-Silico Drug Design Teaching Suite - system check")
    print("=" * 58)
    print(f"  {report.machine['system']} ({report.machine['machine']}), "
          f"Python {report.machine['python']}")
    print()

    for f in report.findings:
        print(f"[{_MARK[f.state]:^9}] {f.name:<46} {f.detail}")
        if f.state is not State.OK:
            print(f"{'':>13}{f.meaning}")
            if f.fix:
                print(f"{'':>13}{f.fix}")

    print()
    print(report.plain_summary())
    print()


if __name__ == "__main__":
    rep = run(check_network="--offline" not in sys.argv)
    print_report(rep)
    raise SystemExit(0 if rep.can_start else 1)

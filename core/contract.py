"""
Module contract for the In-Silico Drug Design Teaching Suite.

Every stage of the workflow -- whether it runs locally, calls a web API,
hands the user off to an external website, or exports a cloud job --
is one subclass of Stage implementing the same five methods.

This uniformity is what allows a new tool to be added in an afternoon
rather than a month. See specification Section 5.2.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any


# ----------------------------------------------------------------------
# Tiers -- specification Section 5.1
# ----------------------------------------------------------------------

class Tier(str, Enum):
    """Where a stage actually executes."""

    BUNDLED = "bundled"   # runs offline on the user's laptop
    API = "api"           # free public endpoint, called from inside
    STATION = "station"   # hand the user to a website; ingest the result
    EXPORT = "export"     # write a job package for cloud compute


# ----------------------------------------------------------------------
# Readiness -- what check() reports
# ----------------------------------------------------------------------

class Readiness(str, Enum):
    READY = "ready"             # go ahead
    BLOCKED = "blocked"         # cannot proceed; reason explains why
    NEEDS_INPUT = "needs_input"  # an earlier stage has not been completed
    DEGRADED = "degraded"       # will work, but something is reduced


@dataclass
class CheckResult:
    """
    The outcome of check(). Always carries plain English.

    A novice user must never see a stack trace. If something is wrong,
    `reason` says what in clinical English and `fix` says what to do
    about it. `technical` holds the real error for the log file only.
    """

    readiness: Readiness
    reason: str = ""
    fix: str = ""
    technical: str = ""

    @property
    def ok(self) -> bool:
        return self.readiness in (Readiness.READY, Readiness.DEGRADED)


# ----------------------------------------------------------------------
# Teaching content -- specification Section 7.1
# ----------------------------------------------------------------------

@dataclass
class Teaching:
    """
    Four depths, always available, never computed.

    Held as data rather than code so that the content can be reviewed
    and corrected by a colleague who does not read Python.

    The `simple` layer exists because the first version did not have
    one. A pharmacology professor read the protein-preparation step and
    said that even MD students would not follow it -- and he was right:
    it was written for someone who already knew what a crystal
    structure was. `simple` assumes none of that. It uses everyday
    comparisons and no jargon at all, and it is the layer shown first.
    Everything else is for the reader who wants more.
    """

    button: str                       # one sentence, on the button itself
    panel: str                        # 150-250 words, "Explain this step"
    not_mean: list[str]               # "What this does not mean" -- Section 7.2
    simple: str = ""                  # plain words, no jargon, shown first
    rule: str = ""                    # THE RULE -- see below
    rule_detail: list[str] = field(default_factory=list)
    # The decision rule for this step, and the detail under it.
    #
    # Added last, and it should have been first. Every other field here
    # explains what the step IS. None of them said how to CHOOSE, and
    # every step of this workflow is a choice: which structure, which
    # chains, where the box goes, which molecules.
    #
    # Found the hard way. A pharmacology professor working through SV2A
    # picked the apo structure at step 2 - no ligand in it. That one
    # choice cost him the known binding site, forced automatic pocket
    # detection, and put the box 13 angstroms from where levetiracetam
    # actually binds. Everything after it was confidently wrong. The
    # teaching text had told him what a crystal structure was. It had
    # not told him to prefer one with a drug already inside.
    #
    # `rule` is one sentence, imperative, and it is the first thing on
    # the screen - above everything, never collapsed. `rule_detail` is
    # the ordering, the thresholds and the exceptions.
    learn_more: list[tuple[str, str]] = field(default_factory=list)
    # learn_more entries are (label, url)


# ----------------------------------------------------------------------
# Provenance -- specification Section 5.3 / P5
# ----------------------------------------------------------------------

@dataclass
class ProvenanceEntry:
    """
    One line of the reproducibility record.

    Every number the software produces must be traceable to the tool
    and parameters that produced it. This is what lets a supervisor or
    an examiner check the work, and it is what makes the software
    citable.
    """

    stage: str
    tool: str
    version: str
    parameters: dict[str, Any] = field(default_factory=dict)
    seed: int | None = None
    source_release: str = ""          # e.g. "PDB 2026-10-21", "ChEMBL 35"
    citation: str = ""
    timestamp: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat(timespec="seconds")
    )


# ----------------------------------------------------------------------
# Stage results
# ----------------------------------------------------------------------

@dataclass
class StageResult:
    """
    What a stage produced.

    `control` is deliberately separate from `values`. Specification
    Section 9 requires that controls cannot be switched off or dropped
    from the report, so they are a first-class field rather than
    something optional mixed in with the results.
    """

    success: bool
    summary: str = ""                              # one line for the user
    values: dict[str, Any] = field(default_factory=dict)
    artefacts: list[Path] = field(default_factory=list)
    control: "ControlOutcome | None" = None
    message: str = ""                              # shown if success is False


class Verdict(str, Enum):
    PASS = "pass"
    WARN = "warn"
    FAIL = "fail"
    NOT_APPLICABLE = "not_applicable"


@dataclass
class ControlOutcome:
    """
    The check that tells the user whether the step is trustworthy.

    Examples: re-docking RMSD against the crystal pose; a decoy
    baseline for a docking score; an applicability-domain indicator
    for an ADMET prediction.

    `plain_verdict` is what the user reads. It must state the
    conclusion, not the number -- "this docking setup is not reliable
    for this pocket", not "RMSD = 3.4".
    """

    name: str
    verdict: Verdict
    plain_verdict: str
    detail: dict[str, Any] = field(default_factory=dict)
    blocks_progress: bool = False
    what_to_do: str = ""
    # What the user should actually DO about it.
    #
    # Added because a check that says "the resolution is 3.31
    # angstroms, which is too coarse to place side chains reliably"
    # tells an MD student nothing they can act on. They do not know
    # whether to stop, go back, pick something else, or carry on. A
    # verdict without a next action is not a check -- it is a worry.


# ----------------------------------------------------------------------
# The contract itself
# ----------------------------------------------------------------------

class Stage(ABC):
    """
    Base class for every stage of the workflow.

    Subclasses set the class attributes and implement the five methods.
    Nothing else in the application needs to know whether a stage runs
    locally, calls an API, opens a website or writes a cloud job.
    """

    stage_id: str = ""        # e.g. "02_structure"
    title: str = ""           # e.g. "Obtain the structure"
    tier: Tier = Tier.BUNDLED
    requires: list[str] = []  # stage_ids that must complete first

    def __init__(self, project: "Project") -> None:  # noqa: F821
        self.project = project
        self._provenance: list[ProvenanceEntry] = []

    # -- the five methods ------------------------------------------------

    @abstractmethod
    def check(self) -> CheckResult:
        """Are the inputs present? Is the service reachable?"""

    @abstractmethod
    def explain(self) -> Teaching:
        """Teaching content. Pure data; no computation, no network."""

    @abstractmethod
    def run(self, **kwargs: Any) -> StageResult:
        """
        Execute the step.

        For a STATION tier stage, this prepares the input and opens the
        handoff screen rather than computing anything.
        """

    def ingest(self, path: Path) -> StageResult:
        """
        Accept a result file the user brought back from an external
        service, validate it, parse it and merge it into the project.

        Only STATION and EXPORT stages need this. The default refuses
        politely rather than raising.
        """
        return StageResult(
            success=False,
            message=(
                f"{self.title} does not take a file from outside - "
                "it runs inside this software."
            ),
        )

    def record(self, entry: ProvenanceEntry) -> None:
        """Append to the reproducibility record for this stage."""
        self._provenance.append(entry)

    @property
    def provenance(self) -> list[ProvenanceEntry]:
        return list(self._provenance)

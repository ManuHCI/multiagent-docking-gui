"""
Stage 4 -- Identify the binding site.

Docking needs a box: somewhere to look. Choosing it is the step that
most often decides, quietly, whether the whole exercise means anything.

Three routes, in order of how much they can be trusted:

  (a) the site a co-crystallised ligand already occupies -- by far the
      best, because the pocket is proven to open and you get a known
      pose to validate against;
  (b) residues named in the literature;
  (c) geometric detection -- find the biggest crevice and hope.

The control measures two things about the chosen box. Its size,
because a big box does not find more; it finds worse. And how enclosed
it is, because a box floating in solvent or straddling the surface
will still return poses and scores, and they will be meaningless.

Runs offline. Numpy only.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from core.contract import (
    CheckResult,
    ControlOutcome,
    ProvenanceEntry,
    Readiness,
    Stage,
    StageResult,
    Teaching,
    Tier,
    Verdict,
)

# A box bigger than this is doing blind docking whether or not the user
# calls it that. 30 x 30 x 30 angstroms.
LARGE_BOX = 27000.0
# AutoDock Vina's own documentation warns against going much past this.
VERY_LARGE_BOX = 64000.0      # 40 x 40 x 40

WATER = {"HOH", "WAT", "DOD", "H2O"}


@dataclass
class Box:
    """A docking search space."""

    center: tuple[float, float, float]
    size: tuple[float, float, float]
    source: str = ""
    note: str = ""

    @property
    def volume(self) -> float:
        return round(self.size[0] * self.size[1] * self.size[2], 1)

    def contains(self, point) -> bool:
        return all(
            abs(point[i] - self.center[i]) <= self.size[i] / 2.0
            for i in range(3)
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "center_x": round(self.center[0], 3),
            "center_y": round(self.center[1], 3),
            "center_z": round(self.center[2], 3),
            "size_x": round(self.size[0], 1),
            "size_y": round(self.size[1], 1),
            "size_z": round(self.size[2], 1),
            "volume_A3": self.volume,
            "source": self.source,
        }


# ----------------------------------------------------------------------
# Reading coordinates
# ----------------------------------------------------------------------

def read_atoms(pdb_text: str) -> dict[str, Any]:
    """Protein heavy atoms, and every het group with its coordinates."""
    protein: list[tuple[float, float, float]] = []
    hets: dict[str, list[tuple[float, float, float]]] = {}
    residues: dict[tuple[str, str], list[tuple[float, float, float]]] = {}

    for line in pdb_text.splitlines():
        rec = line[:6]
        if rec not in ("ATOM  ", "HETATM"):
            continue
        element = line[76:78].strip().upper() or line[12:16].strip()[:1]
        if element == "H":
            continue
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue

        resname = line[17:20].strip().upper()
        chain = line[21]
        resseq = line[22:27].strip()

        if rec == "ATOM  ":
            protein.append(xyz)
            residues.setdefault((chain, resseq), []).append(xyz)
        elif resname not in WATER:
            hets.setdefault(f"{resname} {chain}{resseq}", []).append(xyz)

    return {
        "protein": np.array(protein, dtype=float) if protein
        else np.zeros((0, 3)),
        "hets": {k: np.array(v, dtype=float) for k, v in hets.items()},
        "residues": {k: np.array(v, dtype=float) for k, v in residues.items()},
    }


# ----------------------------------------------------------------------
# Building a box
# ----------------------------------------------------------------------

def box_around(points: np.ndarray, padding: float = 4.0,
               minimum: float = 18.0, source: str = "") -> Box:
    """
    A box enclosing these points, with room for a ligand to move.

    The padding is small on purpose. A common habit is to pad by 10 or
    15 angstroms "to be safe", which does the opposite: it enlarges the
    search space without adding anywhere real to bind, so the same
    sampling effort is spread more thinly and the result gets worse.
    """
    lo = points.min(axis=0) - padding
    hi = points.max(axis=0) + padding
    centre = (lo + hi) / 2.0
    size = np.maximum(hi - lo, minimum)
    return Box(
        center=tuple(float(c) for c in centre),
        size=tuple(float(s) for s in size),
        source=source,
    )


# ----------------------------------------------------------------------
# Grid geometry
#
# Pocket detection needs to ask, thousands of times, "is there protein
# in this direction?". Doing that against the atom list is far too slow
# on a laptop, so the protein is burned once into a binary voxel grid
# and every later question becomes an array lookup.
# ----------------------------------------------------------------------

# The detector must never suggest a box that the control on the next
# screen would then warn about. A suggestion the software immediately
# argues with is worse than no suggestion: the user has no way to tell
# whether the tool or the protein is at fault. So the ceiling here is
# the control's own warning threshold, not a separate number that
# could drift away from it.
MAX_DETECTED_BOX = LARGE_BOX


def _occupancy(protein: np.ndarray, lo: np.ndarray, hi: np.ndarray,
               spacing: float, radius: float) -> np.ndarray:
    """Mark every voxel within `radius` of a protein atom as solid."""
    shape = tuple(int(math.ceil((hi[i] - lo[i]) / spacing)) + 1
                  for i in range(3))
    solid = np.zeros(shape, dtype=bool)

    reach = int(math.ceil(radius / spacing))
    offsets = np.arange(-reach, reach + 1)
    dx, dy, dz = np.meshgrid(offsets, offsets, offsets, indexing="ij")
    within = (dx ** 2 + dy ** 2 + dz ** 2) * spacing ** 2 <= radius ** 2
    shell = np.stack([dx[within], dy[within], dz[within]], axis=1)

    base = np.round((protein - lo) / spacing).astype(int)
    for off in shell:
        p = base + off
        ok = np.all((p >= 0) & (p < np.array(shape)), axis=1)
        p = p[ok]
        solid[p[:, 0], p[:, 1], p[:, 2]] = True
    return solid


def _voxel_index(points: np.ndarray, lo: np.ndarray, spacing: float,
                 shape: tuple[int, ...]) -> np.ndarray:
    idx = np.round((points - lo) / spacing).astype(int)
    for axis in range(3):
        np.clip(idx[:, axis], 0, shape[axis] - 1, out=idx[:, axis])
    return idx


def _within(points: np.ndarray, lo: np.ndarray, spacing: float,
            solid: np.ndarray, inner: float, outer: float) -> np.ndarray:
    """
    True where a point has protein between `inner` and `outer` away.

    Used to keep only the shell of empty space that actually touches
    the protein, rather than bulk water far from it.
    """
    hit = np.zeros(points.shape[0], dtype=bool)
    steps = np.arange(inner, outer + 1e-9, spacing)
    directions = _directions()
    for d in directions:
        for t in steps:
            idx = _voxel_index(points + d * t, lo, spacing, solid.shape)
            hit |= solid[idx[:, 0], idx[:, 1], idx[:, 2]]
        if hit.all():
            break
    return hit


def _directions() -> np.ndarray:
    """The 26 directions from a grid point to its neighbours."""
    out = []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                if (dx, dy, dz) != (0, 0, 0):
                    v = np.array([dx, dy, dz], dtype=float)
                    out.append(v / np.linalg.norm(v))
    return np.array(out)


def _buriedness(points: np.ndarray, lo: np.ndarray, spacing: float,
                solid: np.ndarray, reach: float = 11.0) -> np.ndarray:
    """
    For each point, the fraction of 26 directions that hit protein.

    This is the number that separates a pocket from the open surface,
    and the one the first version of the detector never computed
    before clustering. A buried cavity scores above about 0.8; a flat
    patch of surface scores around 0.4.
    """
    steps = np.arange(spacing, reach + 1e-9, spacing)
    score = np.zeros(points.shape[0])
    for d in _directions():
        blocked = np.zeros(points.shape[0], dtype=bool)
        for t in steps:
            idx = _voxel_index(points + d * t, lo, spacing, solid.shape)
            blocked |= solid[idx[:, 0], idx[:, 1], idx[:, 2]]
        score += blocked
    return score / 26.0


def _cluster(points: np.ndarray, spacing: float) -> list[np.ndarray]:
    """Group points into connected blobs by grid adjacency."""
    keys = {tuple(np.round(p / spacing).astype(int)): p for p in points}
    seen: set = set()
    groups: list[np.ndarray] = []
    neighbours = [(dx, dy, dz)
                  for dx in (-1, 0, 1)
                  for dy in (-1, 0, 1)
                  for dz in (-1, 0, 1)
                  if (dx, dy, dz) != (0, 0, 0)]
    for key in keys:
        if key in seen:
            continue
        stack, group = [key], []
        seen.add(key)
        while stack:
            k = stack.pop()
            group.append(keys[k])
            for dx, dy, dz in neighbours:
                n = (k[0] + dx, k[1] + dy, k[2] + dz)
                if n in keys and n not in seen:
                    seen.add(n)
                    stack.append(n)
        groups.append(np.array(group))
    groups.sort(key=len, reverse=True)
    return groups


def enclosure(box: Box, protein: np.ndarray, spacing: float = 1.5) -> dict:
    """
    How much of this box is a pocket, and how enclosed is it?

    Grid the box. For each point, find the distance to the nearest
    protein heavy atom:

        under 2.4 A   inside the protein -- a ligand atom cannot go here
        2.4 to 6.0 A  usable -- close enough to touch, far enough to fit
        over 6.0 A    open solvent

    Then, for the usable points, cast 26 rays and count how many hit
    protein within 12 A. A real pocket is surrounded; a patch on the
    surface is not. This is what separates a pocket from a dent.
    """
    if protein.shape[0] == 0:
        return {"usable_fraction": 0.0, "solvent_fraction": 1.0,
                "buriedness": 0.0, "points": 0}

    cx, cy, cz = box.center
    sx, sy, sz = box.size
    axes = [
        np.arange(c - s / 2.0, c + s / 2.0 + 1e-9, spacing)
        for c, s in ((cx, sx), (cy, sy), (cz, sz))
    ]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
    if grid.shape[0] > 120_000:            # keep it quick on a laptop
        grid = grid[:: int(math.ceil(grid.shape[0] / 120_000))]

    # Only protein atoms near the box matter.
    half = np.array([sx, sy, sz]) / 2.0 + 8.0
    near = protein[np.all(np.abs(protein - np.array(box.center)) <= half,
                          axis=1)]
    if near.shape[0] == 0:
        return {"usable_fraction": 0.0, "solvent_fraction": 1.0,
                "buriedness": 0.0, "points": int(grid.shape[0])}

    # Nearest-atom distance, in blocks to keep memory sane.
    nearest = np.empty(grid.shape[0])
    block = 4000
    for i in range(0, grid.shape[0], block):
        d = np.linalg.norm(grid[i:i + block, None, :] - near[None, :, :],
                           axis=2)
        nearest[i:i + block] = d.min(axis=1)

    inside = nearest < 2.4
    usable = (nearest >= 2.4) & (nearest <= 6.0)
    solvent = nearest > 6.0

    # Buriedness of the usable points.
    directions = []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                if (dx, dy, dz) != (0, 0, 0):
                    v = np.array([dx, dy, dz], dtype=float)
                    directions.append(v / np.linalg.norm(v))
    directions = np.array(directions)

    pocket_points = grid[usable]
    if pocket_points.shape[0] == 0:
        buried = 0.0
    else:
        sample = pocket_points
        if sample.shape[0] > 400:
            step = int(math.ceil(sample.shape[0] / 400))
            sample = sample[::step]
        steps = np.arange(2.0, 12.1, 1.5)
        hits = np.zeros(sample.shape[0])
        for d in directions:
            probes = (sample[:, None, :]
                      + d[None, None, :] * steps[None, :, None])
            flat = probes.reshape(-1, 3)
            mind = np.empty(flat.shape[0])
            for i in range(0, flat.shape[0], block):
                dd = np.linalg.norm(flat[i:i + block, None, :]
                                    - near[None, :, :], axis=2)
                mind[i:i + block] = dd.min(axis=1)
            blocked = (mind.reshape(sample.shape[0], -1) < 2.4).any(axis=1)
            hits += blocked
        buried = float(np.mean(hits / len(directions)))

    n = float(grid.shape[0])
    return {
        "usable_fraction": round(float(usable.sum()) / n, 3),
        "inside_protein_fraction": round(float(inside.sum()) / n, 3),
        "solvent_fraction": round(float(solvent.sum()) / n, 3),
        "buriedness": round(buried, 3),
        "points": int(n),
        "usable_volume_A3": round(float(usable.sum()) * spacing ** 3, 1),
    }


# ----------------------------------------------------------------------
# The stage
# ----------------------------------------------------------------------

class SiteStage(Stage):
    stage_id = "04_site"
    title = "Identify the binding site"
    tier = Tier.BUNDLED
    requires = ["03_receptor"]

    def check(self) -> CheckResult:
        if not self.project.is_complete("03_receptor"):
            return CheckResult(
                readiness=Readiness.NEEDS_INPUT,
                reason="The protein has not been prepared yet.",
                fix="Complete the protein preparation step first.",
            )
        return CheckResult(readiness=Readiness.READY)

    def explain(self) -> Teaching:
        return Teaching(
            button="Decide where on the protein to dock.",
            panel=(
                "Docking does not search the whole protein. You give it a "
                "box, and it looks inside that box. Choosing the box is "
                "therefore a scientific decision, not a technical one, and "
                "it is where a great many published studies go wrong "
                "without anyone noticing.\n\n"
                "The instinct is to make the box large, to avoid missing "
                "anything. A larger search space can require more search "
                "effort and can make sampling less reliable at unchanged "
                "settings. Use a box that covers the intended site and "
                "allows the ligand to move, then check the result.\n\n"
                "Worse, a big box usually contains the largest crevice on "
                "the protein, which is often not the real site at all. "
                "Docking will fill that crevice happily, and because scores "
                "reward contact, a ligand buried in any deep groove scores "
                "well. If you have a structure with a drug already bound, "
                "use that site. It is the strongest evidence available that "
                "the pocket is real."
            ),
            simple=(
                "A protein is very big. A drug is very small. Binding "
                "can occur at specific regions, and the same protein "
                "can have more than one site. Here we investigate one "
                "site at a time.\n\n"
                "So before docking, you must tell the computer WHERE to "
                "look. You do this by drawing an invisible box around "
                "that groove. The computer will try your drug only "
                "inside this box.\n\n"
                "A common mistake: making the box very large so you do "
                "not miss anything.\n\n"
                "Why this can be a problem: a larger space needs more "
                "searching. With unchanged settings, the computer may "
                "miss useful poses or search an unintended region.\n\n"
                "If your structure already had a drug in it, the "
                "software can use that spot as a starting point. Check "
                "that it is the site you intend to study."
            ),
            not_mean=[
                "A box that encloses a good score is not thereby the right "
                "site. The score function rewards contact, so any deep "
                "groove can score well whether or not anything binds "
                "there in life.",
                "An automatically detected pocket is a geometric "
                "observation, not a biological one. The largest cavity on "
                "a protein is frequently a crystallographic artefact, an "
                "interface between two chains, or just where the fold "
                "happens to be loose.",
                "Choosing the site where a known drug binds does not prove "
                "your molecule binds there too. It provides evidence "
                "for that site's shape in the measured complex.",
                "A box is a search space, not a claim about the protein. "
                "Nothing here tells you whether the site is druggable - "
                "only where to look.",
            ],
            learn_more=[
                ("AutoDock Vina - choosing the search space",
                 "https://autodock-vina.readthedocs.io/en/latest/docking_basic.html"),
                ("Search-space size and sampling guidance",
                 "https://autodock-vina.readthedocs.io/en/latest/faq.html"),
            ],
        )

    # -- the three routes --------------------------------------------------

    def from_ligand(self, atoms: dict, ligand_label: str,
                    padding: float = 4.0) -> Box | None:
        coords = atoms["hets"].get(ligand_label)
        if coords is None or len(coords) == 0:
            return None
        box = box_around(coords, padding=padding,
                         source=f"the site occupied by {ligand_label}")
        box.note = (
            "This is the strongest choice available: something is already "
            "bound here, so the pocket is known to open, and the docking "
            "can be checked against the real pose."
        )
        return box

    def from_residues(self, atoms: dict, chain: str,
                      resseqs: list[str], padding: float = 5.0) -> Box | None:
        picked = [atoms["residues"][(chain, r)]
                  for r in resseqs if (chain, r) in atoms["residues"]]
        if not picked:
            return None
        box = box_around(np.vstack(picked), padding=padding,
                         source=f"residues {', '.join(resseqs)} of chain {chain}")
        box.note = (
            "Built from residues you named. It is only as good as the "
            "source you took them from - check that the numbering matches "
            "this structure, because it often does not."
        )
        return box

    def detect(self, atoms: dict, top: int = 3) -> list[Box]:
        """
        Find enclosed cavities on the protein surface.

        The first version of this was wrong in a way that produced a
        confident answer, which is the worst way to be wrong. It took
        every grid point sitting 2.8-5.5 angstroms from the protein
        and clustered them by adjacency. But that shell of points is
        one continuous layer wrapped around the whole molecule, so
        adjacency merged it into a single cluster covering the entire
        surface. On a membrane transporter it returned a 70 x 70 x 18
        angstrom slab and called it "pocket 1".

        The missing step was buriedness. A pocket is a cavity with
        protein on most sides; a point on the open surface has protein
        on one side and water everywhere else. Filtering on that
        before clustering breaks the shell apart and leaves the actual
        cavities behind.

        It is still crude, and still labelled as such wherever it
        appears. fpocket and P2Rank do this properly and neither
        packages cleanly for Windows. What this will now do reliably
        is find nothing rather than find the whole protein.
        """
        protein = atoms["protein"]
        if protein.shape[0] < 50:
            return []

        spacing = 1.2
        lo = protein.min(axis=0) - 4.0
        hi = protein.max(axis=0) + 4.0
        solid = _occupancy(protein, lo, hi, spacing, radius=2.6)

        axes = [np.arange(lo[i], hi[i], spacing) for i in range(3)]
        grid = np.stack(np.meshgrid(*axes, indexing="ij"),
                        axis=-1).reshape(-1, 3)

        # Candidate points: empty space that a ligand atom could sit
        # in, close enough to the protein to be touching it.
        idx = _voxel_index(grid, lo, spacing, solid.shape)
        free = ~solid[idx[:, 0], idx[:, 1], idx[:, 2]]
        near_surface = _within(grid, lo, spacing, solid, 2.6, 5.6)
        candidate = grid[free & near_surface]
        if candidate.shape[0] == 0:
            return []

        # The step that was missing. A point in a real pocket sees
        # protein in most directions; a point on the open surface does
        # not.
        buried = _buriedness(candidate, lo, spacing, solid, reach=11.0)
        enclosed = candidate[buried >= 0.68]
        if enclosed.shape[0] < 20:
            return []
        buried = buried[buried >= 0.68]

        clusters = _cluster(enclosed, spacing)

        scored: list[tuple[float, Box, dict]] = []
        # Verification grids each candidate, so only the largest few
        # are examined. A pocket that is not among the biggest handful
        # of enclosed cavities is not one this crude a method should
        # be putting in front of anybody.
        for group in clusters[:8]:
            if group.shape[0] < 25:
                continue
            box = box_around(group, padding=2.0, source="")
            # A pocket a drug can occupy is tens of cubic angstroms
            # across, not thousands. Anything larger is a groove, a
            # channel, or the whole surface -- and handing it over as
            # "pocket 1" is what produced the 88,200 A3 box.
            if box.volume > MAX_DETECTED_BOX:
                continue
            longest = max(box.size)
            shortest = max(min(box.size), 1e-6)
            if longest / shortest > 4.0:
                # A slab or a long channel. Docking into it searches a
                # region the ligand cannot occupy as a whole.
                continue
            # Verify the candidate by running the control that will
            # judge it on the next screen, and offer it only if it
            # passes cleanly. Using the control itself rather than a
            # second set of thresholds means the two can never drift
            # apart -- an earlier attempt at this filtered on a
            # stricter usable-space figure than the control used, and
            # threw away a cavity that was 91 per cent enclosed.
            shape = enclosure(box, protein)
            if self._box_control(box, shape, atoms).verdict is not Verdict.PASS:
                continue

            volume = group.shape[0] * spacing ** 3
            scored.append((volume * shape["buriedness"], box,
                           {"points": int(group.shape[0]),
                            "volume_A3": round(volume, 0),
                            "buriedness": shape["buriedness"]}))

        if not scored:
            return []

        scored.sort(key=lambda t: t[0], reverse=True)
        boxes = []
        for i, (_, box, info) in enumerate(scored[:top], start=1):
            box.source = f"automatically detected pocket {i}"
            box.note = (
                f"Found by geometry alone: an enclosed cavity of about "
                f"{info['volume_A3']:.0f} cubic angstroms, "
                f"{info['buriedness']:.0%} surrounded by protein. This "
                "says a cavity exists, not that anything binds in it. "
                "Check it against what is known about your protein "
                "before trusting it."
            )
            boxes.append(box)
        return boxes

    # -- run ---------------------------------------------------------------

    def run(
        self,
        pdb_path: str | Path = "",
        ligand_label: str = "",
        chain: str = "",
        residues: list[str] | None = None,
        padding: float = 4.0,
        box: dict | None = None,
        **kwargs: Any,
    ) -> StageResult:
        path = Path(pdb_path)
        if not path.exists():
            return StageResult(
                success=False,
                message=f"Could not find the prepared structure {path.name}.",
            )

        atoms = read_atoms(path.read_text(errors="ignore"))
        alternatives: list[Box] = []

        if box:
            chosen = Box(
                center=(box["center_x"], box["center_y"], box["center_z"]),
                size=(box["size_x"], box["size_y"], box["size_z"]),
                source="entered by you",
            )
        elif ligand_label:
            chosen = self.from_ligand(atoms, ligand_label, padding)
            if chosen is None:
                return StageResult(
                    success=False,
                    message=(
                        f"No group called {ligand_label} was found in the "
                        "prepared structure. It may have been removed "
                        "during preparation."
                    ),
                )
        elif residues:
            chosen = self.from_residues(atoms, chain or "A", residues, padding)
            if chosen is None:
                return StageResult(
                    success=False,
                    message=(
                        "None of those residue numbers were found in chain "
                        f"{chain or 'A'}. Structure numbering often differs "
                        "from the numbering in a paper - check against the "
                        "file."
                    ),
                )
        else:
            found = self.detect(atoms)
            if not found:
                return StageResult(
                    success=False,
                    message=(
                        "No enclosed pocket was found on this structure. "
                        "That is a real answer, not a crash, and it is "
                        "common - this detector only finds cavities that "
                        "are surrounded by protein on most sides, and "
                        "plenty of real binding sites are shallower than "
                        "that, or sit between two chains, or only open "
                        "up when a drug arrives.\n\n"
                        "You have three ways on, best first:\n\n"
                        "1. Use a bound drug. Go back to the Structure "
                        "step and pick an entry that has a molecule in "
                        "it; that molecule's position is the pocket, "
                        "measured rather than guessed.\n\n"
                        "2. Name the residues. Change 'Where should the "
                        "computer look?' to 'I know which residues line "
                        "the pocket' and type them in. Any paper on "
                        "your protein's active site lists them.\n\n"
                        "3. Type the centre. If you have coordinates "
                        "from PyMOL or from a paper, choose 'I will "
                        "type the box myself'."
                    ),
                )
            chosen, alternatives = found[0], found[1:]

        shape = enclosure(chosen, atoms["protein"])

        # How big is the molecule that defined this box? It turns a
        # bare percentage into something a pharmacologist can act on:
        # "11% usable" means nothing on its own, but "the cavity is
        # about three times the size of levetiracetam" is a sentence
        # about whether a drug fits.
        occupant = ""
        occupant_atoms = 0
        if ligand_label and ligand_label in atoms["hets"]:
            occupant = ligand_label
            occupant_atoms = int(len(atoms["hets"][ligand_label]))

        self.record(
            ProvenanceEntry(
                stage=self.stage_id,
                tool="binding site selection",
                version="built in",
                parameters={**chosen.as_dict(), "padding": padding},
                citation="",
            )
        )
        for e in self.provenance:
            self.project.add_provenance(e)

        return StageResult(
            success=True,
            summary=(
                f"Box {chosen.size[0]:.0f} x {chosen.size[1]:.0f} x "
                f"{chosen.size[2]:.0f} A ({chosen.volume:.0f} A3) from "
                f"{chosen.source}"
            ),
            values={
                "box": chosen.as_dict(),
                "note": chosen.note,
                "shape": shape,
                "alternatives": [b.as_dict() for b in alternatives],
                "available_ligands": sorted(atoms["hets"].keys()),
            },
            control=self._box_control(chosen, shape, atoms,
                                      occupant, occupant_atoms),
        )

    # -- the control -------------------------------------------------------

    @staticmethod
    def _box_control(box: Box, shape: dict, atoms: dict,
                     occupant: str = "",
                     occupant_atoms: int = 0) -> ControlOutcome:
        """
        Specification Stage 4: refuse an oversized box without an
        explicit acknowledgement.

        Two separate questions, because a box can fail either way. Is
        it too big -- spreading the search until nothing is examined
        properly? And is it in the right kind of place -- a pocket
        rather than open water or solid protein? A box can be a
        perfectly reasonable size and still be floating beside the
        protein rather than in it.
        """
        # What the usable space means for a molecule of this size.
        #
        # "11% of the box is space a ligand could occupy" is a true
        # statement that tells a pharmacologist nothing. The number
        # that matters is how that cavity compares with the drug: a
        # pocket two or three times the size of its occupant is a
        # snug, real binding site, and the same cavity is far too
        # small for a molecule twice as big. That second half is what
        # decides whether an analogue is worth docking here at all.
        fit = ""
        cavity = box.volume * shape["usable_fraction"]
        if occupant_atoms >= 6 and cavity > 0:
            # About 20 cubic angstroms per heavy atom is the usual
            # rough figure for the volume an organic molecule occupies.
            occupant_volume = occupant_atoms * 20.0
            ratio = cavity / occupant_volume
            name = occupant.split()[0] if occupant else "the bound molecule"
            if ratio < 1.5:
                fit = (
                    f" The cavity is about {cavity:.0f} cubic angstroms - "
                    f"only {ratio:.1f} times the size of {name} itself, so "
                    "this is a tight pocket. Anything much larger will not "
                    "fit in this conformation."
                )
            elif ratio < 5:
                fit = (
                    f" The cavity holds about {cavity:.0f} cubic angstroms, "
                    f"roughly {ratio:.1f} times the size of {name}. That is "
                    "a snug fit, which is what a genuine binding site looks "
                    "like. A molecule more than about "
                    f"{int(occupant_atoms * ratio)} heavy atoms would "
                    "struggle to fit here."
                )
            else:
                fit = (
                    f" The cavity is about {cavity:.0f} cubic angstroms, "
                    f"{ratio:.0f} times the size of {name}. There is room "
                    "for a considerably larger molecule, and room for "
                    "yours to sit in more than one place - expect the "
                    "poses to disagree with each other."
                )

        detail = {
            "box": f"{box.size[0]:.0f} x {box.size[1]:.0f} x "
                   f"{box.size[2]:.0f} A",
            "volume_A3": box.volume,
            "source": box.source,
            "space a ligand could occupy":
                f"{shape['usable_fraction']:.0%} of the box "
                f"({cavity:.0f} A3)",
            "solid protein": f"{shape.get('inside_protein_fraction', 0):.0%}",
            "open solvent": f"{shape['solvent_fraction']:.0%}",
            "buriedness": f"{shape['buriedness']:.0%}",
        }

        problems: list[str] = []
        blocking = False
        fatal = False          # wrong on its own, not merely worth checking

        if box.volume > VERY_LARGE_BOX:
            problems.append(
                f"this box is {box.volume:.0f} cubic angstroms, which is "
                "blind docking across much of the protein whether or not "
                "you call it that - the search is spread so thin that the "
                "scores stop meaning anything"
            )
            blocking = True
            fatal = True
        elif box.volume > LARGE_BOX:
            problems.append(
                f"at {box.volume:.0f} cubic angstroms this box is large; "
                "the same search effort is spread over more space, so "
                "every pose is examined less carefully"
            )
            blocking = True

        if shape["solvent_fraction"] > 0.5:
            problems.append(
                f"{shape['solvent_fraction']:.0%} of the box is open water "
                "with no protein nearby - the box is mostly beside the "
                "protein rather than in a pocket"
            )
            blocking = True
            fatal = True

        if shape["buriedness"] < 0.35 and shape["usable_fraction"] > 0:
            problems.append(
                f"the space in this box is only {shape['buriedness']:.0%} "
                "enclosed, which looks like a surface dent rather than a "
                "pocket"
            )
            blocking = True

        if shape["usable_fraction"] < 0.05:
            problems.append(
                "almost none of this box is space a ligand could occupy - "
                "it is nearly all solid protein or empty water"
            )
            blocking = True
            fatal = True

        if problems:
            return ControlOutcome(
                name="Is this a sensible place to dock?",
                verdict=Verdict.FAIL if (fatal or len(problems) > 1)
                else Verdict.WARN,
                plain_verdict=(
                    "Check the box before going on: "
                    + "; and ".join(problems) + "."
                ),
                what_to_do=(
                    "Do not dock with this box. Set it again using "
                    "'Where should the computer look?' above, taking "
                    "the first of these that you can:\n\n"
                    "1. If the structure has a drug bound in it, "
                    "choose that option - the software will place a "
                    "correct box for you with nothing to decide. If "
                    "this structure has none, go back to the Structure "
                    "step and look for an entry that does.\n\n"
                    "2. Choose 'I know which residues line the pocket' "
                    "and type in the residue numbers from any paper on "
                    "your protein's active site, or from its UniProt "
                    "entry under 'Binding site'.\n\n"
                    "3. Choose 'I will type the box myself' if you "
                    "already have coordinates, with a size of about "
                    "22 x 22 x 22 angstroms."
                ),
                detail=detail,
                blocks_progress=blocking,
            )

        return ControlOutcome(
            name="Is this a sensible place to dock?",
            verdict=Verdict.PASS,
            plain_verdict=(
                f"A sensible box: {box.volume:.0f} cubic angstroms, "
                f"{shape['usable_fraction']:.0%} of it space a ligand could "
                f"occupy, and {shape['buriedness']:.0%} enclosed by protein. "
                "That is the shape of a real pocket."
                + fit
            ),
            what_to_do=(
                "Nothing to change - press Next. 'Enclosed by protein' is "
                "the important number here: a real drug pocket is a cave "
                "with walls on most sides. A flat patch on the surface "
                "would give a low number and would not hold a drug."
            ),
            detail=detail,
        )

"""
Stage 3 -- Prepare the protein.

A crystal structure is not ready for docking. It carries the chemistry
of the experiment as well as the chemistry of the protein: cryo-
protectant, buffer, precipitant, sometimes several copies of the whole
molecule. It has no hydrogens. It usually has gaps.

Cleaning that up means deleting things, and deleting the wrong thing is
the commonest silent error in published docking work. A zinc in the
active site of a metalloprotease, the heme in a cytochrome P450, the
magnesium that positions ATP in a kinase -- each looks exactly like a
stray ion in the file, and each is part of the binding site. Remove it
and the docking still runs, still produces a score, and the score is
worthless.

So this stage never deletes a metal or cofactor quietly. It classifies
everything it finds, measures how close each one sits to the ligand,
and puts the list in front of the user. Specification Stage 3.

Runs offline. PDBFixer and OpenMM.
"""

from __future__ import annotations

import io
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from core.install import how_to_install
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

# ----------------------------------------------------------------------
# What the things in a structure file actually are
#
# These lists are deliberately conservative. Anything not recognised is
# treated as possibly the ligand and kept, because a false keep costs
# the user a moment and a false delete costs them their result.
# ----------------------------------------------------------------------

WATER = {"HOH", "WAT", "DOD", "H2O"}

# Present because of how the crystal was grown or frozen. Safe to
# remove -- but still listed, by name, every time.
CRYSTALLISATION = {
    "GOL": "glycerol (antifreeze)",
    "EDO": "ethylene glycol (antifreeze)",
    "PEG": "polyethylene glycol (precipitant)",
    "PG4": "polyethylene glycol (precipitant)",
    "PGE": "polyethylene glycol (precipitant)",
    "MPD": "MPD (precipitant)",
    "SO4": "sulfate (from the buffer)",
    "PO4": "phosphate (from the buffer)",
    "ACT": "acetate (from the buffer)",
    "FMT": "formate (from the buffer)",
    "CIT": "citrate (from the buffer)",
    "TRS": "Tris (buffer)",
    "EPE": "HEPES (buffer)",
    "IMD": "imidazole (from purification)",
    "DMS": "DMSO (solvent)",
    "BME": "beta-mercaptoethanol (reducing agent)",
    "NO3": "nitrate (from the buffer)",
    "AZI": "azide (preservative)",
    "CL": "chloride ion",
    "NA": "sodium ion",
    "K": "potassium ion",
    "BR": "bromide ion",
    "IOD": "iodide ion",
    "CS": "caesium ion",

    # Lipids, detergents and sugars. These matter more than they look.
    # A cryo-EM structure of a membrane protein routinely has a dozen
    # of them modelled, and several are larger than the drug. Without
    # this list, "the largest organic molecule in the file" -- which is
    # how the bound drug is identified -- picks out a detergent
    # molecule or a cholesterol sitting on the outside of a
    # transmembrane helix, and the binding site is then placed in the
    # lipid bilayer. Everything after that is confident and wrong.
    "CLR": "cholesterol (membrane lipid)",
    "CHS": "cholesteryl hemisuccinate (detergent)",
    "Y01": "cholesteryl hemisuccinate (detergent)",
    "OLC": "monoolein (lipid, from the crystallisation matrix)",
    "OLA": "oleic acid (lipid)",
    "PLM": "palmitic acid (lipid)",
    "MYR": "myristic acid (lipid)",
    "STE": "stearic acid (lipid)",
    "POV": "POPC phospholipid (membrane)",
    "PTY": "phosphatidylethanolamine (membrane lipid)",
    "PEE": "phosphatidylethanolamine (membrane lipid)",
    "PC1": "phosphatidylcholine (membrane lipid)",
    "PCW": "phosphatidylcholine (membrane lipid)",
    "LMT": "lauryl maltoside (detergent)",
    "LMN": "lauryl maltose neopentyl glycol (detergent)",
    "DDQ": "decylamine-N,N-dimethyl oxide (detergent)",
    "D12": "dodecane (detergent tail)",
    "UND": "undecane (detergent tail)",
    "HEX": "hexane (detergent tail)",
    "DGA": "diacylglycerol (lipid)",
    "NAG": "N-acetylglucosamine (a sugar the protein is glycosylated with)",
    "NDG": "N-acetylglucosamine (glycosylation)",
    "BMA": "mannose (glycosylation)",
    "MAN": "mannose (glycosylation)",
    "FUC": "fucose (glycosylation)",
    "GAL": "galactose (glycosylation)",
    "SIA": "sialic acid (glycosylation)",
}

# Metals that are frequently catalytic. Never removed without the user
# being told, and never at all when they sit in the binding site.
METALS = {
    "ZN": "zinc", "MG": "magnesium", "MN": "manganese", "FE": "iron",
    "FE2": "iron(II)", "CA": "calcium", "CU": "copper", "CU1": "copper(I)",
    "NI": "nickel", "CO": "cobalt", "CD": "cadmium", "MO": "molybdenum",
    "W": "tungsten", "V": "vanadium", "MGF": "magnesium trifluoride",
    "ALF": "aluminium fluoride", "BEF": "beryllium fluoride",
}

# Cofactors and prosthetic groups. Part of the protein's chemistry.
COFACTORS = {
    "HEM": "heme", "HEC": "heme C", "HEA": "heme A", "SRM": "siroheme",
    "NAD": "NAD+", "NAI": "NADH", "NAP": "NADP+", "NDP": "NADPH",
    "FAD": "FAD", "FMN": "FMN", "SAM": "S-adenosylmethionine",
    "SAH": "S-adenosylhomocysteine", "ATP": "ATP", "ADP": "ADP",
    "AMP": "AMP", "ANP": "AMP-PNP (non-hydrolysable ATP)",
    "ACP": "AMP-PCP (non-hydrolysable ATP)", "GTP": "GTP", "GDP": "GDP",
    "GNP": "GMP-PNP", "PLP": "pyridoxal phosphate", "COA": "coenzyme A",
    "TPP": "thiamine pyrophosphate", "B12": "cobalamin",
    "BTN": "biotin", "THF": "tetrahydrofolate", "MTE": "molybdopterin",
    "PQQ": "PQQ", "F43": "coenzyme F430", "UDP": "UDP", "UTP": "UTP",
}

KIND_WATER = "water"
KIND_CRYSTAL = "crystallisation"
KIND_METAL = "metal"
KIND_COFACTOR = "cofactor"
KIND_LIGAND = "ligand"

# Anything within this of the ligand is treated as part of the site.
SITE_RADIUS = 6.0      # angstroms


@dataclass
class Het:
    """One non-protein thing found in the structure."""

    resname: str
    chain: str
    resseq: str
    kind: str
    description: str
    n_atoms: int
    coords: list[tuple[float, float, float]] = field(default_factory=list)
    distance_to_ligand: float | None = None

    @property
    def label(self) -> str:
        return f"{self.resname} {self.chain}{self.resseq}".strip()

    @property
    def in_site(self) -> bool:
        return (self.distance_to_ligand is not None
                and self.distance_to_ligand <= SITE_RADIUS)


# ----------------------------------------------------------------------
# Reading a structure without needing a heavyweight parser
# ----------------------------------------------------------------------

def gap_distances(fixer: Any, missing: dict,
                  ligand_coords: list[tuple[float, float, float]]
                  ) -> dict[tuple, float | None]:
    """
    How far each hole in the chain is from the bound drug.

    This is the question the gap warning used to hand back to the
    user: "see whether any missing piece is near your pocket". They
    cannot see that. It needs the coordinates of residues that are
    not in the file, opened in a program most of them do not have.

    But the software can answer it. A missing stretch has no
    coordinates by definition -- that is what missing means -- yet the
    residues on either side of it do, and they are covalently attached
    to it. If both ends of a hole are forty angstroms from the drug,
    nothing that belongs in that hole reaches the pocket. If an end is
    six angstroms away, the hole is part of the pocket wall.

    Returns the distance in angstroms for each gap, keyed the way
    PDBFixer keys them, or None where there is no ligand to measure
    against.
    """
    if not ligand_coords:
        return {key: None for key in missing}

    chains = list(fixer.topology.chains())
    positions = fixer.positions
    out: dict[tuple, float | None] = {}

    for key in missing:
        chain_index, insert_at = key[0], key[1]
        if chain_index >= len(chains):
            out[key] = None
            continue
        residues = list(chains[chain_index].residues())
        # The residues either side of the hole: the one the gap comes
        # after, and the one it comes before.
        flanks = [i for i in (insert_at - 1, insert_at)
                  if 0 <= i < len(residues)]
        best: float | None = None
        for i in flanks:
            for atom in residues[i].atoms():
                try:
                    p = positions[atom.index]
                    px, py, pz = (p.x * 10.0, p.y * 10.0, p.z * 10.0)
                except Exception:  # noqa: BLE001
                    continue
                for lx, ly, lz in ligand_coords:
                    d = ((px - lx) ** 2 + (py - ly) ** 2
                         + (pz - lz) ** 2) ** 0.5
                    if best is None or d < best:
                        best = d
        out[key] = round(best, 1) if best is not None else None
    return out


def classify(resname: str) -> tuple[str, str]:
    r = resname.strip().upper()
    if r in WATER:
        return KIND_WATER, "water"
    if r in METALS:
        return KIND_METAL, METALS[r]
    if r in COFACTORS:
        return KIND_COFACTOR, COFACTORS[r]
    if r in CRYSTALLISATION:
        return KIND_CRYSTAL, CRYSTALLISATION[r]
    return KIND_LIGAND, "possibly your ligand, or something unrecognised"


def survey(pdb_text: str, restrict_chains: list[str] | None = None) -> dict[str, Any]:
    """
    Read the file and say what is in it.

    Only HETATM residue names are classified. An ATOM record's atom
    name can be 'CA' for alpha-carbon, which is also the code for
    calcium -- confusing the two would be an entertaining way to delete
    an entire protein backbone.

    `restrict_chains` matters more than it looks. A crystal often holds
    several copies of the protein, each with its own ligand and metals.
    Pick the largest ligand across the whole file while keeping only
    chain A, and you can end up keeping chain B's cofactor: a receptor
    with its magnesium 70 angstroms away in empty space, and a deletion
    report that calls the real site "away from the binding site". That
    happened during development on this very structure, and the control
    said it was probably fine. Hence this argument.
    """
    chains: dict[str, int] = {}
    hets: dict[tuple, Het] = {}
    missing_residues = 0
    models = 0

    for line in pdb_text.splitlines():
        rec = line[:6]
        if rec == "MODEL ":
            models += 1
        elif rec == "ATOM  ":
            ch = line[21]
            chains[ch] = chains.get(ch, 0) + 1
        elif rec == "HETATM":
            resname = line[17:20].strip().upper()
            key = (resname, line[21], line[22:27].strip())
            if key not in hets:
                kind, desc = classify(resname)
                hets[key] = Het(resname=resname, chain=line[21],
                                resseq=line[22:27].strip(),
                                kind=kind, description=desc, n_atoms=0)
            hets[key].n_atoms += 1
            try:
                hets[key].coords.append((float(line[30:38]),
                                         float(line[38:46]),
                                         float(line[46:54])))
            except ValueError:
                pass
        elif line.startswith("REMARK 465") and len(line) > 20:
            tail = line[10:].strip()
            # Data lines look like "ALA A   123"; headers do not.
            parts = tail.split()
            if len(parts) >= 3 and len(parts[0]) == 3 and parts[0].isalpha():
                missing_residues += 1

    het_list = list(hets.values())

    # Confine everything to the chains being kept, if any were named.
    wanted = set(restrict_chains or [])
    if wanted:
        het_list = [h for h in het_list if h.chain in wanted]
        chains = {c: n for c, n in chains.items() if c in wanted}

    # The ligand: the largest unrecognised or cofactor-sized organic
    # molecule. Big enough to be drug-like, not a lone ion.
    candidates = [h for h in het_list
                  if h.kind in (KIND_LIGAND, KIND_COFACTOR) and h.n_atoms >= 10]
    ligand = max(candidates, key=lambda h: h.n_atoms) if candidates else None

    if ligand and ligand.coords:
        for h in het_list:
            if h is ligand or not h.coords:
                continue
            h.distance_to_ligand = round(min(
                ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)
                ** 0.5
                for a in h.coords for b in ligand.coords
            ), 2)

    return {
        "chains": chains,
        "hets": het_list,
        "ligand": ligand,
        "missing_residues": missing_residues,
        "models": models,
    }


# ----------------------------------------------------------------------
# The stage
# ----------------------------------------------------------------------

class ReceptorStage(Stage):
    stage_id = "03_receptor"
    title = "Prepare the protein"
    tier = Tier.BUNDLED
    requires = ["02_structure"]

    PH = 7.4

    def check(self) -> CheckResult:
        try:
            import pdbfixer  # noqa: F401
            import openmm  # noqa: F401
        except ImportError:
            return CheckResult(
                readiness=Readiness.BLOCKED,
                reason="The structure-repair tools are not installed "
                       "on this computer.",
                fix="Open Anaconda Prompt and run:  "
                    + how_to_install("openmm", "pdbfixer"),
            )
        return CheckResult(readiness=Readiness.READY)

    def explain(self) -> Teaching:
        return Teaching(
            button="Clean up the structure and get it ready for docking.",
            panel=(
                "What comes out of a crystallography experiment is a record "
                "of that experiment, not a tidy picture of a protein. It "
                "contains whatever was in the drop: glycerol to stop ice "
                "forming, sulfate or phosphate from the buffer, imidazole "
                "left over from purification. It often contains several "
                "copies of the protein, because that is how the crystal "
                "packed. Many X-ray structures omit most hydrogens because "
                "they are difficult to resolve. There may also be gaps "
                "where part of the chain "
                "was too mobile to be seen.\n\n"
                "All of that has to be sorted out, and sorting it out means "
                "deleting things. That is where the danger is. The zinc in "
                "a metalloprotease, the heme in a cytochrome P450, the "
                "magnesium that holds ATP in place in a kinase - in the file "
                "each of these looks exactly like a stray ion. Delete one "
                "and the docking still runs and still gives you a number. "
                "The number is meaningless, and nothing tells you.\n\n"
                "So nothing here is removed quietly. Everything found is "
                "named, with what it probably is and how far it sits from "
                "your ligand. Anything in the binding site is kept unless "
                "you say otherwise."
            ),
            simple=(
                "The structure file you just downloaded is like a "
                "photograph taken in a laboratory. It shows your protein "
                "- but it also shows extra things that were in the test "
                "tube and are not in the human body.\n\n"
                "What is in the file that should not be:\n\n"
                "- Antifreeze (glycerol), used to stop ice forming\n"
                "- Salt and buffer left over from the experiment\n"
                "- Often 2, 3 or 4 copies of the same protein\n\n"
                "What is missing from the file:\n\n"
                "- Many hydrogen atoms, which are difficult to locate "
                "in routine X-ray structures.\n"
                "- Some parts of the protein that were moving too much "
                "to photograph. These show up as holes.\n\n"
                "This step cleans all that up.\n\n"
                "ONE IMPORTANT WARNING: some proteins need a metal to "
                "work. ACE needs zinc. Cytochrome P450 needs iron. In "
                "the file, that metal looks just like leftover salt. If "
                "we throw it away, the docking still runs and still "
                "gives a nice-looking number - but the number is "
                "rubbish. So this software never removes a metal quietly. "
                "Anything removed is shown to you by name."
            ),
            not_mean=[
                "Adding hydrogens at pH 7.4 is an assumption, not a "
                "measurement. The pKa of a residue inside a pocket can be "
                "shifted by two units or more by its surroundings, and if a "
                "histidine in your site is charged when the software thinks "
                "it is neutral, the docking is wrong.",
                "A repaired loop is a guess at where missing atoms went. If "
                "the gap was near your binding site, treat anything built "
                "there as invented rather than observed.",
                "Minimisation removes bad contacts. It does not correct a "
                "wrong structure, and it will not turn a poor model into a "
                "good one.",
                "Keeping a water is a decision too. A water bridging the "
                "ligand and the protein is part of the site; removing it "
                "changes what can bind. Most docking runs delete all of "
                "them, which is a choice people rarely notice making.",
            ],
            learn_more=[
                ("Protein preparation - what it involves and why",
                 "https://pubmed.ncbi.nlm.nih.gov/23579614/"),
                ("PDBFixer documentation",
                 "https://github.com/openmm/pdbfixer"),
            ],
        )

    # -- run ---------------------------------------------------------------

    def run(
        self,
        pdb_path: str | Path = "",
        keep_chains: list[str] | None = None,
        keep_hets: list[str] | None = None,
        remove_waters: bool = True,
        repair_gaps: bool = True,
        max_gap_length: int = 3,
        minimise: bool = False,
        **kwargs: Any,
    ) -> StageResult:
        path = Path(pdb_path)
        if not path.exists():
            return StageResult(
                success=False,
                message=f"Could not find the structure file {path.name}.",
            )

        text = path.read_text(errors="ignore")
        found = survey(text, restrict_chains=keep_chains)

        if keep_chains and not found["chains"]:
            return StageResult(
                success=False,
                message=(
                    f"No chain named {', '.join(keep_chains)} was found in "
                    f"{path.name}. Check the chain letters in the structure."
                ),
            )

        # Default: keep metals and cofactors that sit in the binding
        # site, plus the ligand itself. Everything else recognised as
        # crystallisation chemistry goes.
        if keep_hets is None:
            keep_hets = [
                h.label for h in found["hets"]
                if h.kind in (KIND_METAL, KIND_COFACTOR) and h.in_site
            ]

        kept, removed = [], []
        for h in found["hets"]:
            if h.kind == KIND_WATER:
                (removed if remove_waters else kept).append(h)
            elif found["ligand"] is not None and h is found["ligand"]:
                kept.append(h)
            elif h.label in keep_hets:
                kept.append(h)
            else:
                removed.append(h)

        folder = self.project.step_folder(self.stage_id)
        out_path = folder / f"{path.stem}_prepared.pdb"

        try:
            report = self._fix(text, out_path, keep_chains,
                               [h.label for h in kept], repair_gaps,
                               max_gap_length,
                               ligand_coords=(found["ligand"].coords
                                              if found["ligand"] else None))
        except ModuleNotFoundError as exc:
            # Tell the truth about which kind of failure this is.
            #
            # An earlier version caught every exception here and reported
            # "the structure could not be prepared, it may be unusual or
            # damaged" -- for a missing Python package, about a pristine
            # file straight from the Protein Data Bank. A wrong diagnosis
            # delivered confidently is worse than no diagnosis: it sends
            # the user to inspect a structure that was never the problem.
            missing = getattr(exc, "name", "") or "a required tool"
            return StageResult(
                success=False,
                message=(
                    f"The tool '{missing}' is not installed, so the "
                    "structure cannot be prepared. This is nothing to do "
                    "with your structure. Install it with:\n"
                    f"    conda install -c conda-forge {missing}"
                ),
            )
        except Exception as exc:  # noqa: BLE001
            return StageResult(
                success=False,
                message=(
                    "The structure could not be prepared. It may be unusual "
                    f"or damaged. ({type(exc).__name__}: {exc})"
                ),
            )

        self.record(
            ProvenanceEntry(
                stage=self.stage_id,
                tool="PDBFixer",
                version=self._version(),
                parameters={
                    "pH": self.PH,
                    "waters_removed": remove_waters,
                    "gaps_repaired": repair_gaps,
                    "max_gap_length_filled": max_gap_length,
                    "kept": [h.label for h in kept],
                    "removed": [h.label for h in removed
                                if h.kind != KIND_WATER],
                    "chains_kept": keep_chains or "all",
                },
                citation="Eastman P et al. OpenMM 8. J Phys Chem B. "
                         "2024;128:109-116.",
            )
        )
        for e in self.provenance:
            self.project.add_provenance(e)

        return StageResult(
            success=True,
            summary=(
                f"Prepared {path.stem}: kept "
                f"{len([h for h in kept if h.kind != KIND_WATER])} "
                f"non-water component(s), removed "
                f"{len([h for h in removed if h.kind != KIND_WATER])}, "
                f"added {report.get('hydrogens_added', 0)} hydrogens"
            ),
            values={
                "prepared_file": str(out_path),
                "chains": found["chains"],
                "ligand": found["ligand"].label if found["ligand"] else "",
                "kept": [self._describe(h) for h in kept],
                "removed": [self._describe(h) for h in removed
                            if h.kind != KIND_WATER],
                "waters_removed": sum(1 for h in removed
                                      if h.kind == KIND_WATER),
                "missing_residues_in_file": found["missing_residues"],
                **report,
            },
            artefacts=[out_path],
            control=self._soundness_control(found, kept, removed, report),
        )

    # -- the actual repair -------------------------------------------------

    def _fix(self, text: str, out_path: Path, keep_chains: list[str] | None,
             keep_labels: list[str], repair_gaps: bool,
             max_gap_length: int,
             ligand_coords: list[tuple[float, float, float]] | None = None,
             ) -> dict[str, Any]:
        from openmm.app import PDBFile
        from pdbfixer import PDBFixer

        fixer = PDBFixer(pdbfile=io.StringIO(text))

        if keep_chains:
            wanted = set(keep_chains)
            drop = [i for i, c in enumerate(fixer.topology.chains())
                    if c.id not in wanted]
            if drop:
                fixer.removeChains(drop)

        fixer.findMissingResidues()
        all_gaps = dict(fixer.missingResidues)
        gaps = sum(len(v) for v in all_gaps.values())

        # Measure now, while the original residues are still in the
        # topology. After addMissingResidues the holes are filled and
        # the question can no longer be asked.
        to_ligand = gap_distances(fixer, all_gaps, ligand_coords or [])

        # Only short gaps are filled, and this is a decision about
        # honesty before it is one about speed.
        #
        # Rebuilding a three-residue break is interpolation; the ends
        # are pinned and there is not much room for the middle to be
        # wrong. Rebuilding a twenty-residue disordered loop is
        # invention. The residues are missing from the structure
        # because they were too mobile to see, and a model of them is
        # a guess presented in the same file format as a measurement.
        #
        # It is also what makes this stage usable. On the test
        # structure PDBFixer was asked to rebuild 116 gaps and had not
        # finished after fifteen minutes. A resident watching a frozen
        # window will kill it and conclude the software is broken.
        if repair_gaps:
            short = {k: v for k, v in all_gaps.items()
                     if len(v) <= max_gap_length}
            fixer.missingResidues = short
            filled = sum(len(v) for v in short.values())
            left_open = [
                {"length": len(v), "after_residue_index": k[1],
                 "angstroms_from_ligand": to_ligand.get(k)}
                for k, v in all_gaps.items() if len(v) > max_gap_length
            ]
        else:
            fixer.missingResidues = {}
            filled = 0
            left_open = [{"length": len(v), "after_residue_index": k[1],
                           "angstroms_from_ligand": to_ligand.get(k)}
                         for k, v in all_gaps.items()]

        # Delete only what is NOT on the keep list. PDBFixer's own
        # removeHeterogens is all-or-nothing, which is exactly the
        # behaviour that loses people their cofactor, so the selection
        # is made here instead.
        keep_set = set(keep_labels)
        to_delete = []
        for chain in fixer.topology.chains():
            for residue in chain.residues():
                name = residue.name.strip().upper()
                if name in WATER:
                    continue   # handled by the caller's water setting
                kind, _ = classify(name)
                if kind == "protein":
                    continue
                label = f"{name} {chain.id}{residue.id}".strip()
                if name in WATER or label in keep_set:
                    continue
                if _is_standard_residue(name):
                    continue
                to_delete.append(residue)

        if to_delete:
            modeller_ok = True
            try:
                from openmm.app import Modeller
                m = Modeller(fixer.topology, fixer.positions)
                m.delete(to_delete)
                fixer.topology, fixer.positions = m.topology, m.positions
            except Exception:  # noqa: BLE001
                modeller_ok = False
        else:
            modeller_ok = True

        fixer.findMissingAtoms()
        n_missing_atoms = sum(len(v) for v in fixer.missingAtoms.values())
        fixer.addMissingAtoms()

        before_h = fixer.topology.getNumAtoms()
        fixer.addMissingHydrogens(self.PH)
        after_h = fixer.topology.getNumAtoms()

        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as fh:
            PDBFile.writeFile(fixer.topology, fixer.positions, fh,
                              keepIds=True)

        return {
            "gaps_found": gaps,
            "gaps_repaired": filled,
            "gaps_left_open": left_open,
            "longest_gap_left": max([g["length"] for g in left_open], default=0),
            "nearest_gap_to_ligand": min(
                [g["angstroms_from_ligand"] for g in left_open
                 if g.get("angstroms_from_ligand") is not None],
                default=None),
            "missing_atoms_added": n_missing_atoms,
            "hydrogens_added": after_h - before_h,
            "final_atoms": after_h,
            "selective_deletion": modeller_ok,
        }

    # -- the control -------------------------------------------------------

    @staticmethod
    def _describe(h: Het) -> dict[str, Any]:
        return {
            "label": h.label, "what": h.description, "kind": h.kind,
            "atoms": h.n_atoms,
            "distance_to_ligand": h.distance_to_ligand,
        }

    @classmethod
    def _soundness_control(cls, found: dict, kept: list[Het],
                           removed: list[Het],
                           report: dict[str, Any]) -> ControlOutcome:
        """
        Is the receptor sound after preparation?

        This stage can damage a structure in two independent ways, and
        an earlier version of this control only looked for one of them.
        It passed a receptor with a 55-residue hole in it because
        nothing important had been *deleted* -- the question it was
        named after. The hole was in the report, and the control did
        not read it.

        A control has to cover every way its stage can fail, not the
        failure it is named after. So both are checked here and the
        worse verdict wins.
        """
        deletion = cls._deletion_control(found, kept, removed)
        holes = cls._gap_control(report)

        rank = {Verdict.PASS: 0, Verdict.NOT_APPLICABLE: 0,
                Verdict.WARN: 1, Verdict.FAIL: 2}
        worse = max((deletion, holes), key=lambda c: rank[c.verdict])
        other = holes if worse is deletion else deletion

        if rank[other.verdict] == 0 and rank[worse.verdict] == 0:
            return ControlOutcome(
                name="Is this structure sound after preparation?",
                verdict=Verdict.PASS,
                plain_verdict=f"{deletion.plain_verdict} {holes.plain_verdict}",
                # Both sub-checks passing means both want to say "press
                # Next", and printing that twice reads like a stutter.
                # When everything is fine, one sentence is the whole of
                # what the user needs.
                what_to_do=(
                    "Nothing to change - press Next. Both things that "
                    "go wrong at this step have been checked: nothing "
                    "the protein needs was deleted, and no hole in the "
                    "chain is near the pocket."
                ),
                detail={**deletion.detail, **holes.detail},
            )

        parts = [c.plain_verdict for c in (worse, other)
                 if rank[c.verdict] > 0]
        return ControlOutcome(
            name="Is this structure sound after preparation?",
            verdict=worse.verdict,
            plain_verdict=" ".join(parts),
            what_to_do=" ".join(
                c.what_to_do for c in (worse, other)
                if rank[c.verdict] > 0 and c.what_to_do
            ),
            detail={**deletion.detail, **holes.detail},
            blocks_progress=deletion.blocks_progress or holes.blocks_progress,
        )

    @staticmethod
    def _gap_control(report: dict[str, Any]) -> ControlOutcome:
        """
        What holes are left in the chain, and do they matter?

        A gap is not automatically fatal -- a disordered surface loop
        far from the pocket changes nothing. A gap beside the binding
        site changes everything, because docking will happily place a
        ligand in the empty space where the missing residues should be
        and report an excellent score for it.
        """
        left = report.get("gaps_left_open", []) or []
        longest = report.get("longest_gap_left", 0) or 0
        filled = report.get("gaps_repaired", 0) or 0
        total = report.get("gaps_found", 0) or 0

        nearest = report.get("nearest_gap_to_ligand")
        detail = {
            "amino acids missing from the file": total,
            "short holes rebuilt": filled,
            "holes left open": len(left),
            "longest hole left open": f"{longest} amino acids",
            "nearest hole to the bound drug":
                f"{nearest} A" if nearest is not None
                else "no bound drug to measure against",
        }
        if left:
            detail["each hole (length, distance from the drug)"] = "; ".join(
                f"{g['length']} aa at "
                + (f"{g['angstroms_from_ligand']} A"
                   if g.get("angstroms_from_ligand") is not None
                   else "unknown distance")
                for g in left[:8]
            )

        if not left:
            return ControlOutcome(
                name="Holes in the chain",
                verdict=Verdict.PASS,
                plain_verdict=(
                    "The chain is continuous."
                    if not total else
                    f"All {total} missing residue(s) were short enough to "
                    "rebuild reliably."
                ),
                what_to_do=(
                    "Nothing to do. The protein chain has no holes in it."
                ),
                detail=detail,
            )

        if longest >= 10:
            # The question is always the same: are the holes near the
            # pocket? An earlier version of this handed that question
            # back to the user, which was no help at all - answering it
            # means knowing where residues that are not in the file
            # would have been. The software can do better, because the
            # residues on either side of each hole are in the file and
            # are attached to it. So it measures, and says.
            far = 15.0
            if nearest is None:
                return ControlOutcome(
                    name="Holes in the chain",
                    verdict=Verdict.WARN,
                    plain_verdict=(
                        f"This structure has {len(left)} hole(s) in the "
                        f"chain, the largest {longest} amino acids long. "
                        "They have been left open rather than invented. "
                        "There is no bound drug in this structure, so the "
                        "software cannot work out whether the holes are "
                        "near your pocket."
                    ),
                    what_to_do=(
                        "Carry on to the Binding Site step and set your "
                        "box. Then come back and compare: the hole "
                        "positions are listed below. If one is beside "
                        "your box, pick a different structure - docking "
                        "will happily put a drug into the empty space "
                        "where the missing residues belong and award it a "
                        "good score."
                    ),
                    detail=detail,
                    blocks_progress=True,
                )

            if nearest >= far:
                return ControlOutcome(
                    name="Holes in the chain",
                    verdict=Verdict.PASS,
                    plain_verdict=(
                        f"There are {len(left)} hole(s) in the chain, the "
                        f"largest {longest} amino acids - but the nearest "
                        f"one is {nearest} angstroms from the bound drug, "
                        "which is a long way in protein terms. These are "
                        "floppy loops and tails on the outside, not part "
                        "of the pocket."
                    ),
                    what_to_do=(
                        "Nothing to do - press Next. Holes like these are "
                        "normal and they do not affect docking. What "
                        "would matter is a hole at the pocket itself, and "
                        "there is none."
                    ),
                    detail=detail,
                )

            return ControlOutcome(
                name="Holes in the chain",
                verdict=Verdict.FAIL if nearest < 8.0 else Verdict.WARN,
                plain_verdict=(
                    f"One of the holes in this chain comes within "
                    f"{nearest} angstroms of the bound drug. That is "
                    "close enough to be part of the pocket wall. Docking "
                    "into a pocket with a missing wall gives confident "
                    "scores for poses that sit in a space the real "
                    "protein fills."
                ),
                what_to_do=(
                    "Go back to the Structure step and pick a different "
                    "entry - look for one with fewer missing amino acids, "
                    "even if the picture is slightly blurrier. If this is "
                    "the only structure there is, you can continue, but "
                    "say plainly in your write-up that the binding site "
                    "was incompletely resolved."
                ),
                detail=detail,
                blocks_progress=True,
            )

        return ControlOutcome(
            name="Holes in the chain",
            verdict=Verdict.WARN,
            plain_verdict=(
                f"{len(left)} short hole(s) remain, the largest {longest} "
                "residues. Probably harmless unless one is at your binding "
                "site."
            ),
            what_to_do=(
                f"Carry on. The nearest hole is {nearest} angstroms from "
                "the bound drug, which is well clear of the pocket."
                if nearest is not None and nearest >= 15.0 else
                "Carry on, but check the hole positions below against "
                "your box at the next step."
            ),
            detail=detail,
        )

    @staticmethod
    def _deletion_control(found: dict, kept: list[Het],
                          removed: list[Het]) -> ControlOutcome:
        """
        Specification Stage 3: every deleted molecule listed by name.

        Progress is blocked when something that may belong in the
        binding site is on the removal list, because that is the error
        the whole stage exists to prevent.
        """
        important_removed = [
            h for h in removed
            if h.kind in (KIND_METAL, KIND_COFACTOR)
        ]
        in_site_removed = [h for h in important_removed if h.in_site]

        kept_important = [h for h in kept
                          if h.kind in (KIND_METAL, KIND_COFACTOR)]

        detail = {
            "ligand": found["ligand"].label if found["ligand"] else "none found",
            "kept": [f"{h.label} - {h.description}" for h in kept
                     if h.kind != KIND_WATER],
            "removed": [f"{h.label} - {h.description}" for h in removed
                        if h.kind != KIND_WATER],
            "waters removed": sum(1 for h in removed if h.kind == KIND_WATER),
            "chains": list(found["chains"].keys()),
        }

        if in_site_removed:
            names = ", ".join(f"{h.description} ({h.label}, "
                              f"{h.distance_to_ligand} A from the ligand)"
                              for h in in_site_removed)
            return ControlOutcome(
                name="Was anything important deleted?",
                verdict=Verdict.FAIL,
                plain_verdict=(
                    f"Stop. {names} would be removed, and it sits inside the "
                    "binding site. If it is part of how this protein works - "
                    "a catalytic metal, a cofactor - then docking without it "
                    "will still give you a score, and that score will be "
                    "meaningless."
                ),
                what_to_do=(
                    "Go back, tick this molecule in the 'keep' list, and run "
                    "this step again. Only leave it out if you know for "
                    "certain it does not belong there - for example a "
                    "sulphate ion that happened to land in the pocket during "
                    "crystallisation. When unsure, keep it."
                ),
                detail=detail,
                blocks_progress=True,
            )

        if important_removed:
            names = ", ".join(f"{h.description} ({h.label})"
                              for h in important_removed)
            return ControlOutcome(
                name="Was anything important deleted?",
                verdict=Verdict.WARN,
                plain_verdict=(
                    f"{names} was removed. It is away from the binding site, "
                    "so this is probably right - but check that it is not "
                    "structural, because some metals hold a fold together "
                    "rather than doing chemistry."
                ),
                what_to_do=(
                    "Look up what this molecule does in your protein - one "
                    "line in a textbook or UniProt is enough. If it is "
                    "catalytic or structural, go back and keep it. If it is "
                    "just crystallisation leftover, carry on."
                ),
                detail=detail,
                blocks_progress=True,
            )

        if kept_important:
            names = ", ".join(f"{h.description} ({h.label})"
                              for h in kept_important)
            return ControlOutcome(
                name="Was anything important deleted?",
                verdict=Verdict.PASS,
                plain_verdict=(
                    f"Nothing important was lost. {names} was found in or "
                    "near the binding site and has been kept, which is "
                    "almost certainly what you want."
                ),
                what_to_do=(
                    "Nothing to change - press Next. This is the mistake "
                    "that quietly ruins most student docking projects, and "
                    "it has been avoided here."
                ),
                detail=detail,
            )

        return ControlOutcome(
            name="Was anything important deleted?",
            verdict=Verdict.PASS,
            plain_verdict=(
                "Only water and crystallisation chemistry were removed. No "
                "metal or cofactor was found in this structure. Read the "
                "list below anyway - it is short, and this is the step "
                "where a quiet deletion costs you the whole result."
            ),
            what_to_do=(
                "Read the 'removed' list below once, then press Next. You "
                "are looking for anything that is not water, glycerol, or a "
                "simple salt. If everything in that list is laboratory "
                "leftover, you are fine."
            ),
            detail=detail,
        )

    @staticmethod
    def _version() -> str:
        try:
            import pdbfixer
            return getattr(pdbfixer, "__version__", "installed")
        except Exception:  # noqa: BLE001
            return "unknown"


_STANDARD = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "CYX", "SEC", "PYL", "MSE", "ASX", "GLX", "UNK",
    "A", "C", "G", "U", "I", "DA", "DC", "DG", "DT", "DI", "N",
}


def _is_standard_residue(name: str) -> bool:
    return name.strip().upper() in _STANDARD

"""
Stage 2 -- Obtain the structure.

Finds the best experimental structure for the target, or falls back to
the AlphaFold model, and puts a quality card in front of the user
before anything is docked into it.

The field paths below were read off real responses for EGFR (UniProt
P00533, PDB 1M17, AlphaFold AF-P00533-F1) captured by
verify_endpoints.py, not guessed.

What the ranking optimises for, in order:

  1. A bound ligand. A structure solved with something in the pocket
     is worth far more for docking than a higher-resolution empty one,
     because it proves the pocket exists in a bindable conformation
     and gives a pose to validate against. Most newcomers rank by
     resolution alone and this is the single most useful correction
     this stage makes.
  2. Resolution.
  3. Completeness -- unmodelled residues near the site are a problem.
"""

from __future__ import annotations

import json
import urllib.parse
from pathlib import Path
from typing import Any

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
from core.mmcif import CifError, cif_to_pdb, dropped_chains
from core.web import Web

PDB_DATA = "https://data.rcsb.org/rest/v1/core/entry"
PDB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
PDB_FILES = "https://files.rcsb.org/download"

# Below this many residues of the target protein, a structure is not a
# structure OF your protein -- it is a structure of something else with
# a piece of your protein stuck to it.
#
# Found the hard way. A search for SV2A returned 4V11, "Structure of
# Synaptotagmin-1 with SV2A peptide": 1.95 angstroms, no missing
# residues, calcium ions the metadata counted as ligands. It ranked top
# and the quality card called it "good for this purpose". It is a
# structure of a different protein.
#
# The test is absolute, not proportional: 1M17 holds only 27% of
# full-length EGFR and is a perfectly good kinase domain. A domain is
# hundreds of residues; a peptide is tens.
MIN_TARGET_RESIDUES = 50
ALPHAFOLD = "https://alphafold.ebi.ac.uk/api/prediction"


def _dig(blob: Any, *path: Any, default: Any = None) -> Any:
    cur = blob
    for key in path:
        try:
            cur = cur[key]
        except (KeyError, IndexError, TypeError):
            return default
    return cur


# ----------------------------------------------------------------------
# Parsing -- paths confirmed against captured responses
# ----------------------------------------------------------------------

def parse_pdb_entry(blob: dict) -> dict[str, Any]:
    res = _dig(blob, "rcsb_entry_info", "resolution_combined", 0)
    unmodelled = _dig(
        blob, "rcsb_entry_info", "deposited_unmodeled_polymer_monomer_count",
        default=0,
    ) or 0
    modelled = _dig(
        blob, "rcsb_entry_info", "deposited_modeled_polymer_monomer_count",
        default=0,
    ) or 0

    return {
        "source": "experimental",
        "id": blob.get("rcsb_id", ""),
        "title": _dig(blob, "struct", "title", default="") or "",
        "method": _dig(blob, "rcsb_entry_info", "experimental_method",
                       default="") or "",
        "resolution": float(res) if res is not None else None,
        "released": (_dig(blob, "rcsb_accession_info",
                          "initial_release_date", default="") or "")[:10],
        "r_free": _dig(blob, "refine", 0, "ls_R_factor_R_free"),
        "r_work": _dig(blob, "refine", 0, "ls_R_factor_R_work"),
        "clashscore": _dig(blob, "pdbx_vrpt_summary_geometry", 0, "clashscore"),
        "n_ligands": _dig(blob, "rcsb_entry_info", "nonpolymer_entity_count",
                          default=0) or 0,
        "ligand_atoms": _dig(blob, "refine_hist", 0,
                             "pdbx_number_atoms_ligand", default=0) or 0,
        "residues_modelled": modelled,
        "residues_missing": unmodelled,
        "protein_chains": _dig(blob, "rcsb_entry_info",
                               "polymer_entity_count_protein", default=0) or 0,
    }


def parse_alphafold(blob: list) -> dict[str, Any]:
    """AlphaFold returns a list; the first element is the current model."""
    e = blob[0] if isinstance(blob, list) and blob else {}
    return {
        "source": "predicted",
        "id": e.get("entryId", ""),
        "title": f"AlphaFold model of {e.get('uniprotDescription', '')}".strip(),
        "method": e.get("toolUsed", "AlphaFold"),
        "mean_plddt": e.get("globalMetricValue"),
        "version": e.get("latestVersion"),
        "created": (e.get("modelCreatedDate") or "")[:10],
        # The pLDDT distribution, which is what actually matters and is
        # given here without needing to download the structure at all.
        "frac_very_high": e.get("fractionPlddtVeryHigh"),
        "frac_confident": e.get("fractionPlddtConfident"),
        "frac_low": e.get("fractionPlddtLow"),
        "frac_very_low": e.get("fractionPlddtVeryLow"),
        "length": e.get("uniprotEnd"),
        "pdb_url": e.get("pdbUrl", ""),
        "cif_url": e.get("cifUrl", ""),
        "pae_url": e.get("paeDocUrl", ""),
        "plddt_url": e.get("plddtDocUrl", ""),
    }



def _spread(items: list[str], count: int) -> list[str]:
    """
    Take `count` items spread evenly across a list, not the first few.

    Used to choose which structures to look at. Taking the head of a
    list sorted by coverage returns one kind of structure -- the
    biggest -- and for docking that is the wrong kind. An even spread
    returns large, medium and small, and lets the scoring decide.
    """
    if len(items) <= count or count <= 0:
        return list(items)
    step = len(items) / count
    picked = [items[min(int(i * step), len(items) - 1)] for i in range(count)]
    out: list[str] = []
    for item in picked:
        if item not in out:
            out.append(item)
    return out


def real_ligands(pdb_text: str) -> dict[str, Any]:
    """
    What is actually bound in this file -- a drug, or only sugar?

    The Protein Data Bank's own ligand count is
    `nonpolymer_entity_count`, and it counts everything that is not
    protein: glycerol, sulfate, a chloride ion, and above all the
    N-acetylglucosamine that decorates every glycoprotein. A
    heavily glycosylated extracellular domain with nothing in its
    binding site reports a cheerful handful of "ligands".

    So step 2 would rank such a structure as having a drug bound,
    announce "with a bound ligand", and hand it to step 3 -- which
    applies the real classification, finds no drug, and reports none.
    The two steps disagreed about what a ligand is, and the user was
    left with a binding site placed by guesswork and no way to see
    why. Recounting from the coordinates with the same rules step 3
    uses is what stops that.
    """
    from stages.receptor import KIND_COFACTOR, KIND_LIGAND, classify

    groups: dict[str, int] = {}
    kinds: dict[str, str] = {}
    for line in pdb_text.splitlines():
        if not line.startswith("HETATM"):
            continue
        if (line[76:78].strip().upper() or "C") == "H":
            continue
        resname = line[17:20].strip().upper()
        label = f"{resname} {line[21]}{line[22:27].strip()}"
        groups[label] = groups.get(label, 0) + 1
        kinds[label] = classify(resname)[0]

    drugs = {label: n for label, n in groups.items()
             if kinds[label] in (KIND_LIGAND, KIND_COFACTOR) and n >= 10}
    other = sorted({label.split()[0] for label in groups
                    if label not in drugs})
    best = max(drugs, key=lambda k: drugs[k]) if drugs else ""
    return {
        "n_ligands": len(drugs),
        "ligand_atoms": drugs.get(best, 0),
        "best_ligand": best,
        "non_drug_components": other[:12],
    }


def score_candidate(entry: dict) -> float:
    """
    Rank structures for docking. Higher is better.

    A bound ligand is worth more than a whole Angstrom of resolution,
    which is the opposite of how most people choose.
    """
    score = 0.0
    if entry.get("n_ligands", 0) > 0 and entry.get("ligand_atoms", 0) >= 10:
        score += 10.0           # a real drug-sized ligand, not a buffer ion
    elif entry.get("n_ligands", 0) > 0:
        score += 2.0            # something bound, possibly just a sulfate

    # Resolution, with a cliff rather than a slope.
    #
    # The old version was a straight line that still awarded points at
    # 4.5 angstroms, so a 4.53 A whole-receptor reconstruction scored
    # positively and won. At that resolution you cannot see side
    # chains at all -- the backbone is a tube and the residues that
    # line the pocket are modelled, not observed. Docking into it is
    # docking into somebody's guess. Under 2.5 is sharp, 3.0 is the
    # usual working limit, and past 3.5 the structure is not a
    # docking target whatever else is good about it.
    res = entry.get("resolution")
    if res:
        res = float(res)
        if res <= 3.0:
            score += (3.0 - res) * 4.0 + 6.0
        elif res <= 3.5:
            score += (3.5 - res) * 8.0
        else:
            score -= (res - 3.5) * 12.0

    missing = entry.get("residues_missing", 0)
    modelled = entry.get("residues_modelled", 1) or 1
    score -= (missing / max(modelled, 1)) * 8.0

    if (entry.get("method") or "").upper().startswith("X-RAY"):
        score += 1.0

    # A structure that barely contains the target is not a candidate at
    # all, whatever else is good about it.
    # Coverage is a floor, not a prize. Enough of the protein to hold
    # a pocket is required; more than that is not better, and
    # rewarding it is what put whole-receptor reconstructions at the
    # top of the list.
    covered = entry.get("target_residues")
    if covered is not None:
        if covered < MIN_TARGET_RESIDUES:
            score -= 1000.0
        elif covered < 100:
            score -= 5.0
    return round(score, 2)


# ----------------------------------------------------------------------
# The stage
# ----------------------------------------------------------------------

class StructureStage(Stage):
    stage_id = "02_structure"
    title = "Obtain the structure"
    tier = Tier.API
    requires = ["01_target"]

    def __init__(self, project, web: Web | None = None) -> None:
        super().__init__(project)
        self.web = web or Web(project.folder / "cache")

    # -- contract ---------------------------------------------------------

    def check(self) -> CheckResult:
        if not self.project.is_complete("01_target"):
            return CheckResult(
                readiness=Readiness.NEEDS_INPUT,
                reason="No target chosen yet.",
                fix="Complete the target step first.",
            )
        return CheckResult(readiness=Readiness.READY)

    def explain(self) -> Teaching:
        return Teaching(
            button="Find the best available 3D structure of your target.",
            panel=(
                "Docking needs coordinates. They come from one of two "
                "places and the difference matters.\n\n"
                "An experimental structure is a measurement - usually X-ray "
                "crystallography, sometimes cryo-electron microscopy. "
                "Resolution tells you how fine the detail is; below about "
                "2.5 angstroms you can see where side chains point, which "
                "is what docking depends on. A predicted structure, from "
                "AlphaFold, is a very good guess. For a well-studied "
                "protein fold it is often close to right; for a flexible "
                "loop or a disordered region it is not, and it tells you so "
                "through its confidence score.\n\n"
                "The one thing worth unlearning: most people pick the "
                "structure with the best resolution. For docking, a "
                "structure solved with a drug already bound in the pocket "
                "is usually the better choice even at a worse resolution. "
                "It shows a ligand-bound shape and supplies a reference "
                "pose for a docking check. The software can compare a "
                "matching docked ligand with that reference; this tests "
                "pose reproduction, not binding affinity."
            ),
            simple=(
                "To dock a drug into a protein, the computer needs to "
                "have a model of the protein's shape - estimated "
                "positions for its atoms.\n\n"
                "There are two ways to get this shape:\n\n"
                "1. MEASURED - scientists grew a crystal of the protein "
                "and shone X-rays through it. This is a real "
                "measurement.\n\n"
                "2. PREDICTED - a computer program called AlphaFold "
                "guessed the shape. Often very good, but still a "
                "guess.\n\n"
                "The software picks the best one for you. It prefers a "
                "structure that already has a drug sitting inside it, "
                "even if the picture is slightly blurry. Why? Because "
                "it shows a measured ligand-bound shape. It also gives "
                "us a reference pose to test ourselves later. Check "
                "that the binding site is suitable and resolved.\n\n"
                "You usually do not need to change anything here. Just "
                "press Run."
            ),
            not_mean=[
                "A high confidence score from AlphaFold means the fold is "
                "probably right. It does not mean the side chains in your "
                "pocket are in the position a drug would find, and it says "
                "nothing about which of several shapes the protein adopts.",
                "A crystal structure is one snapshot of one form, frozen, "
                "often at a non-physiological pH, and sometimes distorted "
                "by how the crystal packs. It is not 'the' structure of "
                "the protein.",
                "Good resolution does not mean your binding site is well "
                "resolved. Flexible regions are often missing entirely from "
                "an otherwise excellent structure, and if the gap is at "
                "your pocket, it may be unsuitable without careful "
                "assessment or rebuilding.",
                "For a multi-domain AlphaFold model, inspect local "
                "confidence and uncertainty in domain placement. A "
                "confident fold does not establish a ligand-ready pocket.",
            ],
            learn_more=[
                ("How to read AlphaFold confidence",
                 "https://alphafold.ebi.ac.uk/faq"),
                ("RCSB PDB - structure quality",
                 "https://www.rcsb.org/docs/general-help/assessing-structure-quality"),
            ],
        )

    # -- finding structures -----------------------------------------------

    def search_pdb(self, accession: str, limit: int = 10) -> list[str]:
        query = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers"
                                 ".reference_sequence_identifiers"
                                 ".database_accession",
                    "operator": "exact_match",
                    "value": accession,
                },
            },
            "return_type": "entry",
            "request_options": {"paginate": {"start": 0, "rows": limit}},
        }
        url = f"{PDB_SEARCH}?json={urllib.parse.quote(json.dumps(query))}"
        blob, res = self.web.get_json(url, service="the Protein Data Bank")
        if not res.ok or not blob:
            return []
        return [
            hit.get("identifier", "")
            for hit in (_dig(blob, "result_set", default=[]) or [])
            if hit.get("identifier")
        ]

    def describe_pdb(self, pdb_id: str) -> dict | None:
        blob, res = self.web.get_json(
            f"{PDB_DATA}/{pdb_id}", service="the Protein Data Bank"
        )
        if not res.ok or not blob:
            return None
        return parse_pdb_entry(blob)

    def alphafold(self, accession: str) -> dict | None:
        blob, res = self.web.get_json(
            f"{ALPHAFOLD}/{accession}", service="the AlphaFold Database"
        )
        if not res.ok or not isinstance(blob, list) or not blob:
            return None
        return parse_alphafold(blob)

    # -- run ---------------------------------------------------------------

    def run(self, prefer: str = "auto", max_candidates: int = 20,
            pdb_id: str = "", **kwargs: Any) -> StageResult:
        target_step = self.project.get_step("01_target")
        if not target_step:
            return StageResult(success=False,
                               message="Choose a target first.")
        accession = _dig(json.loads(target_step["values_json"] or "{}"),
                         "target", "accession", default="")
        if not accession:
            return StageResult(
                success=False,
                message="The target step did not record a UniProt accession.",
            )

        target = json.loads(target_step["values_json"] or "{}").get("target", {})
        coverage = target.get("pdb_coverage", {}) or {}
        protein_length = target.get("length", 0) or 0

        candidates: list[dict] = []
        if prefer in ("auto", "experimental"):
            # Prefer the structures UniProt says contain most of this
            # protein. Without this the search returns anything the
            # accession is cross-referenced from, including complexes
            # where the target is a short bound peptide.
            # Which structures to look at, and this is where the
            # previous version went wrong.
            #
            # It sorted by coverage and took the top few, on the
            # reasoning that a structure containing more of the
            # protein is a better structure. For docking that is
            # backwards. The entries with the most residues of a
            # receptor like EGFR are whole-receptor cryo-EM
            # reconstructions: four or five angstroms, no inhibitor,
            # enormous. The kinase domain with erlotinib in it -- the
            # structure anybody would actually dock into -- has a
            # third of the coverage and never reached the list.
            #
            # Coverage belongs here as a floor, not as a ranking. Above
            # the floor it stops mattering, so the selection spreads
            # across the whole qualifying range rather than crowding at
            # the top, and the real criteria (a bound drug, then
            # resolution) decide between them afterwards.
            qualifying = [
                pid for pid, c in sorted(
                    coverage.items(),
                    key=lambda kv: kv[1].get("residues", 0), reverse=True)
                if (c.get("residues") or 0) >= MIN_TARGET_RESIDUES
            ]
            ranked_ids = _spread(qualifying, max_candidates)

            # The database's own relevance order as well, so a recent
            # or heavily cited entry is not missed because UniProt
            # happens to list it late.
            for pid in self.search_pdb(accession, max_candidates // 2):
                if pid not in ranked_ids:
                    ranked_ids.append(pid)

            ids = ranked_ids or self.search_pdb(accession, max_candidates)

            for candidate_id in ids:
                entry = self.describe_pdb(candidate_id)
                if not entry:
                    continue
                cov = coverage.get(candidate_id)
                if cov:
                    entry["target_residues"] = cov["residues"]
                    entry["target_chains"] = cov["chains"]
                    entry["target_range"] = f"{cov['start']}-{cov['end']}"
                    if protein_length:
                        entry["target_fraction"] = round(
                            cov["residues"] / protein_length, 3)
                entry["score"] = score_candidate(entry)
                candidates.append(entry)
            candidates.sort(key=lambda e: e["score"], reverse=True)

        model = None
        if prefer in ("auto", "predicted") or not candidates:
            model = self.alphafold(accession)

        # The user naming a structure overrides the ranking. They have
        # a reason the software cannot see -- a paper everyone in the
        # field cites, a structure their guide asked for, a mutant --
        # and the ranking here is a default, not a verdict.
        wanted = pdb_id.strip().upper()
        if wanted:
            match = next((e for e in candidates
                          if str(e.get("id", "")).upper() == wanted), None)
            if match is None:
                # Not among the ones examined, which is routine for a
                # protein with hundreds of entries. Ask the database
                # for it directly rather than refusing: a structure
                # the user names on purpose should always be reachable.
                direct = self.describe_pdb(wanted)
                if direct:
                    cov = coverage.get(wanted)
                    if cov:
                        direct["target_residues"] = cov["residues"]
                        direct["target_chains"] = cov["chains"]
                        direct["target_range"] = f"{cov['start']}-{cov['end']}"
                    direct["score"] = score_candidate(direct)
                    candidates.insert(0, direct)
                    match = direct
            if match is None:
                listed = ", ".join(str(e.get("id", "")) for e in
                                   candidates[:8]) or "none"
                return StageResult(
                    success=False,
                    message=(
                        f"{wanted} is not one of the structures listed for "
                        "this protein, so the software has no description "
                        "of it to check.\n\n"
                        f"The ones it found are: {listed}.\n\n"
                        "Two usual reasons. The code may belong to a "
                        "different protein - check it on rcsb.org. Or it "
                        "is a very recent release that the search has not "
                        "indexed yet, in which case leave this box empty "
                        "and use the best of the list."
                    ),
                )
            candidates = [match] + [e for e in candidates if e is not match]

        chosen = candidates[0] if candidates else model
        if chosen is None:
            return StageResult(
                success=False,
                message=(
                    "No experimental structure and no AlphaFold model was "
                    "found for this protein. This is unusual and normally "
                    "means the accession is wrong."
                ),
            )

        # Fetch the coordinates. Recording which structure to use and
        # not fetching it leaves the next step with nothing to open.
        local, why = self._download(chosen)
        if local is None:
            return StageResult(
                success=False,
                message=f"Found {chosen['id']} but could not get the "
                        f"coordinates. {why}",
            )
        chosen["local_file"] = str(local)

        # The database's ligand count is whatever is not protein. Now
        # that the file is here, count what is actually bound using the
        # same rules the next step will apply, so the two cannot
        # disagree.
        if chosen["source"] == "experimental":
            claimed = chosen.get("n_ligands", 0)
            found = real_ligands(local.read_text(errors="ignore"))
            chosen.update(found)
            if claimed and not found["n_ligands"]:
                chosen["ligand_warning"] = (
                    f"The database lists {claimed} non-protein "
                    "component(s) for this entry, but none of them is a "
                    "drug-like molecule - they are "
                    + (", ".join(found["non_drug_components"][:6])
                       or "sugars, ions or buffer")
                    + ". There is nothing bound in the binding site."
                )
        if chosen.get("converted_from_cif"):
            chosen["note"] = (
                "The Protein Data Bank has no old-format file for this "
                "entry, so the modern mmCIF file was downloaded and "
                "converted here. Nothing is lost."
            )

        self.record(
            ProvenanceEntry(
                stage=self.stage_id,
                tool="RCSB PDB" if chosen["source"] == "experimental"
                     else "AlphaFold Database",
                version=str(chosen.get("version", "current")),
                parameters={
                    "accession": accession,
                    "chosen": chosen["id"],
                    "candidates_considered": len(candidates),
                    "ranked_by": "bound ligand, then resolution, "
                                 "then completeness",
                },
                source_release=chosen.get("released")
                               or chosen.get("created", ""),
                citation=(
                    "Berman HM et al. The Protein Data Bank. Nucleic Acids "
                    "Res. 2000;28:235-242."
                    if chosen["source"] == "experimental" else
                    "Varadi M et al. AlphaFold Protein Structure Database. "
                    "Nucleic Acids Res. 2024;52:D368-D375."
                ),
            )
        )
        for e in self.provenance:
            self.project.add_provenance(e)

        return StageResult(
            success=True,
            summary=self._summary(chosen),
            values={
                "chosen": chosen,
                "alternatives": candidates[1:5],
                "alphafold_available": model is not None,
            },
            control=self._quality_card(chosen),
        )

    def _download(self, entry: dict) -> tuple[Path | None, str]:
        """
        Save the coordinates into the project folder.

        Returns the file and, when there isn't one, a sentence saying
        why in terms the user can act on.

        Two things are going on here that an earlier version got
        wrong. First, a failure to download is not the same as a file
        that does not exist, and telling someone to check their
        internet connection when the Protein Data Bank has answered
        perfectly well and said "no such file" sends them to look in
        the wrong place for an hour.

        Second, for a growing number of entries -- recent cryo-EM
        structures above all -- there genuinely is no old-format .pdb
        file. mmCIF is now the archive's master format and the legacy
        copy is only generated where it fits. So when .pdb is missing,
        this fetches the mmCIF and converts it, rather than treating a
        perfectly good structure as unobtainable.
        """
        folder = self.project.step_folder(self.stage_id)

        if entry["source"] != "experimental":
            url = entry.get("pdb_url", "")
            if not url:
                return None, ("The AlphaFold Database did not give a "
                              "download address for this model.")
            got = self.web.get(url, service="the AlphaFold Database",
                               timeout=120)
            if not got.ok or not got.text:
                return None, self._why(got)
            path = folder / f"{entry['id']}.pdb"
            path.write_text(got.text)
            return path, ""

        # The old format first: it needs no conversion, so nothing can
        # go wrong in translation.
        legacy = self.web.get(f"{PDB_FILES}/{entry['id']}.pdb",
                              service="the Protein Data Bank", timeout=120)
        if legacy.ok and legacy.text and "ATOM  " in legacy.text:
            path = folder / f"{entry['id']}.pdb"
            path.write_text(legacy.text)
            return path, ""

        # A missing file is worth a second attempt in the other
        # format. An unreachable server is not: the mmCIF lives at the
        # same address, so trying again only doubles the wait before
        # the same answer.
        if legacy.status == 0:
            return None, self._why(legacy)

        cif = self.web.get(f"{PDB_FILES}/{entry['id']}.cif",
                           service="the Protein Data Bank", timeout=180)
        if not cif.ok or not cif.text:
            if legacy.status == 404 and cif.status == 404:
                return None, (
                    f"The Protein Data Bank has no coordinate file for "
                    f"{entry['id']} in either format, which is unusual - "
                    "the entry may have been withdrawn. Go back and pick "
                    "a different structure from the list."
                )
            return None, self._why(cif)

        try:
            converted = cif_to_pdb(cif.text)
        except CifError as exc:
            return None, str(exc)
        except Exception as exc:                       # noqa: BLE001
            return None, (
                "The modern-format file for this entry downloaded but "
                f"could not be read ({exc}). Go back and pick a "
                "different structure."
            )

        path = folder / f"{entry['id']}.pdb"
        path.write_text(converted)
        entry["converted_from_cif"] = True
        left_out = dropped_chains(cif.text)
        if left_out:
            entry["chains_dropped"] = left_out
        return path, ""

    @staticmethod
    def _why(got: Any) -> str:
        """
        The real reason, in the user's words, never a guess.

        `Fetched` already knows whether the server answered and what it
        said. Throwing that away and substituting "check your internet
        connection" is how a working connection gets blamed for a
        missing file.
        """
        if got.problem:
            return f"{got.problem} {got.suggestion}".strip()
        if got.status:
            return (f"The Protein Data Bank answered with an error "
                    f"(HTTP {got.status}). Try this step again in a "
                    "minute; if it keeps happening, pick a different "
                    "structure.")
        return ("The Protein Data Bank could not be reached. Check your "
                "internet connection - and if you are on a college or "
                "hospital network, note that some of them block file "
                "downloads from outside sites.")

    @staticmethod
    def _summary(e: dict) -> str:
        if e["source"] == "experimental":
            res = f"{e['resolution']} A" if e.get("resolution") else "no resolution"
            lig = "with a bound ligand" if e.get("n_ligands") else "no ligand bound"
            return f"{e['id']} - {e['method']}, {res}, {lig}"
        return (f"{e['id']} - AlphaFold model v{e.get('version')}, "
                f"mean confidence {e.get('mean_plddt')}")

    # -- the control -------------------------------------------------------

    @staticmethod
    def _quality_card(e: dict) -> ControlOutcome:
        """
        Specification Stage 2: a mandatory quality card, and progress is
        blocked when the pocket region may not be resolved.

        Written for a reader who has never seen a PDB entry. Every
        verdict says what was found AND what to do about it, because a
        finding with no action is only a worry.
        """
        if e["source"] == "predicted":
            very_low = e.get("frac_very_low") or 0.0
            very_high = e.get("frac_very_high") or 0.0
            detail = {
                "mean_pLDDT": e.get("mean_plddt"),
                "very high confidence": f"{very_high:.0%}",
                "confident": f"{(e.get('frac_confident') or 0):.0%}",
                "low": f"{(e.get('frac_low') or 0):.0%}",
                "very low confidence": f"{very_low:.0%}",
                "residues": e.get("length"),
            }
            if very_low >= 0.20:
                return ControlOutcome(
                    name="Can we dock into this structure?",
                    verdict=Verdict.WARN,
                    plain_verdict=(
                        "This shape was predicted by AlphaFold, not "
                        f"measured. AlphaFold itself says {very_low:.0%} of "
                        "it is a guess it does not trust. Those parts are "
                        "floppy loops that have no fixed shape in the real "
                        "protein."
                    ),
                    what_to_do=(
                        "Carry on, but at the Prepare Protein step keep only "
                        "the well-predicted part. If your pocket turns out "
                        "to sit inside the untrusted region, stop and look "
                        "for an experimental structure instead - docking "
                        "into a guessed shape gives a confident score for "
                        "nothing."
                    ),
                    detail=detail,
                    blocks_progress=True,
                )
            if very_high < 0.50:
                return ControlOutcome(
                    name="Can we dock into this structure?",
                    verdict=Verdict.WARN,
                    plain_verdict=(
                        "This shape was predicted, not measured. AlphaFold "
                        f"is fully confident about {very_high:.0%} of it. "
                        "That is on the low side."
                    ),
                    what_to_do=(
                        "You can continue. At the Binding Site step, check "
                        "that your pocket sits in a confident region - the "
                        "software will show you. Treat the scores as a "
                        "ranking, not as numbers."
                    ),
                    detail=detail,
                )
            return ControlOutcome(
                name="Can we dock into this structure?",
                verdict=Verdict.PASS,
                plain_verdict=(
                    "This shape was predicted by AlphaFold, not measured in "
                    f"a laboratory - but AlphaFold is fully confident about "
                    f"{very_high:.0%} of it, which is good."
                ),
                what_to_do=(
                    "Go ahead. One thing to remember: a predicted structure "
                    "has no drug bound in it, so the self-test at the "
                    "docking step cannot run. You will have no way of "
                    "proving your docking setup works."
                ),
                detail=detail,
            )

        # Experimental
        covered = e.get("target_residues")
        if covered is not None and covered < MIN_TARGET_RESIDUES:
            return ControlOutcome(
                name="Can we dock into this structure?",
                verdict=Verdict.FAIL,
                plain_verdict=(
                    "This is not a structure of your protein. Only "
                    f"{covered} amino acids of it are in this file "
                    f"(numbers {e.get('target_range', '?')}). That is a tiny "
                    "fragment stuck to some other protein - not your protein "
                    "folded up. Your protein's name is listed on this entry, "
                    "which is why it came up, but it is not what the "
                    "scientists actually solved."
                ),
                what_to_do=(
                    "Do not use this. Go back and pick a different structure "
                    "from the list - look for one with a few hundred amino "
                    "acids of your protein, not a few dozen."
                ),
                detail={
                    "entry": e.get("id"),
                    "title": (e.get("title") or "")[:90],
                    "amino acids of your protein present": covered,
                    "which ones": e.get("target_range"),
                    "resolution": f"{e.get('resolution')} A",
                },
                blocks_progress=True,
            )

        res = e.get("resolution")
        missing = e.get("residues_missing", 0)
        modelled = e.get("residues_modelled", 0)
        total = missing + modelled
        frac_missing = (missing / total) if total else 0.0

        detail = {
            "amino acids of your protein present": covered
            if covered is not None else "not known",
            "resolution (lower = sharper)": f"{res} A" if res else "unknown",
            "method": e.get("method"),
            "R-free": e.get("r_free"),
            "R-work": e.get("r_work"),
            "clashscore": e.get("clashscore"),
            "amino acids visible": modelled,
            "amino acids missing": missing,
            "drug-like molecules bound inside": e.get("n_ligands"),
            "which one": e.get("best_ligand") or "none",
            "other components (not drugs)":
                ", ".join(e.get("non_drug_components", [])) or "none",
            "released": e.get("released"),
        }

        # Each problem: what was found, then what to do about it.
        problems = []
        actions = []
        if res and res > 3.0:
            problems.append(
                f"the picture is blurry ({res} angstroms - below 2.5 is "
                "sharp, above 3.0 is blurry), so the side chains that line "
                "the pocket may be in slightly the wrong place"
            )
            actions.append(
                "check the list for a sharper structure of the same protein "
                "before using this one"
            )
        if frac_missing > 0.05:
            problems.append(
                f"{missing} amino acids ({frac_missing:.0%}) are missing - "
                "they were moving too much to be photographed, so there are "
                "holes in the protein"
            )
            actions.append(
                "at the Prepare Protein step, see whether any missing piece "
                "is near your pocket; small holes can be rebuilt, a hole in "
                "the pocket itself means you should pick another structure"
            )
        if not e.get("n_ligands"):
            problems.append(
                "no drug-like bound ligand was identified by this "
                "workflow, so it has no reference pose for the docking "
                "test. Choose a site using other structural or "
                "experimental evidence"
                + (" - " + e["ligand_warning"].split(" - ", 1)[-1]
                   if e.get("ligand_warning") else "")
            )
            actions.append(
                "prefer a structure with a drug already in it if the list "
                "has one - otherwise the self-test at the docking step "
                "cannot run; record that limitation"
            )

        if frac_missing > 0.05 or (res and res > 3.0):
            return ControlOutcome(
                name="Can we dock into this structure?",
                verdict=Verdict.WARN,
                plain_verdict=(
                    "This is a real laboratory structure, but read this "
                    "before going on: " + "; and ".join(problems) + "."
                ),
                what_to_do=" Then ".join(
                    a[0].upper() + a[1:] for a in actions
                ) + ".",
                detail=detail,
                blocks_progress=True,
            )

        if problems:
            return ControlOutcome(
                name="Can we dock into this structure?",
                verdict=Verdict.WARN,
                plain_verdict=(
                    "Usable, with one point to note: "
                    + "; and ".join(problems) + "."
                ),
                what_to_do=" Then ".join(
                    a[0].upper() + a[1:] for a in actions
                ) + ".",
                detail=detail,
            )

        return ControlOutcome(
            name="Can we dock into this structure?",
            verdict=Verdict.PASS,
            plain_verdict=(
                "This structure passes the initial screening criteria: "
                "it is a real laboratory measurement, not a computer "
                f"prediction; the picture is sharp enough at {res} "
                "angstroms; a drug-like bound ligand was identified; "
                "and few residues are reported missing. These global "
                "checks do not guarantee local pocket quality."
            ),
            what_to_do=(
                "Inspect the chosen site and then continue to preparation. "
                "Include the bound ligand among your docking molecules "
                "to compare its predicted and reference poses. Passing "
                "that test supports reference-pose reproduction; it "
                "does not validate affinity predictions."
            ),
            detail=detail,
        )

"""
Stage 1 -- Define the target.

Takes a gene symbol, protein name or UniProt accession and builds a
dossier: what the protein is, what it does, whether anyone has
successfully drugged it, and what structures exist.

The control reports the evidence returned for small-molecule
tractability. Missing evidence does not establish the absence of a
binding pocket or previous research. Specification Stage 1.

NOTE ON VERIFICATION
    The endpoint addresses and field paths below could not be tested
    where this was written -- that machine has no route to UniProt,
    ChEMBL or Open Targets. Run verify_endpoints.py on a connected
    machine before trusting this stage. Parsing is written to tolerate
    missing fields rather than crash, so a changed field shows up as a
    blank in the dossier rather than an error.
"""

from __future__ import annotations

import json
import urllib.parse
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
from core.web import Web

UNIPROT = "https://rest.uniprot.org/uniprotkb"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
OPENTARGETS = "https://api.platform.opentargets.org/api/v4/graphql"


# ----------------------------------------------------------------------
# Small, defensive readers -- a missing field is a blank, never a crash
# ----------------------------------------------------------------------

def _dig(blob: Any, *path: str | int, default: Any = None) -> Any:
    cur = blob
    for key in path:
        try:
            cur = cur[key]
        except (KeyError, IndexError, TypeError):
            return default
    return cur


def parse_uniprot(blob: dict) -> dict[str, Any]:
    """Pull the handful of fields the dossier shows."""
    genes = _dig(blob, "genes", default=[]) or []
    names = [_dig(g, "geneName", "value") for g in genes]

    function = ""
    for comment in _dig(blob, "comments", default=[]) or []:
        if comment.get("commentType") == "FUNCTION":
            function = _dig(comment, "texts", 0, "value", default="") or ""
            break

    xrefs = _dig(blob, "uniProtKBCrossReferences", default=[]) or []
    pdb_ids = [x.get("id", "") for x in xrefs if x.get("database") == "PDB"]

    # How much of THIS protein each structure actually contains.
    #
    # A UniProt entry cross-references every structure the protein
    # appears in, including ones where it appears as a ten-residue
    # peptide bound to something else entirely. Those are real
    # cross-references and completely wrong to dock into, and nothing
    # about resolution or completeness reveals it.
    #
    # UniProt records the range: "Chains: A=695-1022". Reading it here
    # costs nothing, because the entry has already been fetched.
    coverage: dict[str, dict[str, Any]] = {}
    for x in xrefs:
        if x.get("database") != "PDB":
            continue
        # Method and resolution sit beside the chain range, so the
        # structure step can rank candidates before fetching any of
        # them: "Method: X-ray", "Resolution: 2.60 A".
        props = {str(q.get("key")): str(q.get("value", ""))
                 for q in x.get("properties", []) or []}
        try:
            resolution: float | None = float(
                props.get("Resolution", "").split()[0])
        except (ValueError, IndexError):
            resolution = None
        for prop in x.get("properties", []) or []:
            if prop.get("key") != "Chains":
                continue
            spec = str(prop.get("value", ""))
            if "=" not in spec:
                continue
            chains, _, span = spec.partition("=")
            try:
                lo, _, hi = span.partition("-")
                start, end = int(lo), int(hi)
            except ValueError:
                continue
            coverage[x.get("id", "")] = {
                "chains": [c.strip() for c in chains.split("/") if c.strip()],
                "start": start,
                "end": end,
                "residues": max(0, end - start + 1),
                "method": props.get("Method", ""),
                "resolution": resolution,
            }

    # The Ensembl gene id is needed to ask Open Targets about this
    # protein, and it is already here -- no second lookup required.
    # It arrives versioned ("ENSG00000146648.22"); the version is dropped.
    ensembl_genes: set[str] = set()
    for x in xrefs:
        if x.get("database") != "Ensembl":
            continue
        for prop in x.get("properties", []) or []:
            if prop.get("key") == "GeneId":
                ensembl_genes.add(str(prop.get("value", "")).split(".")[0])

    return {
        "accession": blob.get("primaryAccession", ""),
        "entry_name": blob.get("uniProtkbId", ""),
        "protein": _dig(blob, "proteinDescription", "recommendedName",
                        "fullName", "value", default="") or "",
        "genes": [n for n in names if n],
        "organism": _dig(blob, "organism", "scientificName", default="") or "",
        "length": _dig(blob, "sequence", "length", default=0) or 0,
        "sequence": _dig(blob, "sequence", "value", default="") or "",
        "function": function,
        "pdb_ids": pdb_ids,
        "pdb_coverage": coverage,
        "ensembl_gene": sorted(ensembl_genes)[0] if ensembl_genes else "",
        "reviewed": blob.get("entryType", "").startswith("UniProtKB reviewed"),
    }


# ----------------------------------------------------------------------
# Open Targets tractability -- the source for the mandatory control
#
# ChEMBL was the original choice and has been demoted. Across three
# runs from two machines it timed out twice and returned HTTP 500 on
# its own /status endpoint. A control that specification Section 9
# says cannot be switched off must not sit on a service that behaves
# like that.
#
# Open Targets is better anyway. Instead of a raw count of molecules,
# it gives a curated ladder of evidence, strongest first. Shape
# confirmed against a captured response: a flat list of
# {label, modality, value}, where modality "SM" is small molecules.
# ----------------------------------------------------------------------

# Ordered strongest to weakest. The order is the whole point.
SM_LADDER = [
    "Approved Drug",          # a marketed drug binds this protein
    "Advanced Clinical",      # reached phase 2 or 3
    "Phase 1 Clinical",       # reached humans
    "Structure with Ligand",  # a structure exists with something bound
    "High-Quality Ligand",    # a good chemical probe exists
    "High-Quality Pocket",    # a pocket worth binding has been found
    "Med-Quality Pocket",     # a pocket of some sort
    "Druggable Family",       # only: relatives of this protein are drugged
]
CLINICAL = set(SM_LADDER[:3])
STRUCTURAL = set(SM_LADDER[3:6])
INFERRED = set(SM_LADDER[6:])

TRACTABILITY_QUERY = (
    "query T($id: String!) { target(ensemblId: $id) { "
    "id approvedSymbol approvedName biotype "
    "tractability { label modality value } } }"
)


def parse_tractability(blob: dict) -> dict[str, Any]:
    """Reduce the 28-row tractability table to what the control needs."""
    target = _dig(blob, "data", "target", default=None)
    if not target:
        return {"available": False}

    rows = target.get("tractability") or []
    met = [r.get("label", "") for r in rows
           if r.get("modality") == "SM" and r.get("value")]

    return {
        "available": True,
        "symbol": target.get("approvedSymbol", ""),
        "name": target.get("approvedName", ""),
        "biotype": target.get("biotype", ""),
        "sm_buckets": [b for b in SM_LADDER if b in met],   # keep the order
        "clinical": sorted(CLINICAL & set(met), key=SM_LADDER.index),
        "structural": sorted(STRUCTURAL & set(met), key=SM_LADDER.index),
        "inferred_only": bool(met) and not (CLINICAL | STRUCTURAL) & set(met),
        "other_modalities": sorted({
            r.get("modality", "") for r in rows
            if r.get("value") and r.get("modality") != "SM"
        }),
    }


# ----------------------------------------------------------------------
# The stage
# ----------------------------------------------------------------------

class TargetStage(Stage):
    stage_id = "01_target"
    title = "Define the target"
    tier = Tier.API

    def __init__(self, project, web: Web | None = None) -> None:
        super().__init__(project)
        self.web = web or Web(project.folder / "cache")

    # -- contract ---------------------------------------------------------

    def check(self) -> CheckResult:
        probe = self.web.get(
            f"{UNIPROT}/P00533.json?fields=accession",
            service="UniProt", timeout=8, retries=0,
        )
        if probe.ok:
            return CheckResult(readiness=Readiness.READY)
        if probe.from_cache:
            return CheckResult(
                readiness=Readiness.DEGRADED,
                reason="Working from saved copies - UniProt is not answering.",
                fix="Targets you have looked up before will still work.",
            )
        return CheckResult(
            readiness=Readiness.DEGRADED,
            reason=probe.problem or "UniProt could not be reached.",
            fix=(probe.suggestion or "")
            + " You can still open a worked example offline.",
        )

    def explain(self) -> Teaching:
        return Teaching(
            rule=(
                "Name the protein your drug acts on, using its gene symbol - then read the druggability verdict before going on."
            ),
            rule_detail=[
                "Type the short gene name (EGFR, SV2A, PTGS2), not the drug name and not the disease.",
                'If the verdict says only "Druggable Family" or nothing at all, this search has not established direct evidence for the protein. Check the literature and the type of evidence before continuing.',
                "An approved drug provides useful evidence that a target can be modulated. For small-molecule docking, also check its binding site and mechanism.",
            ],
            button="Look up the protein and find out whether it can be drugged.",
            panel=(
                "Everything downstream depends on this choice, so it is worth "
                "more thought than it usually gets.\n\n"
                "A gene being associated with a disease does not make its "
                "protein a drug target. Three further things have to be true. "
                "There should be evidence that changing its activity can "
                "help the disease. For this small-molecule exercise, we "
                "also need evidence for a suitable binding site. The "
                "intended intervention must have an acceptable balance of "
                "benefit and harm. The absence of a drug does not tell us "
                "which of these questions remains unresolved.\n\n"
                "This step gathers what is known: what the protein does, "
                "where it sits, which structures exist, and crucially what "
                "has already been tried against it. That last part is the "
                "most useful and the most often skipped."
            ),
            simple=(
                "Many drugs act by interacting with proteins. The "
                "protein we investigate here is called the TARGET.\n\n"
                "In this step you type the name of your protein. The "
                "software looks it up and tells you one very important "
                "thing: has anyone ever made a drug against this protein "
                "before?\n\n"
                "Why that matters:\n\n"
                "If YES - read what kind of drug and evidence exist. "
                "For this small-molecule exercise, a structure with a "
                "known bound molecule is especially useful.\n\n"
                "If NO - the search has not established an existing "
                "drug. That alone does not tell us whether a pocket "
                "exists or why a development programme has not succeeded.\n\n"
                "Example: type EGFR. You will see it has many approved "
                "drugs (gefitinib, erlotinib) and hundreds of structures. "
                "That is a well-studied, drug-friendly target."
            ),
            not_mean=[
                "A strong disease association is not evidence that blocking "
                "the protein helps. Many associations are markers, not "
                "causes, and some protective-looking associations reverse "
                "when you intervene.",
                "No existing drugs does not establish an open opportunity "
                "or the absence of a pocket. Evidence may be incomplete, "
                "and development can fail for many biological, safety or "
                "practical reasons.",
                "Many existing drugs does not mean the field is closed. It "
                "provides useful evidence of tractability, but it does not "
                "predict whether your new molecule will work.",
                "Nothing on this screen tells you whether a new molecule "
                "will work. It tells you whether the question is worth "
                "asking.",
            ],
            learn_more=[
                ("Open Targets - how target evidence is scored",
                 "https://platform-docs.opentargets.org/"),
                ("UniProt entry view",
                 "https://www.uniprot.org/help/uniprotkb"),
            ],
        )

    # -- lookup -----------------------------------------------------------

    def resolve(self, query: str) -> tuple[dict | None, str]:
        """
        Turn whatever the user typed into one UniProt entry.

        Accepts an accession directly; otherwise searches by gene symbol
        among reviewed human entries, which is what a pharmacologist
        almost always means.
        """
        query = query.strip()
        if not query:
            return None, "Nothing was entered."

        looks_like_accession = (
            6 <= len(query) <= 10
            and query[0].isalpha()
            and any(c.isdigit() for c in query)
            and query.isalnum()
            and query.isupper()
        )

        if looks_like_accession:
            blob, res = self.web.get_json(
                f"{UNIPROT}/{query}.json", service="UniProt"
            )
            if blob:
                return parse_uniprot(blob), ""
            if res.status == 404:
                pass        # fall through to a search
            elif not res.ok:
                return None, res.problem

        q = urllib.parse.quote(
            f'(gene:{query} OR protein_name:"{query}") '
            f"AND organism_id:9606 AND reviewed:true"
        )
        blob, res = self.web.get_json(
            f"{UNIPROT}/search?query={q}&format=json&size=5",
            service="UniProt",
        )
        if not res.ok:
            return None, res.problem
        hits = _dig(blob, "results", default=[]) or []
        if not hits:
            return None, (
                f"No reviewed human protein was found for '{query}'. "
                "Try the gene symbol, such as EGFR, or a UniProt accession."
            )
        return parse_uniprot(hits[0]), ""

    def tractability(self, ensembl_gene: str) -> dict[str, Any]:
        """
        Ask Open Targets how druggable this protein is.

        This feeds the mandatory control, so it is the one network call
        in this stage whose failure the user is told about plainly
        rather than silently.
        """
        if not ensembl_gene:
            return {"available": False,
                    "why": "No Ensembl gene id was found for this protein."}

        payload = json.dumps({
            "query": TRACTABILITY_QUERY,
            "variables": {"id": ensembl_gene},
        }).encode()

        result = self.web.post(
            OPENTARGETS, payload, service="Open Targets",
            content_type="application/json",
        )
        if not result.ok:
            return {"available": False, "why": result.problem}

        blob = result.json()
        if not blob:
            return {"available": False,
                    "why": "Open Targets answered with something unreadable."}
        if blob.get("errors"):
            # It answered and refused. Say which, because the two are
            # not the same thing and confusing them wasted a round of
            # work on this project already.
            msg = _dig(blob, "errors", 0, "message", default="") or ""
            return {"available": False,
                    "why": f"Open Targets rejected the request: {msg[:160]}"}

        return parse_tractability(blob)

    def chembl_enrichment(self, accession: str) -> dict[str, Any]:
        """
        Optional extra: which actions drugs take at this target.

        ChEMBL is unreliable -- demoted from the control for that
        reason -- so this is allowed to fail and nothing depends on it.
        One short attempt, no retry, and the user is not interrupted.
        """
        out: dict[str, Any] = {"available": False, "chembl_id": "",
                               "mechanisms": [], "n_mechanism_records": 0}

        blob, res = self.web.get_json(
            f"{CHEMBL}/target.json?target_components__accession={accession}"
            "&limit=5",
            service="ChEMBL", timeout=8, retries=0,
        )
        if not res.ok or not blob:
            return out

        targets = _dig(blob, "targets", default=[]) or []
        single = [t for t in targets if t.get("target_type") == "SINGLE PROTEIN"]
        chembl_id = (single or targets or [{}])[0].get("target_chembl_id", "")
        if not chembl_id:
            return out
        out["chembl_id"] = chembl_id
        out["available"] = True

        mblob, mres = self.web.get_json(
            f"{CHEMBL}/mechanism.json?target_chembl_id={chembl_id}&limit=100",
            service="ChEMBL", timeout=8, retries=0,
        )
        if mres.ok and mblob:
            mechs = _dig(mblob, "mechanisms", default=[]) or []
            out["n_mechanism_records"] = len(mechs)
            seen: set[str] = set()
            for m in mechs:
                a = m.get("action_type", "")
                if a and a not in seen:
                    seen.add(a)
                    out["mechanisms"].append(a)
        return out

    # -- run ---------------------------------------------------------------

    def run(self, query: str = "", **kwargs: Any) -> StageResult:
        entry, problem = self.resolve(query)
        if entry is None:
            return StageResult(success=False, message=problem)

        tract = self.tractability(entry.get("ensembl_gene", ""))
        chembl = self.chembl_enrichment(entry["accession"])
        history = {**tract, "chembl": chembl}

        self.record(
            ProvenanceEntry(
                stage=self.stage_id,
                tool="UniProt REST",
                version="current",
                parameters={"query": query, "accession": entry["accession"]},
                source_release="UniProtKB, retrieved live",
                citation="The UniProt Consortium. Nucleic Acids Res. "
                         "doi:10.1093/nar/gkae1010",
            )
        )
        if tract.get("available"):
            self.record(
                ProvenanceEntry(
                    stage=self.stage_id,
                    tool="Open Targets Platform (GraphQL)",
                    version="v4",
                    parameters={"ensembl_gene": entry.get("ensembl_gene", ""),
                                "used_for": "small-molecule tractability"},
                    citation="Ochoa D et al. Open Targets Platform. Nucleic "
                             "Acids Res. 2023;51:D1353-D1359.",
                )
            )
        if chembl.get("available"):
            self.record(
                ProvenanceEntry(
                    stage=self.stage_id,
                    tool="ChEMBL REST",
                    version="current",
                    parameters={"target_chembl_id": chembl["chembl_id"]},
                    citation="Zdrazil B et al. ChEMBL. Nucleic Acids Res. "
                             "2024;52:D1180-D1192.",
                )
            )
        for e in self.provenance:
            self.project.add_provenance(e)

        return StageResult(
            success=True,
            summary=(
                f"{entry['protein'] or entry['accession']} "
                f"({', '.join(entry['genes']) or 'no gene symbol'}), "
                f"{entry['length']} residues, "
                f"{len(entry['pdb_ids'])} experimental structure(s)"
            ),
            values={"target": entry, "history": history},
            control=self._tractability_control(entry, history),
        )

    # -- the control -------------------------------------------------------

    @staticmethod
    def _tractability_control(entry: dict, history: dict) -> ControlOutcome:
        """
        Specification Stage 1, rebuilt on the Open Targets evidence
        ladder rather than a count of molecules.

        The ladder is the point. "Approved Drug" and "Druggable Family"
        are both evidence, and they are not remotely the same evidence:
        one means a marketed drug binds this protein, the other means
        somebody drugged a relative of it. A count conflates them; an
        ordered ladder does not.

        An empty ladder reports missing evidence in the queried source;
        it cannot establish that no binding pocket or prior work exists.
        """
        n_pdb = len(entry.get("pdb_ids", []))
        chembl = history.get("chembl", {}) or {}
        base = {
            "pdb_structures": n_pdb,
            "evidence_ladder": history.get("sm_buckets", []),
        }
        if chembl.get("available") and chembl.get("mechanisms"):
            base["drug_actions_recorded"] = chembl["mechanisms"][:6]

        if not history.get("available"):
            return ControlOutcome(
                name="Has anyone drugged this before?",
                verdict=Verdict.NOT_APPLICABLE,
                plain_verdict=(
                    "Druggability could not be checked, so the most "
                    "informative thing you can know about a target is "
                    "missing. "
                    + (history.get("why", "") or "")
                    + " Find this out before committing time to it."
                ),
                what_to_do=(
                    "Check your internet connection and run this step "
                    "again. If it still fails, you can continue - but you "
                    "are now working without the single most useful piece "
                    "of information about your target."
                ),
                detail=base,
            )

        clinical = history.get("clinical", [])
        structural = history.get("structural", [])
        buckets = history.get("sm_buckets", [])

        if clinical:
            top = clinical[0]
            return ControlOutcome(
                name="Has anyone drugged this before?",
                verdict=Verdict.PASS,
                plain_verdict=(
                    f"Yes - this protein has reached the clinic "
                    f"({top.lower()}), and {n_pdb} structure(s) exist. A "
                    "target with clinical evidence is a useful teaching "
                    "starting point. Check that the recorded mechanism is "
                    "relevant to small-molecule docking."
                ),
                what_to_do=(
                    "Good choice of target. Press Next and get its "
                    "3D structure."
                ),
                detail=base,
            )

        if structural:
            return ControlOutcome(
                name="Has anyone drugged this before?",
                verdict=Verdict.PASS,
                plain_verdict=(
                    "No drug has reached the clinic against this protein, but "
                    f"there is real structural evidence of a bindable pocket "
                    f"({', '.join(structural).lower()}). That is a reasonable "
                    "place to start - unproven, not untouched."
                ),
                what_to_do=(
                    "Reasonable to continue. Just keep in mind that nobody "
                    "has turned this pocket into a medicine yet, so you are "
                    "on less certain ground than with a proven target."
                ),
                detail=base,
            )

        if buckets:
            return ControlOutcome(
                name="Has anyone drugged this before?",
                verdict=Verdict.WARN,
                plain_verdict=(
                    "The evidence returned by this search is "
                    f"indirect ({', '.join(buckets).lower()}) - it comes from "
                    "relatives of this protein rather than direct evidence "
                    "for this protein. This does not establish whether it "
                    "has a suitable pocket. Review the literature and "
                    "structures before interpreting docking results."
                ),
                what_to_do=(
                    "Think once before continuing. If this target is your "
                    "own choice, consider a better-studied one instead - "
                    "search for a related protein in the same family. If it "
                    "is fixed (your thesis topic, your guide's project), "
                    "carry on, but you must run the two automatic checks at "
                    "the docking step and report their result honestly in "
                    "your write-up."
                ),
                detail=base,
                blocks_progress=True,
            )

        return ControlOutcome(
            name="Has anyone drugged this before?",
            verdict=Verdict.WARN,
            plain_verdict=(
                "This search returned no small-molecule tractability "
                "evidence for the protein or its family. Missing database "
                "evidence does not show that a pocket is absent or that "
                "previous attempts failed. A docking score alone cannot "
                "fill that evidence gap."
            ),
            what_to_do=(
                "For a first exercise, consider a better-studied target. "
                "If you continue, review the literature and structures, "
                "record the evidence gap, and treat poses as hypotheses "
                "that need experimental testing."
            ),
            detail=base,
            blocks_progress=True,
        )

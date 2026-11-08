"""
The input each stage needs, and how to run it.

Kept apart from the page layout for one reason: these are the only
places where a user decision becomes a parameter, so this is where the
defaults live, and the defaults are most of what protects a novice.

Every default here is the cautious choice, and where a cautious default
costs something the field says so rather than hiding it.

Every field also answers three questions for a reader who has never
done this before: what is this thing, what should I put, and what
happens if I get it wrong. A box labelled "Chains to keep: A" with no
explanation is where a student stops.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

from nicegui import ui


HELP = "font-size:0.85rem; color:#5f6368; line-height:1.5"
WARN = "font-size:0.85rem; color:#b26a00; line-height:1.5"


def _help(text: str, style: str = HELP) -> None:
    """One line of guidance under a field."""
    ui.label(text).style(style).classes("q-mb-sm")


def _more(title: str, text: str) -> None:
    """
    Optional depth, folded away.

    The short line above has to be enough to fill the box in. This is
    for the student who wants to know why -- which should never be the
    price of getting on with the work.
    """
    with ui.expansion(title).classes("w-full q-mb-sm").props("dense"):
        ui.label(text).style("font-size:0.85rem; line-height:1.6")


def chains_in_pdb(path: str | Path) -> list[tuple[str, int]]:
    """
    Which chains are in this file, and how big is each one.

    Needed so the chain question can be asked concretely -- "this file
    has chain A with 327 amino acids and chain B with 4" -- instead of
    asking a novice to type a letter they have no way of choosing.
    """
    try:
        text = Path(path).read_text(errors="ignore")
    except Exception:  # noqa: BLE001
        return []
    seen: dict[str, set[str]] = {}
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        chain = line[21:22].strip() or "?"
        residue = line[22:27]            # number + insertion code
        seen.setdefault(chain, set()).add(residue)
    return sorted(((c, len(r)) for c, r in seen.items()),
                  key=lambda t: -t[1])


# ----------------------------------------------------------------------

class StageForm:
    """Inputs for one stage, and the values they produce."""

    def __init__(self, stage_id: str, project: Any) -> None:
        self.stage_id = stage_id
        self.project = project
        self.fields: dict[str, Any] = {}

    # -- building ---------------------------------------------------------

    def build(self) -> None:
        builder = getattr(self, f"_build_{self.stage_id[3:]}", None)
        if builder is None:
            ui.label("This step takes no input.").style("color:#5f6368")
            return
        builder()

    def values(self) -> dict[str, Any]:
        out: dict[str, Any] = {}
        for key, field in self.fields.items():
            value = field.value if hasattr(field, "value") else field
            out[key] = value
        return self._clean(out)

    def _clean(self, raw: dict[str, Any]) -> dict[str, Any]:
        cleaner = getattr(self, f"_clean_{self.stage_id[3:]}", None)
        return cleaner(raw) if cleaner else raw

    # -- stage 1: target --------------------------------------------------

    def _build_target(self) -> None:
        saved = self.project.get_step("01_target")
        target = json.loads(saved["values_json"] or "{}").get("target", {}) if saved else {}
        self.fields["query"] = ui.input(
            "Name of your protein",
            value=next(iter(target.get("genes", [])), ""),
            placeholder="EGFR",
        ).props("outlined").classes("w-full")
        _help(
            "Type the short gene name - EGFR, ACE, PTGS2, HMGCR. That is "
            "what works best. A full name (\"epidermal growth factor "
            "receptor\") also works. Human proteins are searched first."
        )
        _more(
            "I do not know my protein's gene name",
            "Think of the drug first, then its target. Atorvastatin acts "
            "on HMG-CoA reductase, whose gene name is HMGCR. Aspirin acts "
            "on cyclo-oxygenase, gene name PTGS1 and PTGS2. Enalapril "
            "acts on angiotensin converting enzyme, gene name ACE. If you "
            "are not sure of the short name, type the full name and the "
            "software will show you the matches to choose from."
        )

    # -- stage 2: structure -----------------------------------------------

    def _build_structure(self) -> None:
        self.fields["prefer"] = ui.select(
            {"auto": "Let the software choose (recommended)",
             "experimental": "Only real laboratory structures",
             "predicted": "Only the AlphaFold computer model"},
            value="auto", label="Which 3D structure to use",
        ).props("outlined").classes("w-full")
        _help(
            "Start with the first option. The software prefers a suitable "
            "ligand-bound experimental structure when available. Review "
            "the selected protein, binding site and structure quality."
        )
        self.fields["pdb_id"] = ui.input(
            "Or use a particular structure (leave empty to let the "
            "software choose)",
            placeholder="9OKG",
        ).props("outlined").classes("w-full")
        _help(
            "Only if you already know which one you want - from a paper, "
            "or because your guide asked for it. Otherwise leave it "
            "empty."
        )
        _more(
            "Selecting a suitable structure",
            "Run this step once with the box empty. The software picks "
            "one and, under 'The details', lists the others it "
            "considered with the numbers for each. Then judge them in "
            "this order:\n\n"
            "1. RELEVANT BOUND LIGAND. A measured complex shows a "
            "ligand-bound conformation and provides a reference pose for "
            "re-docking. Check that it is the site and mechanism you "
            "intend to study. An unbound structure can adopt a different "
            "shape, but is not automatically unsuitable.\n\n"
            "2. LOCAL SITE QUALITY. Lower reported resolution "
            "can support finer detail, but a global number does not "
            "guarantee local atom positions. X-ray and cryo-EM use "
            "different quality measures; inspect the local density, "
            "model validation and binding-site residues. The app's "
            "resolution thresholds are screening rules.\n\n"
            "3. MISSING REGIONS NEAR THE SITE. Location matters as "
            "well as the overall fraction missing. A missing distant "
            "tail may be less relevant to a local docking exercise, "
            "but do not assume all distant regions are unimportant. "
            "The preparation step estimates gap proximity to a bound "
            "ligand where possible; inspect unresolved or rebuilt "
            "regions near your site.\n\n"
            "4. PROTEIN IDENTITY AND COVERAGE. The software "
            "already refuses structures that contain only a short "
            "fragment.\n\n"
            "Choose using all of this evidence rather than resolution "
            "alone. A relevant ligand-bound structure can be especially "
            "helpful, provided the local model is suitable."
        )
        _more(
            "Structures with a bound reference molecule",
            "A relevant experimental complex shows a ligand-bound "
            "conformation and provides reference coordinates. Include "
            "that ligand among your docking molecules to test whether "
            "the setup recovers its pose. Passing supports pose "
            "reproduction for that ligand; it does not validate affinity "
            "predictions or rankings of new molecules. Without a usable "
            "reference, record that this test could not be performed.\n\n"
            "Experimental structures can come from X-ray crystallography, "
            "cryo-EM and other methods. The AlphaFold model retrieved by "
            "this app is predicted and does not supply a measured "
            "ligand pose. Inspect its confidence and site suitability."
        )

    def _clean_structure(self, raw: dict[str, Any]) -> dict[str, Any]:
        return {
            "prefer": raw.get("prefer", "auto"),
            "pdb_id": str(raw.get("pdb_id") or "").strip(),
        }

    # -- stage 3: prepare protein -----------------------------------------

    def _build_receptor(self) -> None:
        structure = self.project.get_step("02_structure")
        chosen: dict[str, Any] = {}
        if structure:
            chosen = json.loads(
                structure["values_json"] or "{}").get("chosen", {}) or {}
        pdb_id = chosen.get("id", "")
        ui.label(f"Working on structure: {pdb_id or 'none chosen yet'}").style(
            "font-weight:500").classes("q-mb-sm")

        # -- which chain -------------------------------------------------
        chains = chains_in_pdb(chosen.get("local_file", ""))
        protein_chains = [(c, n) for c, n in chains if n >= 30]

        if protein_chains:
            options: dict[str, str] = {}
            for i, (c, n) in enumerate(protein_chains):
                tag = " - recommended" if i == 0 else ""
                options[c] = f"Chain {c}  ({n} amino acids){tag}"
            if len(protein_chains) > 1:
                options["*"] = "Keep all chains (only if they work together)"
            default = protein_chains[0][0]
            self.fields["chains"] = ui.select(
                options, value=default, label="Which copy of the protein to use",
            ).props("outlined").classes("w-full")

            if len(protein_chains) == 1:
                _help(
                    f"This file has only one copy of the protein, labelled "
                    f"chain {default}. Nothing to decide - carry on."
                )
            else:
                names = ", ".join(c for c, _ in protein_chains)
                _help(
                    f"This file contains {len(protein_chains)} copies of the "
                    f"protein, labelled {names}. Chain {default} is "
                    "the longest and is selected initially. Check the "
                    "biological assembly and whether other chains "
                    "contribute to your binding site."
                )
        else:
            self.fields["chains"] = ui.input(
                "Which copy of the protein to use", value="A",
            ).props("outlined").classes("w-full")
            _help(
                "Type the chain letter, usually A. (The chain list could "
                "not be read from the file, so this box is free text.)"
            )

        _more(
            "Protein chains and chain labels",
            "A chain is a labelled polymer strand in the structure "
            "file. Different chains may be copies of the same protein "
            "or different proteins. A, B and other identifiers are "
            "labels, not quality rankings.\n\n"
            "Choose the assembly that represents the site you want to "
            "study. One chain may be sufficient for a site within that "
            "chain, while an interface site needs its contributing "
            "chains. HIV protease is a familiar two-chain example.\n\n"
            "Check the PDB biological assembly and literature with your "
            "instructor if the choice is unclear. The longest chain is "
            "only a starting default; an enclosure score cannot by "
            "itself determine the correct biological assembly."
        )

        # -- waters ------------------------------------------------------
        self.fields["remove_waters"] = ui.checkbox(
            "Remove water molecules", value=True)
        _help(
            "Water removal is this teaching workflow's default. Review "
            "whether conserved waters are important to your site."
        )
        _more(
            "Water removal and retention",
            "Experimental structures may contain ordered waters. This "
            "workflow normally removes them as a simplification; how "
            "water should be modelled depends on the docking method "
            "and scientific question.\n\n"
            "When it matters: occasionally a single water sits between "
            "the drug and the protein and passes a hydrogen bond between "
            "them, like a bridge. That water is really part of the "
            "pocket. This is well known in HIV protease and some kinases. "
            "If you are working on a protein where the literature "
            "mentions a conserved or structural water, review it with "
            "your instructor before choosing what to retain. Keeping "
            "all waters is not automatically the better model."
        )

        # -- gaps --------------------------------------------------------
        self.fields["max_gap_length"] = ui.number(
            "Repair breaks in the chain up to this length", value=3,
            min=0, max=20, format="%.0f",
        ).props("outlined").classes("w-full")
        _help(
            "Leave it at 3. This means: if up to 3 amino acids are "
            "missing in a row, rebuild them. If more than 3 are missing, "
            "leave the hole and tell me about it."
        )
        _more(
            "Chain breaks and the repair limit",
            "A protein is one long chain of amino acids. In an X-ray "
            "structure, some parts of that chain wobble so much that the "
            "X-rays cannot see them at all - usually floppy loops on the "
            "surface. Those amino acids are simply absent from the file, "
            "leaving a hole in the middle of the chain.\n\n"
            "This workflow attempts to rebuild gaps up to the selected "
            "length. Rebuilt atoms are modelled, not observed, even for "
            "short gaps. A longer gap generally adds more uncertainty; "
            "gaps above the limit are left unresolved and reported.\n\n"
            "Missing or modelled residues near the site can alter its "
            "shape and predicted contacts. Inspect their location before "
            "interpreting a docking pose.\n"
            "\n"
            "Set it to 0 if you want nothing rebuilt at all. Setting it "
            "to 15 or 20 is not recommended."
        )

    def _clean_receptor(self, raw: dict[str, Any]) -> dict[str, Any]:
        structure = self.project.get_step("02_structure")
        values = json.loads(structure["values_json"] or "{}") if structure else {}
        chosen = values.get("chosen", {})
        path = chosen.get("local_file", "")
        value = str(raw.get("chains", "")).strip()
        if value == "*":
            chains: list[str] = []
        else:
            chains = [c.strip() for c in value.split(",") if c.strip()]
        return {
            "pdb_path": path,
            "keep_chains": chains or None,
            "remove_waters": bool(raw.get("remove_waters", True)),
            "max_gap_length": int(raw.get("max_gap_length") or 3),
        }

    # -- stage 4: binding site --------------------------------------------

    def _build_site(self) -> None:
        receptor = self.project.get_step("03_receptor")
        ligand = ""
        chains: list[str] = []
        if receptor:
            values = json.loads(receptor["values_json"] or "{}")
            ligand = values.get("ligand", "")
            chains = sorted((values.get("chains") or {}).keys())

        options: dict[str, str] = {}
        if ligand:
            options["ligand"] = f"Where {ligand} sits  - recommended"
        options["detect"] = "Let the software hunt for pockets"
        options["residues"] = "I know which residues line the pocket"
        options["manual"] = "I will type the box myself"

        self.fields["mode"] = mode = ui.select(
            options, value="ligand" if ligand else "detect",
            label="Binding-site selection method",
        ).props("outlined").classes("w-full")

        if ligand:
            _help(
                f"The molecule {ligand} provides reference coordinates "
                "for this structure. The software can place the box "
                "there; check that this is the site you intend to study."
            )
        else:
            _help(
                "There is no drug bound in this structure, so there is "
                "nothing to copy the position from. Start with the "
                "automatic hunt. If it finds nothing - which happens "
                "often, and is an honest answer rather than a failure - "
                "use one of the other two options below.",
                WARN,
            )

        # -- the fallback that stops the dead end ------------------------
        #
        # Automatic detection fails on shallow sites, on sites between
        # two chains, and on anything that only opens when a drug
        # arrives. Telling the user to "go back and set the box again"
        # when there is nothing to go back to is not advice. These two
        # panels are the way through, and they show only when chosen.
        with ui.column().classes("w-full q-gutter-none").bind_visibility_from(
                mode, "value", lambda v: v == "residues"):
            self.fields["site_chain"] = ui.input(
                "Chain", value=(chains[0] if chains else "A"),
            ).props("outlined dense").classes("w-full")
            self.fields["site_residues"] = ui.input(
                "Residue numbers, separated by commas",
                placeholder="745, 790, 793, 855",
            ).props("outlined").classes("w-full")
            _help(
                "Type the residue numbers that line the pocket. The "
                "software will draw the box around them. Four or five "
                "is plenty; they should be residues that face each "
                "other across the pocket, not a stretch along the "
                "chain."
            )
            _more(
                "Finding binding-site residue numbers",
                "From the literature on your protein. Almost any paper "
                "on its active site names them, usually in a figure "
                "caption: \"the inhibitor hydrogen bonds to Met793 and "
                "Thr790\". UniProt also lists them - open your "
                "protein's entry and look under 'Binding site' and "
                "'Active site'.\n\n"
                "One warning that catches everybody. Papers use the "
                "numbering of the full protein. A crystal structure "
                "often contains only one domain and sometimes renumbers "
                "it. If the software says it cannot find your residue "
                "numbers, that is what has happened - open the "
                "structure's page on rcsb.org and check which numbers "
                "it actually uses.\n\n"
                "After running, inspect the enclosure figure and 3D "
                "view. Enclosure is a geometric diagnostic, not proof "
                "of a binding site. Shallow or interface sites may "
                "have lower enclosure than a buried cavity."
            )

        with ui.column().classes("w-full q-gutter-none").bind_visibility_from(
                mode, "value", lambda v: v == "manual"):
            _help(
                "For coordinates you already have - from PyMOL, from a "
                "paper's methods section, or from a colleague's config "
                "file."
            )
            with ui.row().classes("w-full no-wrap q-gutter-sm"):
                for axis in ("x", "y", "z"):
                    self.fields[f"center_{axis}"] = ui.number(
                        f"Centre {axis}", value=0.0, step=0.5,
                        format="%.2f",
                    ).props("outlined dense").style("flex:1")
            with ui.row().classes("w-full no-wrap q-gutter-sm"):
                for axis in ("x", "y", "z"):
                    self.fields[f"size_{axis}"] = ui.number(
                        f"Size {axis}", value=22.0, min=6, max=40, step=1,
                        format="%.0f",
                    ).props("outlined dense").style("flex:1")
            _help(
                "22 x 22 x 22 angstroms is a sensible default size. In "
                "PyMOL you can get the centre with: "
                "print(cmd.centerofmass('resi 745+790+793'))"
            )

        _more(
            "Docking search box",
            "A protein might be 300 amino acids. A drug is 30 atoms. If "
            "you asked the computer to try the drug everywhere on the "
            "protein, the larger search space could require much more "
            "effort and include unintended sites.\n\n"
            "So you draw an invisible cube around the pocket you care "
            "about, and the computer only searches inside that cube. "
            "The cube is called the grid box or the search box.\n\n"
            "Choosing the box is a scientific decision. A box in the "
            "wrong place can still "
            "gives you a confident-looking number."
        )

        self.fields["padding"] = ui.number(
            "Extra space around the pocket (angstroms)", value=4.0,
            min=0, max=15, step=0.5, format="%.1f",
        ).props("outlined").classes("w-full").bind_visibility_from(
            mode, "value", lambda v: v != "manual")
        _help(
            "Leave it at 4. This adds a 4 angstrom margin around the "
            "pocket so your molecule has room to turn, without making "
            "the box so big that the search gets thin."
        )
        _more(
            "Search-box size and sampling",
            "A larger search space can need more sampling. With "
            "unchanged settings it may miss useful poses or explore an "
            "unintended region. The effect is not a simple fixed amount "
            "per unit volume.\n\n"
            "Cover the intended site and leave enough room for the "
            "ligand to move. Inspect the box in 3D and use the reference "
            "pose check when available. An angstrom is one "
            "ten-billionth of a metre; the padding here is a distance "
            "around the selected site, not a quality guarantee."
        )

    def _clean_site(self, raw: dict[str, Any]) -> dict[str, Any]:
        receptor = self.project.get_step("03_receptor")
        values = json.loads(receptor["values_json"] or "{}") if receptor else {}
        mode = raw.get("mode", "detect")
        out: dict[str, Any] = {
            "pdb_path": values.get("prepared_file", ""),
            "padding": float(raw.get("padding") or 4.0),
        }
        if mode == "ligand":
            out["ligand_label"] = values.get("ligand", "")
        elif mode == "residues":
            out["chain"] = str(raw.get("site_chain") or "A").strip()
            out["residues"] = [
                r.strip()
                for r in str(raw.get("site_residues") or "").replace(
                    ";", ",").split(",")
                if r.strip()
            ]
        elif mode == "manual":
            out["box"] = {
                "center_x": float(raw.get("center_x") or 0.0),
                "center_y": float(raw.get("center_y") or 0.0),
                "center_z": float(raw.get("center_z") or 0.0),
                "size_x": float(raw.get("size_x") or 22.0),
                "size_y": float(raw.get("size_y") or 22.0),
                "size_z": float(raw.get("size_z") or 22.0),
            }
        return out

    # -- stage 5: prepare molecules ---------------------------------------

    def _build_ligands(self) -> None:
        saved = self.project.get_step("05_ligands")
        prepared = json.loads(saved["values_json"] or "{}").get("prepared", []) if saved else []
        self.fields["molecules"] = ui.textarea(
            "Your molecules - one per line, SMILES code then the name",
            value="\n".join(f"{m.get('input_smiles') or m['smiles']} {m['name']}" for m in prepared),
            placeholder=(
                "CC(=O)Oc1ccccc1C(=O)O aspirin\n"
                "CC(C)Cc1ccc(cc1)C(C)C(O)=O ibuprofen"
            ),
        ).props("outlined rows=6").classes("w-full")
        _help(
            "One molecule per line. First the SMILES code, then a space, "
            "then any name you like. The name is only for your own "
            "reading."
        )
        _more(
            "Chemical structures as SMILES",
            "SMILES is a way of writing a chemical structure as a line of "
            "text, so a computer can read it. Aspirin is "
            "CC(=O)Oc1ccccc1C(=O)O. The letters are atoms, the brackets "
            "are branches, the small letters mean an aromatic ring.\n\n"
            "You do not need to learn to write it. Get it like this:\n\n"
            "1. Go to pubchem.ncbi.nlm.nih.gov\n"
            "2. Search your drug by name\n"
            "3. Open the page and look for \"Canonical SMILES\"\n"
            "4. Copy that line and paste it here\n\n"
            "Careful while copying: Cl means chlorine, but CL is wrong. "
            "Do not paste from Word, which inserts invisible characters. "
            "If the software says it cannot read your molecule, copying "
            "it again from the source can help identify copying errors.\n\n"
            "To compare several drugs, put each on its own line."
        )

    def _clean_ligands(self, raw: dict[str, Any]) -> dict[str, Any]:
        molecules = []
        for line in str(raw.get("molecules", "")).splitlines():
            line = line.strip()
            if not line:
                continue
            parts = line.split(None, 1)
            molecules.append({
                "smiles": parts[0],
                "name": parts[1].strip() if len(parts) > 1 else "",
            })
        return {"molecules": molecules}

    # -- stage 6: dock ----------------------------------------------------

    def _build_dock(self) -> None:
        self.fields["exhaustiveness"] = ui.number(
            "Docking search effort", value=8,
            min=1, max=32, format="%.0f",
        ).props("outlined").classes("w-full")
        _help(
            "Start at 8 for the workshop. Higher effort takes longer and "
            "may help the search, but does not guarantee a correct result. "
            "Use the re-docking check to assess this setup."
        )
        _more(
            "Search effort and calculation time",
            "The computer does not try every possible position - there "
            "are far too many. It makes a number of independent attempts, "
            "each starting from a random position, and keeps the best "
            "result.\n\n"
            "Exhaustiveness controls search effort; it is not a direct "
            "count of poses. Higher values take longer and can improve "
            "sampling, but do not guarantee convergence or correctness.\n\n"
            "This app records a fixed seed for repeatability. To assess "
            "search variability, compare runs with different seeds in a "
            "separate validation study; repeating the same fixed-seed "
            "calculation is not an independent check."
        )

        self.fields["n_decoys"] = ui.number(
            "Number of comparison molecules", value=3,
            min=1, max=8, format="%.0f",
        ).props("outlined").classes("w-full")
        _help(
            "Start at 3. The app attempts a small comparison set matched "
            "for size and calculated lipophilicity. These are not "
            "confirmed non-binders, and this exercise does not validate "
            "binding or potency."
        )
        _more(
            "Comparison-molecule checks",
            "A score is easier to interpret alongside a baseline than "
            "on its own. Size, lipophilicity and other molecular "
            "features can influence scores, as can the search and "
            "preparation choices. No score is guaranteed for a particular "
            "kind of molecule or pocket.\n\n"
            "The app picks a few available molecules with similar size "
            "and calculated lipophilicity to the leading compound. "
            "Their activity at this target has not been established. "
            "Failure to separate from this small baseline weakens a "
            "score-based claim; it does not show why the scores are "
            "similar. Passing also does not establish binding.\n\n"
            "Report the comparison set, matching limitations and "
            "outcomes. A larger, independently validated evaluation is "
            "needed to assess screening performance."
        )

    # -- stage 9: molecular dynamics --------------------------------------

    def _build_md(self) -> None:
        dock = self.project.get_step("06_dock")
        results = []
        if dock:
            results = [r for r in json.loads(
                dock["values_json"] or "{}").get("results", [])
                if "score" in r]

        if results:
            best = min(results, key=lambda r: r["score"])
            options = {
                r.get("name") or "molecule":
                f"{r.get('name') or 'molecule'}  ({r['score']} kcal/mol)"
                + ("  - best" if r is best else "")
                for r in results
            }
            self.fields["molecule"] = ui.select(
                options, value=best.get("name") or "molecule",
                label="Which molecule to simulate",
            ).props("outlined").classes("w-full")
            _help(
                "One molecule at a time. Simulating is expensive, so "
                "review the docking checks and pose plausibility before "
                "choosing a starting model."
            )
        else:
            ui.label("Nothing has been docked yet - do step 6 first.").style(
                "font-weight:500; color:#b26a00").classes("q-mb-sm")

        self.fields["nanoseconds"] = ui.number(
            "How long to simulate (nanoseconds)", value=100,
            min=1, max=1000, format="%.0f",
        ).props("outlined").classes("w-full")
        _help(
            "100 ns is the package default, not a guarantee of adequate "
            "sampling. Choose duration, equilibration and replicas for "
            "your question. A shorter pilot can help check the setup."
        )
        _more(
            "Interpreting a simulation",
            "A nanosecond is a billionth of a second. Different "
            "molecular motions occur on different timescales; a fixed "
            "duration cannot guarantee that the motion of interest "
            "will be sampled.\n\n"
            "Inspect displacement, contacts and trajectories after "
            "equilibration, and compare independent replicas. Small "
            "ligand RMSD in one run does not establish binding, "
            "activity or convergence. Large RMSD needs inspection of "
            "the model and trajectory before attributing its cause.\n\n"
            "Report what was observed under the recorded conditions "
            "and over the simulated time. An event not observed in "
            "100 ns has not been shown to be impossible."
        )

        self.fields["pose_rank"] = ui.number(
            "Which pose from the docking", value=1, min=1, max=9,
            format="%.0f",
        ).props("outlined").classes("w-full")
        _help(
            "Leave it at 1 - the best-scoring pose. Change it only if "
            "you looked at the poses and pose 2 made more chemical "
            "sense than pose 1."
        )

    def _clean_md(self, raw: dict[str, Any]) -> dict[str, Any]:
        return {
            "molecule": str(raw.get("molecule") or ""),
            "nanoseconds": float(raw.get("nanoseconds") or 100),
            "pose_rank": int(raw.get("pose_rank") or 1),
        }

    # -- stage 8: pharmacokinetics ----------------------------------------

    def _build_admet(self) -> None:
        ligands = self.project.get_step("05_ligands")
        n = 0
        if ligands:
            n = len(json.loads(
                ligands["values_json"] or "{}").get("prepared", []))
        ui.label(
            f"{n} molecule(s) ready to profile."
            if n else "No molecules prepared yet - do step 5 first."
        ).style("font-weight:500").classes("q-mb-sm")
        _help(
            "Nothing to fill in. Press Run and the software will work out "
            "the drug-like properties of your molecules from their "
            "structure."
        )
        _more(
            "Interpreting calculated descriptors",
            "Start with the structure and method used for each value. "
            "Molecular weight and atom counts follow the entered "
            "structure. Hydrogen-bond counts and polar surface area "
            "depend on descriptor rules. Calculated logP estimates "
            "lipophilicity with an empirical method; it is not a "
            "measured partition coefficient. Ionisation and tautomer "
            "choices can change the inputs and results.\n\n"
            "Optional absorption, brain-penetration and enzyme "
            "predictions come from separate models when available. "
            "Read each model's validation, applicability information "
            "and uncertainty before using its output. These predictions "
            "do not establish clinical exposure or safety.\n\n"
            "The reference-range check only flags unusual descriptor "
            "values under this exercise's rules. Passing does not "
            "establish that any prediction model applies to the "
            "molecule or that the molecule is a good drug."
        )

    def _clean_admet(self, raw: dict[str, Any]) -> dict[str, Any]:
        return {}

    def _clean_dock(self, raw: dict[str, Any]) -> dict[str, Any]:
        return {
            "exhaustiveness": int(raw.get("exhaustiveness") or 8),
            "n_decoys": int(raw.get("n_decoys") or 3),
        }

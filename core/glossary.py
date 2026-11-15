"""
The words a student meets in this application, in plain English.

Every entry is written for somebody who has never docked anything and
has no chemistry beyond an MBBS or BPharm course. Two rules were
followed throughout: say what the thing IS before saying what it is
for, and never define a term using another term from this same list
without explaining that one too.

Where a word is routinely over-read -- score, RMSD, exhaustiveness --
the entry says what it does NOT tell you. That is the part students
need and textbooks skip.
"""

from __future__ import annotations

GLOSSARY: list[tuple[str, str, str]] = [
    # (term, short plain meaning, the bit that is usually got wrong)
    ("Affinity",
     "How tightly a drug sticks to its target, usually written as Kd or "
     "IC50. A smaller number means stronger binding.",
     "Docking does not measure this. It estimates it, with an error "
     "larger than the range it is trying to measure."),

    ("Apo and holo",
     "Apo means the structure was solved with nothing bound in the "
     "pocket. Holo means a molecule was bound when the structure was "
     "solved.",
     "Prefer holo. The pocket is already in the shape that accepts a "
     "drug, and the bound molecule gives you a free way to test your "
     "setup."),

    ("Binding site, or pocket",
     "The dent in the protein's surface where the drug sits. A good one "
     "is deep enough to surround a drug-sized molecule on most sides.",
     "Only about one in ten human proteins has such a pocket. The rest "
     "are called undruggable, and docking against them still produces "
     "numbers."),

    ("Cofactor",
     "A non-protein helper the enzyme needs in order to work: haem, NAD, "
     "FAD, ATP, or a metal such as zinc.",
     "Keep it. An enzyme without its cofactor has a different pocket, "
     "and deleting a catalytic zinc is the single most damaging mistake "
     "in structure preparation."),

    ("Decoy",
     "A molecule you are confident does NOT bind your target, docked on "
     "purpose to see what a meaningless score looks like.",
     "It must be about the same size as your real compound. A big "
     "inhibitor beating small aspirin proves nothing at all."),

    ("Docking score",
     "The program's estimate of binding free energy, in kcal/mol. More "
     "negative means it thinks the binding is stronger.",
     "Treat differences under about 3 kcal/mol as no difference. The "
     "error bar is roughly 2.5."),

    ("Exhaustiveness",
     "How many separate attempts the search makes, each starting from a "
     "random position. The usual setting is 8.",
     "If two runs of the same job give scores more than about 0.2 "
     "kcal/mol apart, the search has not settled. Raise this and make "
     "the box smaller."),

    ("Force field",
     "The set of rules a simulation uses for how stiff each bond is and "
     "how strongly each pair of atoms attracts or repels.",
     "It is fitted to measurements, not derived from first principles, "
     "so it is accurate inside the chemistry it was fitted to."),

    ("Heavy atom",
     "Any atom that is not hydrogen. The usual way of describing how big "
     "a molecule is in this field.",
     "Docking works best between about 20 and 40 heavy atoms. Below 20, "
     "all the scores bunch together and ranking becomes hopeless."),

    ("Hydrogen bond",
     "A weak attraction between a hydrogen already attached to nitrogen "
     "or oxygen, and another nitrogen or oxygen nearby.",
     "It needs the right angle, not just the right distance. A contact "
     "list based on distance alone can only call it possible."),

    ("Induced fit",
     "The pocket changing shape to accommodate whatever binds in it.",
     "Ordinary docking cannot model this: the protein is held rigid, so "
     "you only find molecules that suit the one shape you downloaded."),

    ("kcal/mol",
     "Kilocalories per mole, the energy unit used for binding. Negative "
     "means favourable.",
     "At body temperature, a ten-fold change in potency is worth only "
     "1.42 kcal/mol. The whole range from an excellent drug to an "
     "inactive one is about 7."),

    ("Ligand",
     "The small molecule that binds. In this application, your drug.",
     "A crystal structure's 'ligand' may also be a buffer molecule or a "
     "sugar that happened to be in the solution."),

    ("Ligand efficiency",
     "The binding energy divided by the number of heavy atoms. Above "
     "about 0.3 is drug-like.",
     "Use it whenever you compare molecules of different sizes, because "
     "raw score rewards bulk."),

    ("Molecular dynamics (MD)",
     "A simulation that moves every atom forward in tiny steps of time, "
     "so you can watch whether a pose holds together.",
     "A simulation that runs without crashing is not evidence that the "
     "drug binds. It is evidence that the model did not fall apart."),

    ("PDB",
     "The Protein Data Bank at rcsb.org, the public archive of "
     "experimentally determined protein structures.",
     "Each entry is one photograph of one molecule under one set of "
     "conditions, not a portrait of the protein in general."),

    ("PDBQT",
     "The file format AutoDock and Vina use. An ordinary coordinate file "
     "with two extra columns: partial charge and atom type.",
     "Most 3D viewers read it poorly. Use the complex.pdb in your export "
     "to look at a result."),

    ("Pose",
     "One specific position, orientation and internal shape of the drug "
     "inside the pocket.",
     "This is the trustworthy output of docking. Usually right to within "
     "2 angstroms when the setup has been checked."),

    ("Re-docking",
     "Taking the molecule that was in the crystal structure, removing "
     "it, docking it back, and measuring how far off the prediction is.",
     "The single most useful check you can run. Without it you have a "
     "number, not a result."),

    ("Resolution",
     "How much detail was visible in the experimental map the structure "
     "was built from, in angstroms. A smaller number is better.",
     "Under 2.5 is comfortable. Beyond 3.0 the side chains lining your "
     "pocket were modelled from the sequence rather than seen."),

    ("RMSD",
     "Root mean square deviation: the average distance between matching "
     "atoms of two structures, in angstroms.",
     "Under 2.0 is the usual mark for a correct pose. It must be "
     "measured in place, not after superimposing the two molecules."),

    ("SMILES",
     "A molecule written as one line of text. Aspirin is "
     "CC(=O)Oc1ccccc1C(=O)O.",
     "A wrong SMILES usually reads perfectly well, because it is a valid "
     "description of a different molecule. Always check the formula."),

    ("Scoring function",
     "The formula that turns a pose into a number. Vina's has six terms "
     "fitted to about 1,300 crystal structures.",
     "It has no term for charge attracting charge, and none for the cost "
     "of stripping water off a polar group."),

    ("Search box",
     "The cube inside which the program is allowed to look. Nothing "
     "outside it is ever considered.",
     "Doubling each side multiplies the volume by eight while the amount "
     "of searching stays the same, so a big box searches worse."),
]


def entries() -> list[dict[str, str]]:
    """Alphabetical, for display."""
    return [{"term": term, "meaning": meaning, "caution": caution}
            for term, meaning, caution in sorted(GLOSSARY,
                                                 key=lambda e: e[0].lower())]


def as_text() -> str:
    """The same list as plain text, for export beside the references."""
    lines = ["GLOSSARY", "=" * 70, "",
             "Plain-English meanings of the words used in this application.",
             "Each entry ends with the part that is most often got wrong.", ""]
    for entry in entries():
        lines.append(entry["term"].upper())
        lines += ["    " + chunk for chunk in _wrap(entry["meaning"], 66)]
        lines.append("    Watch out:")
        lines += ["      " + chunk for chunk in _wrap(entry["caution"], 64)]
        lines.append("")
    return "\n".join(lines)


def _wrap(text: str, width: int) -> list[str]:
    import textwrap
    return textwrap.wrap(text, width) or [""]

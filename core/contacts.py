"""
Which amino acids the docked molecule actually touches.

This is the part of a docking result that is worth trusting. The score
is an estimate with an error bar wider than the range it is trying to
measure; the pose, when the setup has been validated, is usually right
to within about 2 angstroms. So "the molecule sits against Met793 and
Thr854" is a statement a student can reason about, take to the
literature, and test by mutagenesis. "-8.4 kcal/mol" is not.

WHAT THIS IS NOT. Everything here is geometry: how far one atom is
from another. A real hydrogen bond needs the donor, the hydrogen and
the acceptor to be roughly in line, and the hydrogen positions a
docking program writes are not reliable enough to judge that. So a
nitrogen or oxygen close to another nitrogen or oxygen is reported as
a POSSIBLE hydrogen bond, never as a hydrogen bond, and the exported
ChimeraX script runs that program's own `hbonds` command, which does
check the angles. Under-claiming here is deliberate: a student who is
taught that a distance is a bond has been taught something false.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# Any heavy atom of the ligand within this of any heavy atom of a
# residue makes that residue a contact. 4.5 A is the usual convention
# for a van der Waals contact between carbons.
CONTACT = 4.5
# Nitrogen or oxygen pairs this close are geometrically capable of a
# hydrogen bond. The angle is not checked -- see the module note.
HBOND = 3.5
# A charged group pair this close is geometrically capable of a salt
# bridge.
SALT = 4.0
# A metal this close to a ligand heteroatom is being coordinated.
METAL_BOND = 3.0

POLAR = {"N", "O"}
METALS = {"ZN", "MG", "MN", "FE", "CA", "CU", "NI", "CO", "NA", "K"}
WATER = {"HOH", "WAT", "DOD", "H2O"}

# Side-chain atoms that carry a formal charge at pH 7.4.
ANIONIC = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}
CATIONIC = {("LYS", "NZ"), ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"),
            ("HIS", "ND1"), ("HIS", "NE2")}

# Plain-English names, because "ILE" means nothing to a second-year
# student and the point of the table is that they can read it.
RESIDUE_NAMES = {
    "ALA": "alanine", "ARG": "arginine", "ASN": "asparagine",
    "ASP": "aspartate", "CYS": "cysteine", "GLN": "glutamine",
    "GLU": "glutamate", "GLY": "glycine", "HIS": "histidine",
    "ILE": "isoleucine", "LEU": "leucine", "LYS": "lysine",
    "MET": "methionine", "PHE": "phenylalanine", "PRO": "proline",
    "SER": "serine", "THR": "threonine", "TRP": "tryptophan",
    "TYR": "tyrosine", "VAL": "valine",
    "HIE": "histidine", "HID": "histidine", "HIP": "histidine",
    "CYX": "cysteine", "SEC": "selenocysteine", "MSE": "selenomethionine",
}

HYDROPHOBIC_RESIDUES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP",
                        "PRO", "TYR", "CYS", "MSE", "CYX"}

# What each kind of contact is, said once, in the words the interface
# and the exported file both use.
KIND_TEXT = {
    "hbond": ("possible hydrogen bond",
              "A nitrogen or oxygen of the drug is close enough to one of "
              "the protein to form a hydrogen bond. The angle has not been "
              "checked here, so this is a possibility, not a measurement."),
    "salt": ("possible salt bridge",
             "A charged group of the drug sits against an oppositely "
             "charged side chain. AutoDock Vina has no electrostatic term, "
             "so this interaction is barely counted in the score."),
    "metal": ("metal coordination",
              "A nitrogen, oxygen or sulphur of the drug is close to a "
              "metal ion in the protein. In a metalloenzyme this is often "
              "the whole basis of binding."),
    "hydrophobic": ("hydrophobic contact",
                    "An oily part of the drug is resting against an oily "
                    "part of the protein. This is what most of a docking "
                    "score is made of."),
    "contact": ("close contact",
                "The atoms are touching, without falling into one of the "
                "categories above."),
}

# Most informative first, so a residue making both a hydrogen bond and
# a hydrophobic contact is listed under the hydrogen bond.
PRIORITY = ["metal", "salt", "hbond", "hydrophobic", "contact"]


@dataclass
class Residue:
    """One amino acid (or cofactor, or metal) touching the ligand."""
    chain: str
    number: str
    code: str                      # three-letter code, e.g. MET
    distance: float                # closest heavy-atom approach, A
    kinds: set = field(default_factory=set)
    atoms: int = 0                 # how many atom pairs are in contact
    nearest: tuple = ("", "")      # (protein atom, ligand atom)

    @property
    def name(self) -> str:
        # A cofactor or metal keeps its own code: "HEM" and "ZN" mean
        # something, "Hem" and "Zn" read like typing mistakes.
        return RESIDUE_NAMES.get(self.code.upper(), self.code.upper())

    @property
    def label(self) -> str:
        """What a paper would call it: Met793."""
        pretty = self.code.title() if self.code.isalpha() else self.code
        return f"{pretty}{self.number}"

    @property
    def kind(self) -> str:
        for key in PRIORITY:
            if key in self.kinds:
                return key
        return "contact"

    @property
    def kind_text(self) -> str:
        return KIND_TEXT[self.kind][0]


def _atoms(pdb_text: str, want_hetero: bool | None = None,
           exclude: str | list[str] = "") -> list[dict]:
    """Heavy atoms from PDB text. Hydrogens and water are left out."""
    skip = {exclude} if isinstance(exclude, str) else set(exclude)
    skip = {s.strip().upper() for s in skip if s and s.strip()}
    out = []
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if want_hetero is True and not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip().upper()
        if resname in WATER:
            continue
        # The project records a het as "AQ4 A999": name, then chain
        # and number run together.
        if skip and f"{resname} {line[21]}{line[22:27].strip()}".upper() in skip:
            continue
        element = (line[76:78].strip() or line[12:16].strip()[:1]).upper()
        if element in ("H", "D"):
            continue
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        out.append({"name": line[12:16].strip(), "resname": resname,
                    "chain": line[21], "resseq": line[22:27].strip(),
                    "element": element, "xyz": xyz})
    return out


def _classify(protein: dict, ligand: dict, distance: float) -> str | None:
    p_el, l_el = protein["element"], ligand["element"]
    key = (protein["resname"], protein["name"])

    if p_el in METALS and l_el in POLAR | {"S"} and distance <= METAL_BOND:
        return "metal"
    if l_el in METALS and p_el in POLAR | {"S"} and distance <= METAL_BOND:
        return "metal"
    if distance <= SALT and (key in ANIONIC or key in CATIONIC) \
            and l_el in POLAR:
        return "salt"
    if distance <= HBOND and p_el in POLAR and l_el in POLAR:
        return "hbond"
    if distance <= CONTACT and p_el in ("C", "S") and l_el in ("C", "S", "CL",
                                                               "BR", "I", "F"):
        return "hydrophobic"
    if distance <= CONTACT:
        return "contact"
    return None


def find(receptor_pdb: str, ligand_pdb: str,
         cutoff: float = CONTACT, exclude: str | list[str] = "") -> list[Residue]:
    """
    Residues of the receptor within `cutoff` of the docked molecule.

    `exclude` takes labels in the project's own form, "AQ4 A999". The
    crystal ligand belongs there: a holo structure keeps it, the docked
    pose lands on top of it, and without this the molecule is reported
    as being in contact with itself at 0.3 angstroms -- which is true,
    and useless, and would be the first row of the table.

    Returned closest first, which is usually the order of how much the
    contact matters.
    """
    protein = _atoms(receptor_pdb, exclude=exclude)
    ligand = _atoms(ligand_pdb, want_hetero=True) or _atoms(ligand_pdb)
    if not protein or not ligand:
        return []

    p_xyz = np.array([a["xyz"] for a in protein], dtype=float)
    l_xyz = np.array([a["xyz"] for a in ligand], dtype=float)

    # Throw away everything far from the ligand's bounding box before
    # the full distance calculation: a receptor is thousands of atoms
    # and all but a few dozen are nowhere near.
    low = l_xyz.min(axis=0) - cutoff
    high = l_xyz.max(axis=0) + cutoff
    near = np.all((p_xyz >= low) & (p_xyz <= high), axis=1)
    if not near.any():
        return []
    indices = np.flatnonzero(near)
    distances = np.linalg.norm(
        p_xyz[indices][:, None, :] - l_xyz[None, :, :], axis=2)

    found: dict[tuple, Residue] = {}
    for row, p_index in enumerate(indices):
        close = np.flatnonzero(distances[row] <= cutoff)
        if close.size == 0:
            continue
        p_atom = protein[p_index]
        key = (p_atom["chain"], p_atom["resseq"], p_atom["resname"])
        residue = found.get(key)
        if residue is None:
            residue = Residue(chain=p_atom["chain"], number=p_atom["resseq"],
                              code=p_atom["resname"], distance=float("inf"))
            found[key] = residue
        for l_index in close:
            separation = float(distances[row][l_index])
            kind = _classify(p_atom, ligand[l_index], separation)
            if kind is None:
                continue
            residue.kinds.add(kind)
            residue.atoms += 1
            if separation < residue.distance:
                residue.distance = separation
                residue.nearest = (p_atom["name"], ligand[l_index]["name"])

    residues = [r for r in found.values() if r.atoms]
    for residue in residues:
        residue.distance = round(residue.distance, 2)
        # A whole hydrophobic side chain resting on the ligand is a
        # hydrophobic contact even when its backbone oxygen happens to
        # be the nearest atom; one stray polar pair should not relabel
        # the residue.
        if residue.kinds == {"contact"} \
                and residue.code.upper() in HYDROPHOBIC_RESIDUES:
            residue.kinds = {"hydrophobic"}
    residues.sort(key=lambda r: (r.distance, r.number))
    return residues


def rows(residues: list[Residue]) -> list[dict]:
    """Table rows for the interface."""
    return [{"residue": r.label, "name": r.name, "chain": r.chain,
             "kind": r.kind_text, "distance": f"{r.distance:.2f}",
             "atoms": r.atoms} for r in residues]


def summary(residues: list[Residue]) -> str:
    """One sentence a student can put in a write-up."""
    if not residues:
        return ("No residue is within 4.5 angstroms of the molecule. Either "
                "the pose is outside the pocket or the wrong files were "
                "compared.")
    bonds = [r for r in residues if r.kind in ("hbond", "salt", "metal")]
    text = (f"The molecule touches {len(residues)} residue"
            f"{'s' if len(residues) != 1 else ''}, closest "
            f"{residues[0].label} at {residues[0].distance:.2f} angstroms.")
    if bonds:
        named = ", ".join(r.label for r in bonds[:4])
        text += (f" Polar contacts close enough for a hydrogen bond or salt "
                 f"bridge: {named}"
                 f"{' and others' if len(bonds) > 4 else ''}.")
    else:
        text += (" No polar contact is close enough for a hydrogen bond; "
                 "this pose is held by shape and hydrophobic contact alone.")
    return text


def report(residues: list[Residue], target: str = "", drug: str = "",
           score: float | None = None) -> str:
    """The contacts.txt that ships in the export."""
    width = 74
    lines = ["Contacts between the docked molecule and the protein",
             "=" * width, ""]
    if target or drug:
        lines += [f"Protein : {target}", f"Drug    : {drug}"]
    if score is not None:
        lines.append(f"Score   : {score} kcal/mol  "
                     "(read it with the error bar; see Docking_report.txt)")
    lines += ["", summary(residues), "",
              f"Any residue with a heavy atom within {CONTACT} angstroms of a "
              "heavy atom of the", "drug is listed. Closest first.", "",
              f"{'Residue':<12}{'Amino acid':<18}{'Chain':<7}"
              f"{'Closest':<10}{'Contact type':<28}",
              "-" * width]
    for residue in residues:
        lines.append(f"{residue.label:<12}{residue.name:<18}"
                     f"{residue.chain:<7}{residue.distance:<10.2f}"
                     f"{residue.kind_text:<28}")
    lines += ["", "What these mean", "-" * width]
    for key in PRIORITY:
        if any(key in r.kinds for r in residues):
            title, explanation = KIND_TEXT[key]
            lines.append(f"{title}:")
            for chunk in _wrap(explanation, width - 4):
                lines.append("    " + chunk)
            lines.append("")
    lines += [
        "A caution worth reading twice",
        "-" * width,
        *_wrap(
            "Everything above is distance. A hydrogen bond also needs the "
            "donor, hydrogen and acceptor to be roughly in a line, and the "
            "hydrogen positions written by a docking program are not good "
            "enough to judge that. Open view_in_chimerax.cxc to see the "
            "hydrogen bonds worked out properly, with angles, by ChimeraX.",
            width),
        "",
        *_wrap(
            "These contacts are only as good as the pose. If the re-docking "
            "check in Docking_report.txt did not pass, the pose has not been "
            "shown to be reliable and neither has this list.", width),
        "",
    ]
    return "\n".join(lines)


def _wrap(text: str, width: int) -> list[str]:
    import textwrap
    return textwrap.wrap(text, width) or [""]

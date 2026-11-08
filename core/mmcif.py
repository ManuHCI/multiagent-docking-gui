"""
Read an mmCIF coordinate file and write it out as PDB.

Why this exists
---------------
The Protein Data Bank is in the middle of a long migration. mmCIF is
the archive's master format; the legacy .pdb file is a convenience
copy, and for a growing number of entries -- particularly recent
cryo-EM structures -- that copy is simply not generated. Asking
files.rcsb.org for such an entry as .pdb returns a 404.

The first version of the structure step did not know this. It asked
for .pdb, got nothing, and told the user to check their internet
connection. The internet was fine. The file did not exist. That is the
worst kind of error message: it is confident, it is specific, and it
sends the user to look in the wrong place.

So this module gives the software a second thing to try. If the legacy
file is missing, fetch the mmCIF instead and convert it here.

Why not use OpenMM's PDBxFile
-----------------------------
It is already installed and it would work for most entries. But it
rebuilds the topology to do it, which means it drops anything it does
not recognise, renumbers in some cases, and refuses outright on files
with unusual residues. Everything downstream in this software --
chain selection, the cofactor survey, the PDB coverage check -- reads
author chain IDs and author residue numbers and expects them to be
exactly what the depositor wrote. A direct column-to-column
transcription preserves that; a round trip through a topology does
not.

What it does not do
-------------------
It writes the first model only, which is what every structure-based
calculation in this software wants, and it refuses rather than
corrupting when a structure genuinely cannot be held in PDB format
(more than 99,999 atoms, or chain IDs longer than one character).
Those are real limits of the old format, not bugs, and the honest
response is to say so and suggest another entry.
"""

from __future__ import annotations

from typing import Iterator

# Everything the PDB ATOM record needs, in the order it is written,
# with the mmCIF tag that holds it. The auth_ tags are the depositor's
# own numbering -- the numbering in every paper about the structure --
# and are preferred over the label_ tags, which are a separate
# internal scheme that will not match anything the user reads.
PREFERRED = {
    "atom_name": ("auth_atom_id", "label_atom_id"),
    "alt_loc": ("label_alt_id",),
    "res_name": ("auth_comp_id", "label_comp_id"),
    "chain": ("auth_asym_id", "label_asym_id"),
    "res_seq": ("auth_seq_id", "label_seq_id"),
    "ins_code": ("pdbx_PDB_ins_code",),
    "x": ("Cartn_x",),
    "y": ("Cartn_y",),
    "z": ("Cartn_z",),
    "occupancy": ("occupancy",),
    "b_factor": ("B_iso_or_equiv",),
    "element": ("type_symbol",),
    "charge": ("pdbx_formal_charge",),
    "group": ("group_PDB",),
    "model": ("pdbx_PDB_model_num",),
}


class CifError(Exception):
    """The file cannot be converted, with a reason a user can act on."""


def _tokenise(line: str) -> list[str]:
    """
    Split one row of an mmCIF loop into its values.

    Values may be bare, single-quoted or double-quoted, and quoted ones
    contain spaces -- "O5'" and 'HIS A' both occur in real files. A
    plain .split() silently shifts every column after the first quoted
    value, which puts coordinates in the occupancy field and is not
    detectable downstream.
    """
    out: list[str] = []
    i = 0
    n = len(line)
    while i < n:
        while i < n and line[i] in " \t":
            i += 1
        if i >= n:
            break
        if line[i] in "'\"":
            quote = line[i]
            i += 1
            start = i
            while i < n and not (line[i] == quote
                                 and (i + 1 >= n or line[i + 1] in " \t")):
                i += 1
            out.append(line[start:i])
            i += 1
        else:
            start = i
            while i < n and line[i] not in " \t":
                i += 1
            out.append(line[start:i])
    return out


def _atom_site_rows(text: str) -> Iterator[tuple[dict[str, int], list[str]]]:
    """Walk the _atom_site loop, yielding the column map and each row."""
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        if lines[i].strip() != "loop_":
            i += 1
            continue
        # Collect the tag names that follow 'loop_'.
        j = i + 1
        tags: list[str] = []
        while j < len(lines) and lines[j].lstrip().startswith("_"):
            tags.append(lines[j].strip().split()[0])
            j += 1
        if not tags or not tags[0].startswith("_atom_site."):
            i = j
            continue
        columns = {tag.split(".", 1)[1]: k for k, tag in enumerate(tags)}
        while j < len(lines):
            row = lines[j]
            stripped = row.strip()
            if not stripped or stripped.startswith("#") \
                    or stripped.startswith("loop_") \
                    or stripped.startswith("_"):
                break
            yield columns, _tokenise(row)
            j += 1
        return


def _value(columns: dict[str, int], row: list[str], field: str) -> str:
    for tag in PREFERRED[field]:
        k = columns.get(tag)
        if k is None or k >= len(row):
            continue
        value = row[k]
        if value in (".", "?"):
            return ""
        return value
    return ""


def cif_to_pdb(text: str) -> str:
    """
    Convert mmCIF coordinates to a PDB file, first model only.

    Raises CifError with a sentence the user can act on when the
    structure genuinely will not fit the old format.
    """
    out: list[str] = []
    serial = 0
    first_model = ""
    wide_chains: set[str] = set()
    last_chain = ""
    saw_any = False

    for columns, row in _atom_site_rows(text):
        saw_any = True
        model = _value(columns, row, "model")
        if not first_model:
            first_model = model
        if model and model != first_model:
            break

        chain = _value(columns, row, "chain") or "A"
        if len(chain) > 1:
            wide_chains.add(chain)
            continue

        if last_chain and chain != last_chain:
            out.append(f"TER   {serial + 1:5d}")
            serial += 1
        last_chain = chain

        serial += 1
        if serial > 99_999:
            raise CifError(
                "This structure has more than 99,999 atoms, which the "
                "older PDB file format cannot hold. That is a limit of "
                "the format, not a problem with the download. Pick a "
                "smaller structure of the same protein - usually a "
                "single domain rather than a whole complex."
            )

        element = (_value(columns, row, "element") or "C").strip()
        raw_name = _value(columns, row, "atom_name")
        # PDB convention: a one-letter element starts in column 14.
        name = raw_name.ljust(4) if len(raw_name) >= 4 else f" {raw_name:<3}"

        try:
            x = float(_value(columns, row, "x"))
            y = float(_value(columns, row, "y"))
            z = float(_value(columns, row, "z"))
        except ValueError:
            continue

        try:
            occupancy = float(_value(columns, row, "occupancy") or 1.0)
        except ValueError:
            occupancy = 1.0
        try:
            b_factor = float(_value(columns, row, "b_factor") or 0.0)
        except ValueError:
            b_factor = 0.0

        res_seq_text = _value(columns, row, "res_seq") or "1"
        try:
            res_seq = int(res_seq_text)
        except ValueError:
            continue

        record = _value(columns, row, "group") or "ATOM"
        # mmCIF writes the formal charge as a signed integer ("2",
        # "-1"); PDB writes it sign-last ("2+", "1-"). Some files
        # already use the PDB form, so both are accepted and the field
        # is dropped entirely if it will not fit the two columns --
        # a third character there runs into the next record.
        charge = _value(columns, row, "charge").strip()
        digits = charge.strip("+-")
        charge_field = ""
        if digits.isdigit() and digits != "0":
            sign = "-" if "-" in charge else "+"
            charge_field = f"{digits}{sign}" if len(digits) == 1 else ""

        out.append(
            f"{record:<6}{serial:5d} {name}"
            f"{_value(columns, row, 'alt_loc')[:1] or ' '}"
            f"{_value(columns, row, 'res_name')[:3]:>3} "
            f"{chain}{res_seq:4d}"
            f"{_value(columns, row, 'ins_code')[:1] or ' '}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
            f"{occupancy:6.2f}{b_factor:6.2f}          "
            f"{element[:2]:>2}{charge_field:<2}"
        )

    if not saw_any:
        raise CifError(
            "No coordinates were found in the file that was downloaded. "
            "It may have been truncated. Try running this step again."
        )
    if not out:
        if wide_chains:
            raise CifError(
                "Every chain in this structure has a label longer than "
                "one letter, which the older PDB file format cannot "
                "hold. This happens with very large assemblies. Pick a "
                "smaller structure of the same protein."
            )
        raise CifError(
            "The coordinates in this file could not be read. Pick a "
            "different structure."
        )

    out.append("END")
    return "\n".join(out) + "\n"


def dropped_chains(text: str) -> list[str]:
    """Chains a conversion had to leave out, so the user can be told."""
    wide: set[str] = set()
    first_model = ""
    for columns, row in _atom_site_rows(text):
        model = _value(columns, row, "model")
        if not first_model:
            first_model = model
        if model and model != first_model:
            break
        chain = _value(columns, row, "chain")
        if len(chain) > 1:
            wide.add(chain)
    return sorted(wide)

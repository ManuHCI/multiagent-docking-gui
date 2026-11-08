"""Source identity, resolver errors and reference-export ownership checks."""
from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from urllib.parse import parse_qs, urlsplit

from core.contract import ProvenanceEntry
from core.project import Project
from core.references import (DOWNLOAD_LINKS, drug_links, references,
                             render_references, resolve_drug, structure_links)
from core.web import Fetched


class ReferenceLinkTests(unittest.TestCase):
    def test_exact_identity_and_reserved_characters_survive_download_links(self):
        protein = structure_links("1m17")
        self.assertEqual(protein[1]["url"], "https://files.rcsb.org/download/1M17.cif")
        self.assertEqual(structure_links("../1M17?other=1"), [])
        drug = drug_links("Display label", "C/C=C\\C.[Na+]")
        query = parse_qs(urlsplit(drug[1]["url"]).query)
        self.assertEqual(query["smiles"], ["C/C=C\\C.[Na+]"])
        self.assertEqual(query["record_type"], ["2d"])
        cid = drug_links("An unrelated label", "CCO", 2244)
        self.assertEqual(cid[0]["url"], "https://pubchem.ncbi.nlm.nih.gov/compound/2244")
        self.assertIn("/cid/2244/SDF?record_type=3d", cid[2]["url"])

    def test_catalogue_has_required_sources_and_returns_independent_records(self):
        first = references()
        second = references()
        ids = {record["id"] for record in first}
        self.assertEqual(len(ids), len(first))
        self.assertTrue({"vina12", "rdkit", "meeko", "pdbfixer", "pdb", "pubchem", "lipinski", "veber", "3dmol"}.issubset(ids))
        self.assertTrue(all(record["url"].startswith("https://") and record["citation"] for record in first))
        first[0]["citation"] = "Changed by a UI consumer"
        self.assertNotEqual(first[0]["citation"], second[0]["citation"])
        self.assertTrue(all(record["title"] and record["description"] for record in DOWNLOAD_LINKS))


class PubChemResolutionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="reference_lookup_")
        self.cache = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    @staticmethod
    def response(properties):
        return Fetched(ok=True, url="https://pubchem.ncbi.nlm.nih.gov/test",
                       text=json.dumps({"PropertyTable": {"Properties": properties}}),
                       retrieved="2026-11-07T12:00:00+00:00", from_cache=True)

    def test_current_pubchem_stereochemical_smiles_and_source_are_preserved(self):
        stereo = "N[C@@H](C)C(=O)O"
        answer = self.response([{"CID": 5950, "Title": "L-alanine", "SMILES": stereo,
                                 "ConnectivitySMILES": "CC(C(=O)O)N"}])
        with patch("core.web.Web.get", return_value=answer) as request:
            drug = resolve_drug("CID: 5950", self.cache)
        self.assertEqual(drug["smiles"], stereo)
        self.assertEqual(drug["cid"], 5950)
        self.assertTrue(drug["from_cache"])
        self.assertEqual(drug["retrieved"], answer.retrieved)
        self.assertEqual(drug["source_url"], request.call_args.args[0])
        self.assertIn("/cid/5950/property/SMILES,Title/JSON", drug["source_url"])

    def test_cached_older_property_names_keep_isomeric_structure(self):
        stereo = "C[C@H](O)C(=O)O"
        answer = self.response([{"CID": 107689, "Title": "Lactic acid", "IsomericSMILES": stereo,
                                 "CanonicalSMILES": "CC(O)C(=O)O"}])
        with patch("core.web.Web.get", return_value=answer):
            drug = resolve_drug("lactic acid", self.cache)
        self.assertEqual(drug["smiles"], stereo)
        self.assertIn("/name/lactic%20acid/", drug["source_url"])

    def test_unavailable_or_ambiguous_data_does_not_invent_a_compound(self):
        offline = Fetched(ok=False, url="test", problem="PubChem could not be reached.",
                          suggestion="Enter a verified SMILES structure.")
        with patch("core.web.Web.get", return_value=offline):
            with self.assertRaisesRegex(ValueError, "could not be reached"):
                resolve_drug("aspirin", self.cache)
        ambiguous = self.response([{"CID": 1, "SMILES": "CC"}, {"CID": 2, "SMILES": "CCC"}])
        with patch("core.web.Web.get", return_value=ambiguous):
            with self.assertRaisesRegex(ValueError, "several compound records"):
                resolve_drug("ambiguous", self.cache)
        broken = self.response([{"CID": "invalid", "SMILES": "CC"}])
        with patch("core.web.Web.get", return_value=broken):
            with self.assertRaisesRegex(ValueError, "no usable compound structure"):
                resolve_drug("invalid data", self.cache)


class RunReferenceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="run_references_")
        self.folder = Path(self.temporary.name) / "Target one"
        self.project = Project.create(self.folder, "Target one")
        self.source_url = "https://files.rcsb.org/download/1M17.cif"
        self.project.add_provenance(ProvenanceEntry(
            stage="02_structure", tool="RCSB PDB", version="source release",
            parameters={"pdb_id": "1M17", "source_url": self.source_url},
            seed=23, source_release="test release", citation="Recorded structure citation"))

    def tearDown(self):
        self.project.close()
        self.temporary.cleanup()

    def test_export_keeps_actual_settings_sources_and_borrowed_project_open(self):
        with patch.object(self.project, "close", wraps=self.project.close) as close:
            text = render_references(self.project)
            close.assert_not_called()
        self.assertIn("PROJECT: Target one", text)
        self.assertIn(self.source_url, text)
        self.assertIn("Random seed: 23", text)
        self.assertIn("Source release: test release", text)
        self.assertIn("Recorded structure citation", text)
        self.assertEqual(len(self.project.provenance()), 1)

    def test_folder_export_closes_only_the_handle_it_opens(self):
        owned = Project.open(self.folder)
        with patch("core.project.Project.open", return_value=owned), patch.object(owned, "close", wraps=owned.close) as close:
            text = render_references([self.folder])
            close.assert_called_once()
        self.assertIn("Recorded structure citation", text)
        self.assertEqual(len(self.project.provenance()), 1)


if __name__ == "__main__":
    unittest.main()

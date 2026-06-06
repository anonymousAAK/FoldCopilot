"""Tests for AlphaMissense and AlphaFill annotation tools."""


from foldcopilot.clients.alphamissense_client import (
    _classify_score,
    _compute_fraction,
    _extract_position,
)
from foldcopilot.clients.alphafill_client import _classify_compound
from foldcopilot.tools.annotations import _find_cofactor_pathogenicity_hotspots


class TestAlphaMissenseHelpers:
    def test_extract_position_standard(self):
        assert _extract_position("M1A") == 1
        assert _extract_position("G123R") == 123
        assert _extract_position("W456F") == 456

    def test_extract_position_edge_cases(self):
        assert _extract_position("") is None
        assert _extract_position("ABC") is None

    def test_classify_pathogenic(self):
        assert _classify_score(0.9) == "likely_pathogenic"
        assert _classify_score(0.564) == "likely_pathogenic"

    def test_classify_benign(self):
        assert _classify_score(0.1) == "likely_benign"
        assert _classify_score(0.34) == "likely_benign"

    def test_classify_ambiguous(self):
        assert _classify_score(0.45) == "ambiguous"
        assert _classify_score(0.5) == "ambiguous"

    def test_compute_fraction_empty(self):
        assert _compute_fraction([], "pathogenic") == 0.0

    def test_compute_fraction_all(self):
        variants = [
            {"classification": "likely_pathogenic"},
            {"classification": "likely_pathogenic"},
        ]
        assert _compute_fraction(variants, "pathogenic") == 1.0

    def test_compute_fraction_mixed(self):
        variants = [
            {"classification": "likely_pathogenic"},
            {"classification": "likely_benign"},
            {"classification": "ambiguous"},
            {"classification": "likely_pathogenic"},
        ]
        assert _compute_fraction(variants, "pathogenic") == 0.5


class TestAlphaFillHelpers:
    def test_classify_cofactor(self):
        assert _classify_compound({"compound_id": "NAD"}) == "cofactor"
        assert _classify_compound({"compound_id": "FAD"}) == "cofactor"
        assert _classify_compound({"compound_id": "HEM"}) == "cofactor"

    def test_classify_metal(self):
        assert _classify_compound({"compound_id": "ZN"}) == "metal_ion"
        assert _classify_compound({"compound_id": "MG"}) == "metal_ion"
        assert _classify_compound({"compound_id": "CA"}) == "metal_ion"

    def test_classify_nucleotide(self):
        assert _classify_compound({"compound_id": "ATP"}) == "nucleotide"
        assert _classify_compound({"compound_id": "GTP"}) == "nucleotide"

    def test_classify_solvent(self):
        assert _classify_compound({"compound_id": "HOH"}) == "solvent"

    def test_classify_ligand_default(self):
        assert _classify_compound({"compound_id": "XYZ"}) == "ligand"

    def test_classify_by_name(self):
        assert _classify_compound(
            {"compound_id": "ABC", "compound_name": "kinase inhibitor"}
        ) == "ligand"


class TestCofactorPathogenicityHotspots:
    def test_no_data(self):
        result = _find_cofactor_pathogenicity_hotspots(
            {"available": False}, {"available": False}
        )
        assert result == []

    def test_no_transplants(self):
        missense = {
            "available": True,
            "residue_landscape": {"1": {"classification": "likely_pathogenic", "mean_pathogenicity": 0.9}},
        }
        cofactors = {"available": True, "transplants": []}
        result = _find_cofactor_pathogenicity_hotspots(missense, cofactors)
        assert result == []

    def test_finds_hotspots(self):
        missense = {
            "available": True,
            "residue_landscape": {
                "10": {"classification": "likely_pathogenic", "mean_pathogenicity": 0.95},
                "20": {"classification": "likely_benign", "mean_pathogenicity": 0.1},
                "30": {"classification": "likely_pathogenic", "mean_pathogenicity": 0.8},
            },
        }
        cofactors = {
            "available": True,
            "transplants": [{"compound_id": "ATP", "source_pdb": "1ATP"}],
        }
        result = _find_cofactor_pathogenicity_hotspots(missense, cofactors)
        assert len(result) == 2  # positions 10 and 30
        # Sorted by pathogenicity, descending
        assert result[0]["residue_position"] == 10
        assert result[0]["mean_pathogenicity"] == 0.95

    def test_max_20_hotspots(self):
        missense = {
            "available": True,
            "residue_landscape": {
                str(i): {"classification": "likely_pathogenic", "mean_pathogenicity": 0.9}
                for i in range(50)
            },
        }
        cofactors = {
            "available": True,
            "transplants": [{"compound_id": "NAD"}],
        }
        result = _find_cofactor_pathogenicity_hotspots(missense, cofactors)
        assert len(result) <= 20

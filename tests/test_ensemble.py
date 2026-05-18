"""Tests for ensemble comparison and cross-model disagreement detection."""

import numpy as np
import pytest

from foldcopilot.models.ensemble import AgreementLevel, EnsembleReport
from foldcopilot.tools.ensemble import (
    classify_agreement,
    compare_structures,
    compute_ca_distances,
    compute_rmsd,
    find_disagreement_spans,
    parse_ca_coords,
    parse_plddt_from_pdb,
    ResidueDisagreement,
)


# Helper: generate a minimal PDB line for a CA atom
def _pdb_ca_line(
    serial: int, res_name: str, res_seq: int, x: float, y: float, z: float, b: float
) -> str:
    return (
        f"ATOM  {serial:5d}  CA  {res_name:3s} A{res_seq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{b:6.2f}           C"
    )


def _make_pdb(residues: list[tuple[str, int, float, float, float, float]]) -> str:
    """Create minimal PDB text. Each tuple: (res_name, res_seq, x, y, z, bfactor)."""
    lines = []
    for i, (name, seq, x, y, z, b) in enumerate(residues, 1):
        lines.append(_pdb_ca_line(i, name, seq, x, y, z, b))
    lines.append("END")
    return "\n".join(lines)


class TestParseCaCoords:
    def test_basic(self):
        pdb = _make_pdb([("ALA", 1, 1.0, 2.0, 3.0, 85.0)])
        coords = parse_ca_coords(pdb)
        assert len(coords) == 1
        assert coords[0][0] == 1
        assert coords[0][1] == "ALA"
        np.testing.assert_array_almost_equal(coords[0][2], [1.0, 2.0, 3.0])

    def test_multiple_residues(self):
        pdb = _make_pdb([
            ("ALA", 1, 1.0, 0.0, 0.0, 90.0),
            ("GLY", 2, 4.0, 0.0, 0.0, 80.0),
            ("VAL", 3, 7.0, 0.0, 0.0, 50.0),
        ])
        coords = parse_ca_coords(pdb)
        assert len(coords) == 3


class TestParsePlddt:
    def test_basic(self):
        pdb = _make_pdb([
            ("ALA", 1, 0.0, 0.0, 0.0, 85.5),
            ("GLY", 2, 0.0, 0.0, 0.0, 92.3),
        ])
        plddt = parse_plddt_from_pdb(pdb)
        assert plddt[1] == 85.5
        assert plddt[2] == 92.3


class TestCaDistances:
    def test_identical(self):
        coords = [(1, "ALA", np.array([0.0, 0.0, 0.0]))]
        distances = compute_ca_distances(coords, coords)
        assert distances[1] == 0.0

    def test_known_distance(self):
        a = [(1, "ALA", np.array([0.0, 0.0, 0.0]))]
        b = [(1, "ALA", np.array([3.0, 4.0, 0.0]))]
        distances = compute_ca_distances(a, b)
        assert abs(distances[1] - 5.0) < 0.001

    def test_missing_residue(self):
        a = [(1, "ALA", np.array([0.0, 0.0, 0.0])), (2, "GLY", np.array([1.0, 0.0, 0.0]))]
        b = [(1, "ALA", np.array([0.0, 0.0, 0.0]))]  # missing residue 2
        distances = compute_ca_distances(a, b)
        assert 1 in distances
        assert 2 not in distances


class TestRMSD:
    def test_zero(self):
        assert compute_rmsd({1: 0.0, 2: 0.0}) == 0.0

    def test_known(self):
        # RMSD of [3, 4] = sqrt((9+16)/2) = sqrt(12.5) ~ 3.536
        rmsd = compute_rmsd({1: 3.0, 2: 4.0})
        assert abs(rmsd - 3.536) < 0.01

    def test_empty(self):
        assert compute_rmsd({}) == float("inf")


class TestClassifyAgreement:
    def test_strong_agree(self):
        result = classify_agreement(90.0, 85.0, ca_distance=1.5)
        assert result == AgreementLevel.STRONG_AGREE

    def test_disagree(self):
        result = classify_agreement(90.0, 85.0, ca_distance=5.0)
        assert result == AgreementLevel.DISAGREE

    def test_both_uncertain(self):
        result = classify_agreement(40.0, 35.0, ca_distance=1.0)
        assert result == AgreementLevel.BOTH_UNCERTAIN

    def test_moderate_agree(self):
        # One confident, one not, structures close
        result = classify_agreement(85.0, 50.0, ca_distance=2.0)
        assert result == AgreementLevel.MODERATE_AGREE

    def test_no_distance(self):
        result = classify_agreement(90.0, 85.0, ca_distance=None)
        assert result == AgreementLevel.MODERATE_AGREE

    def test_custom_thresholds(self):
        result = classify_agreement(
            60.0, 65.0, ca_distance=1.0, plddt_threshold=50.0
        )
        assert result == AgreementLevel.STRONG_AGREE


class TestDisagreementSpans:
    def _make_residues(self, agreements: list[AgreementLevel]) -> list[ResidueDisagreement]:
        return [
            ResidueDisagreement(
                residue_index=i + 1,
                model_a="a", model_b="b",
                plddt_a=85.0, plddt_b=80.0,
                ca_distance=1.0 if a == AgreementLevel.STRONG_AGREE else 6.0,
                agreement=a,
            )
            for i, a in enumerate(agreements)
        ]

    def test_no_disagreement(self):
        residues = self._make_residues([AgreementLevel.STRONG_AGREE] * 10)
        spans = find_disagreement_spans(residues, AgreementLevel.DISAGREE)
        assert spans == []

    def test_single_span(self):
        agreements = (
            [AgreementLevel.STRONG_AGREE] * 5
            + [AgreementLevel.DISAGREE] * 4
            + [AgreementLevel.STRONG_AGREE] * 5
        )
        residues = self._make_residues(agreements)
        spans = find_disagreement_spans(residues, AgreementLevel.DISAGREE)
        assert len(spans) == 1
        assert spans[0].start == 6
        assert spans[0].end == 9
        assert spans[0].length == 4

    def test_min_length_filter(self):
        agreements = (
            [AgreementLevel.STRONG_AGREE] * 5
            + [AgreementLevel.DISAGREE] * 2  # too short
            + [AgreementLevel.STRONG_AGREE] * 5
        )
        residues = self._make_residues(agreements)
        spans = find_disagreement_spans(residues, AgreementLevel.DISAGREE, min_length=3)
        assert len(spans) == 0

    def test_span_at_end(self):
        agreements = [AgreementLevel.STRONG_AGREE] * 5 + [AgreementLevel.DISAGREE] * 4
        residues = self._make_residues(agreements)
        spans = find_disagreement_spans(residues, AgreementLevel.DISAGREE)
        assert len(spans) == 1
        assert spans[0].end == 9

    def test_find_agree_spans(self):
        agreements = [AgreementLevel.STRONG_AGREE] * 10
        residues = self._make_residues(agreements)
        spans = find_disagreement_spans(
            residues, AgreementLevel.STRONG_AGREE, min_length=5
        )
        assert len(spans) == 1
        assert spans[0].length == 10


class TestCompareStructures:
    @pytest.mark.asyncio
    async def test_identical_structures(self):
        pdb = _make_pdb([
            ("ALA", 1, 1.0, 0.0, 0.0, 90.0),
            ("GLY", 2, 4.0, 0.0, 0.0, 85.0),
            ("VAL", 3, 7.0, 0.0, 0.0, 80.0),
            ("LEU", 4, 10.0, 0.0, 0.0, 88.0),
            ("ILE", 5, 13.0, 0.0, 0.0, 92.0),
        ])
        result = await compare_structures(pdb, pdb, "boltz2", "af3")
        assert result["mean_ca_rmsd"] == 0.0
        assert result["agreement_fraction"] == 1.0
        assert result["disagree_fraction"] == 0.0
        assert result["strong_agree_fraction"] == 1.0

    @pytest.mark.asyncio
    async def test_different_structures(self):
        pdb_a = _make_pdb([
            ("ALA", 1, 0.0, 0.0, 0.0, 90.0),
            ("GLY", 2, 3.0, 0.0, 0.0, 85.0),
            ("VAL", 3, 6.0, 0.0, 0.0, 80.0),
            ("LEU", 4, 9.0, 0.0, 0.0, 88.0),
            ("ILE", 5, 12.0, 0.0, 0.0, 92.0),
        ])
        # Shift all coords by 10A — clear disagreement
        pdb_b = _make_pdb([
            ("ALA", 1, 10.0, 10.0, 10.0, 90.0),
            ("GLY", 2, 13.0, 10.0, 10.0, 85.0),
            ("VAL", 3, 16.0, 10.0, 10.0, 80.0),
            ("LEU", 4, 19.0, 10.0, 10.0, 88.0),
            ("ILE", 5, 22.0, 10.0, 10.0, 92.0),
        ])
        result = await compare_structures(pdb_a, pdb_b, "boltz2", "af3")
        assert result["disagree_fraction"] > 0.5
        assert result["mean_ca_rmsd"] > 5.0

    @pytest.mark.asyncio
    async def test_empty_pdb(self):
        result = await compare_structures("", "ATOM ...", "a", "b")
        assert "error" in result

    @pytest.mark.asyncio
    async def test_report_has_caveats(self):
        pdb = _make_pdb([
            ("ALA", 1, 1.0, 0.0, 0.0, 90.0),
            ("GLY", 2, 4.0, 0.0, 0.0, 85.0),
            ("VAL", 3, 7.0, 0.0, 0.0, 80.0),
        ])
        result = await compare_structures(pdb, pdb, "a", "b")
        assert len(result["caveats"]) >= 3
        assert any("research use" in c for c in result["caveats"])


class TestEnsembleReport:
    def test_standard_caveats(self):
        report = EnsembleReport(
            models=[],
            residue_count=100,
            agreement_fraction=0.9,
            strong_agree_fraction=0.7,
            moderate_agree_fraction=0.2,
            disagree_fraction=0.05,
            both_uncertain_fraction=0.05,
            interpretation="test",
        )
        report.add_standard_caveats()
        assert len(report.caveats) == 4
        assert any("training data" in c for c in report.caveats)

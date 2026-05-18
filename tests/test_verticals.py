"""Tests for therapeutic vertical packs."""

import pytest

from foldcopilot.tools.verticals import _identify_cdrs


class TestCDRIdentification:
    def test_heavy_chain_cdrs(self):
        # Typical heavy chain ~120 residues
        heavy = "A" * 130
        cdrs = _identify_cdrs(heavy, chain_type="heavy")
        assert len(cdrs) == 3
        names = [c["cdr"] for c in cdrs]
        assert "CDR-H1" in names
        assert "CDR-H2" in names
        assert "CDR-H3" in names

    def test_light_chain_cdrs(self):
        light = "A" * 110
        cdrs = _identify_cdrs(light, chain_type="light")
        assert len(cdrs) == 3
        names = [c["cdr"] for c in cdrs]
        assert "CDR-L1" in names
        assert "CDR-L2" in names
        assert "CDR-L3" in names

    def test_short_sequence_no_cdrs(self):
        short = "A" * 50
        cdrs = _identify_cdrs(short, chain_type="heavy")
        assert len(cdrs) == 0

    def test_cdr_h3_warning(self):
        heavy = "A" * 130
        cdrs = _identify_cdrs(heavy, chain_type="heavy")
        h3 = [c for c in cdrs if c["cdr"] == "CDR-H3"][0]
        assert "variable" in h3["note"].lower() or "low plddt" in h3["note"].lower()

    def test_cdr_positions_are_valid(self):
        heavy = "A" * 130
        cdrs = _identify_cdrs(heavy, chain_type="heavy")
        for cdr in cdrs:
            assert cdr["approximate_start"] > 0
            assert cdr["approximate_end"] > cdr["approximate_start"]
            assert cdr["approximate_end"] <= len(heavy)


class TestAntibodyAnalysis:
    @pytest.mark.asyncio
    async def test_nanobody_detection(self):
        from foldcopilot.tools.verticals import antibody_analysis
        result = await antibody_analysis("A" * 130)
        assert result["is_nanobody"] is True
        assert result["light_chain_length"] == 0
        assert result["pack"] == "antibody"

    @pytest.mark.asyncio
    async def test_full_antibody(self):
        from foldcopilot.tools.verticals import antibody_analysis
        result = await antibody_analysis("A" * 130, "G" * 110)
        assert result["is_nanobody"] is False
        assert result["light_chain_length"] == 110
        # Should have both heavy and light CDRs
        assert len(result["predicted_cdrs"]) == 6

    @pytest.mark.asyncio
    async def test_has_warnings(self):
        from foldcopilot.tools.verticals import antibody_analysis
        result = await antibody_analysis("A" * 130)
        assert len(result["cdr_warnings"]) >= 1
        assert len(result["recommendations"]) >= 1

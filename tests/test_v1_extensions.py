"""Tests for v1 extension features: disorder detection, reproducibility manifests,
confidence report fields, and Protenix backend support."""

from __future__ import annotations

import pytest

from foldcopilot.tools.confidence import _extract_af2_disorder_regions
from foldcopilot.models.prediction import (
    BACKEND_LICENSES,
    LicenseType,
    PredictionBackend,
    ReproducibilityManifest,
)
from foldcopilot.models.confidence import ConfidenceReport


# ---------------------------------------------------------------------------
# 1. TestAF2DisorderRegions
# ---------------------------------------------------------------------------

class TestAF2DisorderRegions:
    """Tests for _extract_af2_disorder_regions helper."""

    def test_no_disorder(self):
        scores = [90.0, 85.0, 80.0, 75.0]  # all above 50
        regions = _extract_af2_disorder_regions(scores)
        assert regions == []

    def test_single_disorder_span(self):
        scores = [90.0, 40.0, 35.0, 45.0, 85.0]  # disorder at positions 2-4
        regions = _extract_af2_disorder_regions(scores)
        assert len(regions) == 1
        assert regions[0]["start"] == 2  # 1-indexed
        assert regions[0]["end"] == 4
        assert regions[0]["length"] == 3

    def test_multiple_disorder_spans(self):
        scores = [30.0, 25.0, 80.0, 90.0, 40.0, 35.0]
        regions = _extract_af2_disorder_regions(scores)
        assert len(regions) == 2

    def test_custom_threshold(self):
        scores = [70.0, 65.0, 60.0]  # below 70 threshold
        regions = _extract_af2_disorder_regions(scores, threshold=70.0)
        assert len(regions) == 1


# ---------------------------------------------------------------------------
# 2. TestReproducibilityManifestProvenance
# ---------------------------------------------------------------------------

class TestReproducibilityManifestProvenance:
    """Tests for training_data_source on ReproducibilityManifest."""

    def test_manifest_with_training_data_source(self):
        m = ReproducibilityManifest.create(
            backend="openfold3",
            sequences=["MKFL"],
            parameters={},
            training_data_source="OpenFold3 training data (AWS Open Data)",
        )
        assert m.training_data_source == "OpenFold3 training data (AWS Open Data)"

    def test_manifest_without_training_data_source(self):
        m = ReproducibilityManifest.create(
            backend="boltz2",
            sequences=["MKFL"],
            parameters={},
        )
        assert m.training_data_source is None


# ---------------------------------------------------------------------------
# 3. TestConfidenceReportFields
# ---------------------------------------------------------------------------

class TestConfidenceReportFields:
    """Tests that new fields exist on the ConfidenceReport model."""

    def test_report_has_af2_disorder_field(self):
        assert "af2_disorder_regions" in ConfidenceReport.model_fields

    def test_report_has_idr_strategy_note(self):
        assert "idr_strategy_note" in ConfidenceReport.model_fields

    def test_report_has_interpretation_note(self):
        assert "confidence_interpretation_note" in ConfidenceReport.model_fields


# ---------------------------------------------------------------------------
# 4. TestProtenixBackend
# ---------------------------------------------------------------------------

class TestProtenixBackend:
    """Tests for Protenix backend enum and license mapping."""

    def test_protenix_in_backends(self):
        assert PredictionBackend.PROTENIX.value == "protenix"

    def test_protenix_commercial_ok(self):
        assert BACKEND_LICENSES[PredictionBackend.PROTENIX] == LicenseType.COMMERCIAL_OK

    def test_six_backends(self):
        assert len(PredictionBackend) == 6

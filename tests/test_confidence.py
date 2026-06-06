"""Tests for confidence interpretation logic."""


from foldcopilot.models.confidence import ConfidenceBucket, ConfidenceReport
from foldcopilot.tools.confidence import (
    _bucket,
    _check_hallucinations,
    _find_low_confidence_spans,
)


class TestBucket:
    def test_very_high(self):
        assert _bucket(95.0) == ConfidenceBucket.VERY_HIGH

    def test_high(self):
        assert _bucket(80.0) == ConfidenceBucket.HIGH

    def test_low(self):
        assert _bucket(60.0) == ConfidenceBucket.LOW

    def test_very_low(self):
        assert _bucket(30.0) == ConfidenceBucket.VERY_LOW

    def test_boundaries(self):
        assert _bucket(90.0) == ConfidenceBucket.HIGH  # <= 90 is HIGH
        assert _bucket(90.1) == ConfidenceBucket.VERY_HIGH
        assert _bucket(70.0) == ConfidenceBucket.LOW   # <= 70 is LOW
        assert _bucket(50.0) == ConfidenceBucket.VERY_LOW  # <= 50 is VERY_LOW


class TestLowConfidenceSpans:
    def test_no_low_regions(self):
        scores = [90.0] * 100
        spans = _find_low_confidence_spans(scores)
        assert spans == []

    def test_single_span(self):
        scores = [90.0] * 10 + [40.0] * 5 + [90.0] * 10
        spans = _find_low_confidence_spans(scores)
        assert len(spans) == 1
        assert spans[0].start == 11  # 1-indexed
        assert spans[0].end == 15
        assert spans[0].length == 5

    def test_multiple_spans(self):
        scores = [90.0] * 5 + [40.0] * 4 + [90.0] * 5 + [30.0] * 6 + [90.0] * 5
        spans = _find_low_confidence_spans(scores)
        assert len(spans) == 2

    def test_min_length_filter(self):
        scores = [90.0] * 5 + [40.0] * 2 + [90.0] * 5  # span too short
        spans = _find_low_confidence_spans(scores, min_length=3)
        assert len(spans) == 0

    def test_span_at_end(self):
        scores = [90.0] * 5 + [40.0] * 5
        spans = _find_low_confidence_spans(scores)
        assert len(spans) == 1
        assert spans[0].end == 10


class TestHallucinationDetection:
    def test_hallucination_high_severity(self):
        # AF is confident (pLDDT 85) in a region DisProt says is disordered
        scores = [85.0] * 100
        idr_regions = [{"start": 20, "end": 40, "annotation": "disorder"}]
        warnings = _check_hallucinations(scores, idr_regions, "disprot")
        assert len(warnings) == 1
        assert warnings[0].severity == "high"
        assert warnings[0].af_mean_plddt > 70

    def test_hallucination_moderate_severity(self):
        # AF is somewhat confident (pLDDT 60) in known IDR
        scores = [60.0] * 100
        idr_regions = [{"start": 10, "end": 30, "annotation": "disorder"}]
        warnings = _check_hallucinations(scores, idr_regions, "mobidb")
        assert len(warnings) == 1
        assert warnings[0].severity == "moderate"

    def test_no_hallucination_when_af_uncertain(self):
        # AF is not confident (pLDDT 30) — agrees with IDR
        scores = [30.0] * 100
        idr_regions = [{"start": 10, "end": 30, "annotation": "disorder"}]
        warnings = _check_hallucinations(scores, idr_regions, "disprot")
        assert len(warnings) == 0

    def test_empty_idr_regions(self):
        scores = [85.0] * 100
        warnings = _check_hallucinations(scores, [], "disprot")
        assert len(warnings) == 0


class TestConfidenceReport:
    def test_standard_caveats(self):
        report = ConfidenceReport(
            source="afdb",
            chain_summaries=[],
            overall_mean_plddt=80.0,
            overall_median_plddt=82.0,
        )
        report.add_standard_caveats()
        assert len(report.caveats) == 4
        assert any("hallucinate" in c for c in report.caveats)
        assert any("research use only" in c for c in report.caveats)

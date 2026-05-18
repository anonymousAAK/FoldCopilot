"""Tests for education mode — plain-language explanations."""

import pytest

from foldcopilot.tools.education import (
    explain_confidence_report,
    explain_hallucination_warning,
    explain_pae,
    explain_plddt,
)


class TestExplainPlddt:
    def test_very_high(self):
        result = explain_plddt(95.0)
        assert result["bucket"] == "very_high"
        assert "backbone AND side-chain" in result["what_it_means"]
        assert result["analogy"]
        assert result["citation"]

    def test_high(self):
        result = explain_plddt(80.0)
        assert result["bucket"] == "high"
        assert "backbone" in result["what_it_means"]

    def test_low(self):
        result = explain_plddt(60.0)
        assert result["bucket"] == "low"
        assert "DO NOT" in result["what_to_do"]

    def test_very_low(self):
        result = explain_plddt(30.0)
        assert result["bucket"] == "very_low"
        assert "IGNORE" in result["what_to_do"]
        assert "spaghetti" in result["analogy"]

    def test_boundary_90(self):
        result = explain_plddt(90.0)
        assert result["bucket"] == "high"  # <= 90 is HIGH

    def test_all_have_required_fields(self):
        for score in [95, 80, 60, 30]:
            result = explain_plddt(score)
            assert "summary" in result
            assert "what_it_means" in result
            assert "what_to_do" in result
            assert "analogy" in result
            assert "citation" in result


class TestExplainPae:
    def test_excellent(self):
        result = explain_pae(3.0)
        assert result["quality"] == "excellent"

    def test_moderate(self):
        result = explain_pae(7.0)
        assert result["quality"] == "moderate"

    def test_poor(self):
        result = explain_pae(15.0)
        assert result["quality"] == "poor"

    def test_very_poor(self):
        result = explain_pae(25.0)
        assert result["quality"] == "very poor"

    def test_interface_context(self):
        result = explain_pae(12.0, context="interface")
        assert result["context_note"] is not None
        assert "interface" in result["context_note"].lower()

    def test_domain_context(self):
        result = explain_pae(12.0, context="domain")
        assert result["context_note"] is not None
        assert "domain" in result["context_note"].lower()

    def test_general_no_context_note(self):
        result = explain_pae(12.0, context="general")
        assert result["context_note"] is None


class TestExplainHallucination:
    def test_high_severity(self):
        result = explain_hallucination_warning(85.0, "disprot", "high")
        assert "CONFIDENTLY WRONG" in result["headline"]
        assert result["severity"] == "high"
        assert "disprot" in result["explanation"]
        assert "background" in result

    def test_moderate_severity(self):
        result = explain_hallucination_warning(60.0, "mobidb", "moderate")
        assert "partially wrong" in result["headline"]
        assert "mobidb" in result["explanation"]

    def test_has_action_items(self):
        result = explain_hallucination_warning(85.0, "disprot", "high")
        assert "1." in result["what_to_do"]  # numbered list


class TestExplainConfidenceReport:
    def _make_report(self, plddt=85.0, hallucinations=None, idr_flags=None):
        return {
            "overall_mean_plddt": plddt,
            "hallucination_warnings": hallucinations or [],
            "idr_flags": idr_flags or [],
            "chain_summaries": [
                {
                    "low_confidence_spans": [
                        {"start": 50, "end": 60, "mean_plddt": 45.0, "length": 11}
                    ]
                }
            ],
            "pae_summary": {"mean_pae": 5.0},
        }

    def test_high_confidence_verdict(self):
        report = self._make_report(plddt=90.0)
        result = explain_confidence_report(report)
        assert "HIGH CONFIDENCE" in result["verdict"]

    def test_moderate_confidence(self):
        report = self._make_report(plddt=75.0)
        result = explain_confidence_report(report)
        assert "MODERATE" in result["verdict"]

    def test_low_confidence(self):
        report = self._make_report(plddt=55.0)
        result = explain_confidence_report(report)
        assert "LOW CONFIDENCE" in result["verdict"]

    def test_very_low_confidence(self):
        report = self._make_report(plddt=35.0)
        result = explain_confidence_report(report)
        assert "VERY LOW" in result["verdict"]

    def test_includes_plddt_explanation(self):
        report = self._make_report(plddt=85.0)
        result = explain_confidence_report(report)
        assert "plddt_explanation" in result
        assert result["plddt_explanation"]["bucket"] == "high"

    def test_includes_low_regions(self):
        report = self._make_report()
        result = explain_confidence_report(report)
        assert len(result["low_confidence_regions"]) > 0
        assert "50-60" in result["low_confidence_regions"][0]

    def test_hallucination_summary(self):
        report = self._make_report(
            plddt=75.0,
            hallucinations=[{
                "start": 20, "end": 40,
                "af_mean_plddt": 85.0,
                "idr_source": "disprot",
                "severity": "high",
            }],
        )
        result = explain_confidence_report(report)
        assert "DISORDERED" in result["hallucination_warnings"][0]

    def test_no_pae(self):
        report = self._make_report()
        report["pae_summary"] = None
        result = explain_confidence_report(report)
        assert result["pae_explanation"] is None

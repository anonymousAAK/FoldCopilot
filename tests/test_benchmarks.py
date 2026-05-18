"""Tests for benchmarking harness."""

from textwrap import dedent

import numpy as np
import pytest

from foldcopilot.tools.benchmarks import (
    _accuracy_by_plddt_bucket,
    evaluate_batch,
    evaluate_structure_pair,
    generate_benchmark_report,
    list_benchmark_datasets,
)


def _pdb_ca_line(serial, res_name, res_seq, x, y, z, b):
    return (
        f"ATOM  {serial:5d}  CA  {res_name:3s} A{res_seq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{b:6.2f}           C"
    )


def _make_pdb(residues):
    lines = []
    for i, (name, seq, x, y, z, b) in enumerate(residues, 1):
        lines.append(_pdb_ca_line(i, name, seq, x, y, z, b))
    lines.append("END")
    return "\n".join(lines)


class TestListBenchmarks:
    def test_has_datasets(self):
        result = list_benchmark_datasets()
        assert "datasets" in result
        assert "disprot_hallucination" in result["datasets"]
        assert "casp16_monomers" in result["datasets"]


class TestEvaluateStructurePair:
    def test_identical(self):
        pdb = _make_pdb([
            ("ALA", 1, 1.0, 0.0, 0.0, 90.0),
            ("GLY", 2, 4.0, 0.0, 0.0, 85.0),
            ("VAL", 3, 7.0, 0.0, 0.0, 80.0),
            ("LEU", 4, 10.0, 0.0, 0.0, 75.0),
            ("ILE", 5, 13.0, 0.0, 0.0, 92.0),
        ])
        result = evaluate_structure_pair(pdb, pdb, "test_identical")
        assert result["ca_rmsd"] == 0.0
        assert result["gdt_ts"] == 1.0
        assert result["gdt_1A"] == 1.0

    def test_shifted(self):
        pred = _make_pdb([
            ("ALA", 1, 0.0, 0.0, 0.0, 90.0),
            ("GLY", 2, 3.0, 0.0, 0.0, 85.0),
            ("VAL", 3, 6.0, 0.0, 0.0, 40.0),  # low pLDDT
        ])
        ref = _make_pdb([
            ("ALA", 1, 1.5, 0.0, 0.0, 90.0),
            ("GLY", 2, 4.5, 0.0, 0.0, 85.0),
            ("VAL", 3, 16.0, 0.0, 0.0, 40.0),  # 10A off
        ])
        result = evaluate_structure_pair(pred, ref, "test_shifted")
        assert result["ca_rmsd"] > 0
        assert result["gdt_ts"] < 1.0
        assert result["aligned_residues"] == 3

    def test_empty_pdb(self):
        result = evaluate_structure_pair("", "ATOM...", "bad")
        assert "error" in result

    def test_gdt_components(self):
        pdb = _make_pdb([
            ("ALA", 1, 0.0, 0.0, 0.0, 90.0),
            ("GLY", 2, 3.0, 0.0, 0.0, 80.0),
        ])
        result = evaluate_structure_pair(pdb, pdb)
        assert result["gdt_1A"] == 1.0
        assert result["gdt_2A"] == 1.0
        assert result["gdt_4A"] == 1.0
        assert result["gdt_8A"] == 1.0

    def test_plddt_accuracy_correlation(self):
        # High pLDDT residue is accurate, low pLDDT is inaccurate
        pred = _make_pdb([
            ("ALA", 1, 0.0, 0.0, 0.0, 95.0),  # high confidence
            ("GLY", 2, 3.0, 0.0, 0.0, 95.0),
            ("VAL", 3, 6.0, 0.0, 0.0, 30.0),  # low confidence
            ("LEU", 4, 9.0, 0.0, 0.0, 30.0),
        ])
        ref = _make_pdb([
            ("ALA", 1, 0.0, 0.0, 0.0, 0.0),   # matches
            ("GLY", 2, 3.0, 0.0, 0.0, 0.0),
            ("VAL", 3, 20.0, 0.0, 0.0, 0.0),  # far off
            ("LEU", 4, 25.0, 0.0, 0.0, 0.0),
        ])
        result = evaluate_structure_pair(pred, ref)
        # pLDDT should positively correlate with accuracy (neg correlate with distance)
        assert result["plddt_accuracy_correlation"] is not None


class TestAccuracyByBucket:
    def test_all_buckets(self):
        distances = {1: 0.5, 2: 1.0, 3: 3.0, 4: 8.0}
        plddt = {1: 95.0, 2: 80.0, 3: 60.0, 4: 30.0}
        result = _accuracy_by_plddt_bucket(distances, plddt)
        assert result["very_high_gt90"]["count"] == 1
        assert result["high_70_90"]["count"] == 1
        assert result["low_50_70"]["count"] == 1
        assert result["very_low_lt50"]["count"] == 1

    def test_empty(self):
        result = _accuracy_by_plddt_bucket({}, {})
        for bucket in result.values():
            assert bucket["count"] == 0


class TestEvaluateBatch:
    def test_batch(self):
        pdb = _make_pdb([
            ("ALA", 1, 1.0, 0.0, 0.0, 90.0),
            ("GLY", 2, 4.0, 0.0, 0.0, 85.0),
            ("VAL", 3, 7.0, 0.0, 0.0, 80.0),
        ])
        pairs = [
            {"predicted_pdb": pdb, "reference_pdb": pdb, "target_name": "t1"},
            {"predicted_pdb": pdb, "reference_pdb": pdb, "target_name": "t2"},
        ]
        result = evaluate_batch(pairs)
        assert result["summary"]["total_targets"] == 2
        assert result["summary"]["successful"] == 2
        assert result["summary"]["mean_rmsd"] == 0.0

    def test_empty_batch(self):
        result = evaluate_batch([])
        assert result["summary"]["total_targets"] == 0


class TestGenerateReport:
    def test_report_structure(self):
        pdb = _make_pdb([
            ("ALA", 1, 1.0, 0.0, 0.0, 90.0),
            ("GLY", 2, 4.0, 0.0, 0.0, 85.0),
            ("VAL", 3, 7.0, 0.0, 0.0, 80.0),
        ])
        batch = evaluate_batch([
            {"predicted_pdb": pdb, "reference_pdb": pdb, "target_name": "t1"},
        ])
        report = generate_benchmark_report(batch, "custom", "boltz2")
        assert report["report_type"] == "benchmark"
        assert report["leaderboard_entry"]["backend"] == "boltz2"
        assert report["leaderboard_entry"]["mean_rmsd"] == 0.0
        assert "citation" in report

"""Tests for notebook export."""

from foldcopilot.tools.notebook_export import (
    export_benchmark_notebook,
    export_confidence_notebook,
)


class TestConfidenceNotebook:
    def test_basic_structure(self):
        report = {"overall_mean_plddt": 85.0, "hallucination_warnings": []}
        result = export_confidence_notebook("P04637", report)
        assert result["format"] == "nbformat_v4"
        assert "P04637" in result["filename"]
        nb = result["notebook"]
        assert nb["nbformat"] == 4
        assert len(nb["cells"]) > 0

    def test_has_code_cells(self):
        report = {"overall_mean_plddt": 85.0}
        result = export_confidence_notebook("P00520", report)
        nb = result["notebook"]
        code_cells = [c for c in nb["cells"] if c["cell_type"] == "code"]
        assert len(code_cells) >= 3

    def test_has_markdown_cells(self):
        report = {}
        result = export_confidence_notebook("P00520", report)
        nb = result["notebook"]
        md_cells = [c for c in nb["cells"] if c["cell_type"] == "markdown"]
        assert len(md_cells) >= 2

    def test_metadata(self):
        report = {}
        result = export_confidence_notebook("P04637", report)
        meta = result["notebook"]["metadata"]
        assert meta["foldcopilot"]["uniprot_id"] == "P04637"
        assert meta["foldcopilot"]["version"] == "0.1.0"


class TestBenchmarkNotebook:
    def test_basic_structure(self):
        batch = {
            "summary": {"total_targets": 2, "mean_rmsd": 1.5},
            "per_target": [
                {"target": "t1", "ca_rmsd": 1.0, "gdt_ts": 0.9},
                {"target": "t2", "ca_rmsd": 2.0, "gdt_ts": 0.8},
            ],
        }
        result = export_benchmark_notebook(batch, "casp16", "boltz2")
        assert "casp16" in result["filename"]
        assert "boltz2" in result["filename"]
        assert result["notebook"]["nbformat"] == 4

    def test_has_cells(self):
        batch = {"summary": {}, "per_target": []}
        result = export_benchmark_notebook(batch)
        assert len(result["notebook"]["cells"]) >= 3

"""Tests for Protenix prediction client helpers."""

import json
from textwrap import dedent


from foldcopilot.clients.protenix_client import (
    _parse_output,
    _write_input_json,
    get_status,
)


class TestParseProtenixOutput:
    def test_empty_dir(self, tmp_path):
        result = _parse_output(tmp_path)
        assert "output_pdb_path" not in result

    def test_finds_pdb(self, tmp_path):
        pdb = dedent("""\
            ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 85.50           C
            END
        """)
        (tmp_path / "result.pdb").write_text(pdb)
        result = _parse_output(tmp_path)
        assert "output_pdb_path" in result
        assert result["output_pdb_path"].endswith(".pdb")

    def test_finds_cif(self, tmp_path):
        (tmp_path / "result.cif").write_text("data_test\n")
        result = _parse_output(tmp_path)
        assert "output_cif_path" in result
        assert result["output_cif_path"].endswith(".cif")

    def test_finds_confidence_json(self, tmp_path):
        conf = {"plddt": [85.0, 90.0, 70.0], "ptm": 0.85}
        (tmp_path / "confidence_scores.json").write_text(json.dumps(conf))
        result = _parse_output(tmp_path)
        assert result["mean_plddt"] is not None
        assert abs(result["mean_plddt"] - 81.7) < 0.1
        assert result["predicted_tm_score"] == 0.85

    def test_finds_iptm(self, tmp_path):
        conf = {"plddt": [80.0], "iptm": 0.72}
        (tmp_path / "confidence_summary.json").write_text(json.dumps(conf))
        result = _parse_output(tmp_path)
        assert result["predicted_iptm"] == 0.72


class TestInputWriters:
    def test_write_input_json(self, tmp_path):
        path = _write_input_json(["MKFL", "GHIJ"], tmp_path)
        assert path.exists()
        data = json.loads(path.read_text())
        assert data["name"] == "foldcopilot_job"
        assert len(data["sequences"]) == 2
        assert data["sequences"][0]["id"] == "A"
        assert data["sequences"][0]["sequence"] == "MKFL"
        assert data["sequences"][0]["type"] == "protein"
        assert data["sequences"][1]["id"] == "B"
        assert data["sequences"][1]["sequence"] == "GHIJ"
        assert data["modelSeeds"] == [42]

    def test_write_single_chain(self, tmp_path):
        path = _write_input_json(["ACDEFG"], tmp_path)
        data = json.loads(path.read_text())
        assert len(data["sequences"]) == 1
        assert data["sequences"][0]["id"] == "A"
        assert data["sequences"][0]["sequence"] == "ACDEFG"


class TestGetStatus:
    def test_returns_expected_keys(self):
        status = get_status()
        assert isinstance(status, dict)
        assert "installed" in status
        assert "path" in status
        assert "gpu_available" in status
        assert "prediction_dir" in status
        assert "setup_instructions" in status
        assert isinstance(status["installed"], bool)
        assert isinstance(status["gpu_available"], bool)

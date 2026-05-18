"""Tests for prediction models, license routing, and Boltz-2 client helpers."""

import json
from pathlib import Path
from textwrap import dedent

import pytest

from foldcopilot.clients.boltz2_client import (
    _parse_boltz_output,
    _parse_plddt_from_pdb,
    _write_input_fasta,
    _write_input_yaml,
)
from foldcopilot.models.prediction import (
    BACKEND_LICENSES,
    LicenseType,
    PredictionBackend,
    PredictionInput,
    ReproducibilityManifest,
)
from foldcopilot.tools.predict import check_license_compatibility, list_backends


class TestLicenseRouting:
    def test_boltz2_commercial_ok(self):
        err = check_license_compatibility(PredictionBackend.BOLTZ2, commercial=True)
        assert err is None

    def test_boltz2_noncommercial_ok(self):
        err = check_license_compatibility(PredictionBackend.BOLTZ2, commercial=False)
        assert err is None

    def test_all_current_backends_commercial_ok(self):
        # All v0.3 backends are MIT/Apache
        for backend in PredictionBackend:
            assert BACKEND_LICENSES[backend] == LicenseType.COMMERCIAL_OK

    def test_list_backends(self):
        result = list_backends()
        assert "backends" in result
        names = [b["name"] for b in result["backends"]]
        assert "boltz2" in names
        assert "openfold3" in names
        assert "chai1" in names
        # All backends now implemented
        for b in result["backends"]:
            assert b["implemented"] is True


class TestReproducibilityManifest:
    def test_create(self):
        manifest = ReproducibilityManifest.create(
            backend="boltz2",
            sequences=["MKFL", "GHIJ"],
            parameters={"recycling_steps": 3, "sampling_steps": 200},
            backend_version="2.0.0",
            seed=42,
        )
        assert manifest.backend == "boltz2"
        assert manifest.backend_version == "2.0.0"
        assert manifest.seed == 42
        assert manifest.input_sequence_hash  # not empty
        assert manifest.runtime_env["python_version"]
        assert manifest.timestamp_utc > 0

    def test_deterministic_hash(self):
        m1 = ReproducibilityManifest.create(
            backend="boltz2", sequences=["MKFL", "GHIJ"], parameters={}
        )
        m2 = ReproducibilityManifest.create(
            backend="boltz2", sequences=["MKFL", "GHIJ"], parameters={}
        )
        assert m1.input_sequence_hash == m2.input_sequence_hash

    def test_different_sequences_different_hash(self):
        m1 = ReproducibilityManifest.create(
            backend="boltz2", sequences=["MKFL"], parameters={}
        )
        m2 = ReproducibilityManifest.create(
            backend="boltz2", sequences=["GHIJ"], parameters={}
        )
        assert m1.input_sequence_hash != m2.input_sequence_hash

    def test_serializable(self):
        manifest = ReproducibilityManifest.create(
            backend="boltz2", sequences=["MKFL"], parameters={"a": 1}
        )
        data = json.loads(manifest.model_dump_json())
        assert data["backend"] == "boltz2"


class TestPredictionInput:
    def test_defaults(self):
        inp = PredictionInput(sequences=["MKFL"])
        assert inp.backend == PredictionBackend.BOLTZ2
        assert inp.recycling_steps == 3
        assert inp.sampling_steps == 200
        assert inp.use_msa is True
        assert inp.predict_affinity is False

    def test_requires_sequence(self):
        with pytest.raises(Exception):
            PredictionInput(sequences=[])


class TestInputWriters:
    def test_write_fasta(self, tmp_path):
        path = _write_input_fasta(["MKFL", "GHIJ"], tmp_path)
        content = path.read_text()
        assert ">chain_A" in content
        assert "MKFL" in content
        assert ">chain_B" in content
        assert "GHIJ" in content

    def test_write_yaml(self, tmp_path):
        path = _write_input_yaml(["MKFL", "GHIJ"], tmp_path)
        content = path.read_text()
        assert "version: 1" in content
        assert "sequence: MKFL" in content
        assert "sequence: GHIJ" in content

    def test_write_yaml_with_affinity(self, tmp_path):
        path = _write_input_yaml(["MKFL", "GHIJ"], tmp_path, predict_affinity=True)
        content = path.read_text()
        assert "affinity" in content
        assert "binder: A" in content
        assert "target: B" in content


class TestParsePlddt:
    def test_parse_from_pdb(self, tmp_path):
        # Minimal PDB with CA atoms and B-factors (pLDDT)
        pdb = dedent("""\
            ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 85.50           C
            ATOM      2  CB  ALA A   1       1.500   2.500   3.500  1.00 80.00           C
            ATOM      3  CA  GLY A   2       4.000   5.000   6.000  1.00 92.30           C
            ATOM      4  CA  VAL A   3       7.000   8.000   9.000  1.00 45.10           C
            END
        """)
        pdb_path = tmp_path / "test.pdb"
        pdb_path.write_text(pdb)
        mean = _parse_plddt_from_pdb(pdb_path)
        assert mean is not None
        # Mean of CA B-factors: (85.5 + 92.3 + 45.1) / 3 = 74.3
        assert abs(mean - 74.3) < 0.1


class TestParseBoltzOutput:
    def test_empty_dir(self, tmp_path):
        result = _parse_boltz_output(tmp_path)
        assert "output_pdb_path" not in result

    def test_finds_pdb(self, tmp_path):
        pdb = dedent("""\
            ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 85.50           C
            END
        """)
        (tmp_path / "result.pdb").write_text(pdb)
        result = _parse_boltz_output(tmp_path)
        assert "output_pdb_path" in result
        assert result["output_pdb_path"].endswith(".pdb")

    def test_finds_confidence_json(self, tmp_path):
        conf = {"plddt": [85.0, 90.0, 70.0], "ptm": 0.85}
        (tmp_path / "confidence_scores.json").write_text(json.dumps(conf))
        result = _parse_boltz_output(tmp_path)
        assert result["mean_plddt"] is not None
        assert abs(result["mean_plddt"] - 81.7) < 0.1
        assert result["predicted_tm_score"] == 0.85

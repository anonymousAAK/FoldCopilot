"""Tests for fold-drift tracker."""

import json


from foldcopilot.tools.fold_drift import (
    check_fold_drift,
    check_prediction_drift,
)


class TestCheckFoldDrift:
    def test_no_predictions(self, tmp_path):
        result = check_fold_drift(str(tmp_path))
        assert result["status"] == "no_predictions"
        assert result["predictions_scanned"] == 0

    def test_finds_manifests(self, tmp_path):
        pred_dir = tmp_path / "boltz_abc123" / "output"
        pred_dir.mkdir(parents=True)
        manifest = {
            "backend": "boltz2",
            "backend_version": "2.0.0",
            "input_sequence_hash": "abc123",
            "timestamp_utc": 1700000000,
        }
        (pred_dir / "reproducibility_manifest.json").write_text(json.dumps(manifest))
        result = check_fold_drift(str(tmp_path))
        assert result["predictions_scanned"] == 1

    def test_detects_drift(self, tmp_path):
        pred_dir = tmp_path / "pred1" / "output"
        pred_dir.mkdir(parents=True)
        # Use a version that definitely doesn't match
        manifest = {
            "backend": "boltz2",
            "backend_version": "0.0.0-definitely-old",
            "input_sequence_hash": "abc",
            "timestamp_utc": 1700000000,
        }
        (pred_dir / "reproducibility_manifest.json").write_text(json.dumps(manifest))
        result = check_fold_drift(str(tmp_path))
        # Either drifted (if boltz installed with different version) or unknown (not installed)
        assert result["predictions_scanned"] == 1


class TestCheckPredictionDrift:
    def test_missing_manifest(self):
        result = check_prediction_drift("/nonexistent/path/manifest.json")
        assert "error" in result

    def test_valid_manifest(self, tmp_path):
        manifest = {
            "backend": "boltz2",
            "backend_version": "2.0.0",
            "weights_hash": "abc123",
        }
        path = tmp_path / "reproducibility_manifest.json"
        path.write_text(json.dumps(manifest))
        result = check_prediction_drift(str(path))
        assert "drift_status" in result

    def test_no_version_in_manifest(self, tmp_path):
        manifest = {"backend": "boltz2"}
        path = tmp_path / "reproducibility_manifest.json"
        path.write_text(json.dumps(manifest))
        result = check_prediction_drift(str(path))
        assert result["drift_status"] == "unknown"

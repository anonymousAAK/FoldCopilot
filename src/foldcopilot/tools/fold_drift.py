"""Fold-drift tracker — detect when backend updates would change stored predictions.

Model-agnostic: tracks backend version + weights hash from reproducibility manifests.
Alerts users when a stored prediction was made with a different backend version than
currently installed.
"""

from __future__ import annotations

import json
from pathlib import Path

PREDICTION_DIR = Path.home() / ".cache" / "foldcopilot" / "predictions"


def _scan_manifests(prediction_dir: Path | None = None) -> list[dict]:
    """Scan prediction directory for reproducibility manifests."""
    scan_dir = prediction_dir or PREDICTION_DIR
    if not scan_dir.exists():
        return []

    manifests = []
    for manifest_path in scan_dir.rglob("reproducibility_manifest.json"):
        try:
            data = json.loads(manifest_path.read_text())
            data["_manifest_path"] = str(manifest_path)
            data["_prediction_dir"] = str(manifest_path.parent)
            manifests.append(data)
        except (json.JSONDecodeError, OSError):
            continue

    return manifests


def _get_current_backend_versions() -> dict[str, str | None]:
    """Get currently installed backend versions."""
    versions: dict[str, str | None] = {}

    # Boltz-2
    try:
        import importlib.metadata
        versions["boltz2"] = importlib.metadata.version("boltz")
    except Exception:
        versions["boltz2"] = None

    # OpenFold3
    try:
        import importlib.metadata
        versions["openfold3"] = importlib.metadata.version("openfold3")
    except Exception:
        versions["openfold3"] = None

    # Chai-1
    try:
        import importlib.metadata
        versions["chai1"] = importlib.metadata.version("chai-lab")
    except Exception:
        versions["chai1"] = None

    return versions


def check_fold_drift(prediction_dir: str | None = None) -> dict:
    """Check all stored predictions for fold drift.

    Compares the backend version in each stored reproducibility manifest
    against the currently installed version. Flags predictions that may
    produce different results if re-run.
    """
    scan_path = Path(prediction_dir) if prediction_dir else None
    manifests = _scan_manifests(scan_path)

    if not manifests:
        return {
            "status": "no_predictions",
            "message": "No stored predictions found.",
            "predictions_scanned": 0,
        }

    current_versions = _get_current_backend_versions()
    drifted = []
    stable = []
    unknown = []

    for m in manifests:
        backend = m.get("backend", "unknown")
        stored_version = m.get("backend_version")
        current_version = current_versions.get(backend)

        entry = {
            "backend": backend,
            "stored_version": stored_version,
            "current_version": current_version,
            "input_sequence_hash": m.get("input_sequence_hash", ""),
            "timestamp_utc": m.get("timestamp_utc", 0),
            "prediction_dir": m.get("_prediction_dir", ""),
        }

        if stored_version is None or current_version is None:
            entry["drift_status"] = "unknown"
            unknown.append(entry)
        elif stored_version != current_version:
            entry["drift_status"] = "drifted"
            entry["message"] = (
                f"{backend} updated: {stored_version} → {current_version}. "
                f"Re-running may produce different results."
            )
            drifted.append(entry)
        else:
            entry["drift_status"] = "stable"
            stable.append(entry)

    return {
        "status": "drift_detected" if drifted else "stable",
        "predictions_scanned": len(manifests),
        "drifted": drifted,
        "stable_count": len(stable),
        "unknown_count": len(unknown),
        "current_backend_versions": current_versions,
        "recommendation": (
            f"{len(drifted)} prediction(s) may produce different results with "
            f"current backend versions. Consider re-running for reproducibility."
            if drifted else "All predictions are consistent with installed backends."
        ),
    }


def check_prediction_drift(manifest_path: str) -> dict:
    """Check a single prediction's reproducibility manifest for drift."""
    path = Path(manifest_path)
    if not path.exists():
        return {"error": f"Manifest not found: {manifest_path}"}

    try:
        manifest = json.loads(path.read_text())
    except json.JSONDecodeError:
        return {"error": f"Invalid JSON in manifest: {manifest_path}"}

    backend = manifest.get("backend", "unknown")
    stored_version = manifest.get("backend_version")
    current_versions = _get_current_backend_versions()
    current_version = current_versions.get(backend)

    if stored_version is None:
        return {"drift_status": "unknown", "reason": "No version recorded in manifest."}

    if current_version is None:
        return {
            "drift_status": "unknown",
            "reason": f"{backend} not currently installed.",
        }

    if stored_version != current_version:
        return {
            "drift_status": "drifted",
            "backend": backend,
            "stored_version": stored_version,
            "current_version": current_version,
            "weights_hash": manifest.get("weights_hash"),
            "recommendation": "Re-run prediction with current version for reproducibility.",
        }

    return {
        "drift_status": "stable",
        "backend": backend,
        "version": stored_version,
    }

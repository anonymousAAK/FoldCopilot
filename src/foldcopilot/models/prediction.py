"""Data models for structure prediction jobs and reproducibility manifests."""

from __future__ import annotations

import hashlib
import platform
import time
from enum import Enum
from typing import Any

from pydantic import BaseModel, Field


class PredictionBackend(str, Enum):
    BOLTZ2 = "boltz2"
    OPENFOLD3 = "openfold3"
    CHAI1 = "chai1"
    PROTENIX = "protenix"
    AF3 = "alphafold3"       # BYO-weights only, non-commercial
    AQAFFINITY = "aqaffinity"  # SandboxAQ, on top of OpenFold3


class LicenseType(str, Enum):
    COMMERCIAL_OK = "commercial_ok"
    NON_COMMERCIAL = "non_commercial"


# Backend license routing table
BACKEND_LICENSES: dict[PredictionBackend, LicenseType] = {
    PredictionBackend.BOLTZ2: LicenseType.COMMERCIAL_OK,       # MIT
    PredictionBackend.OPENFOLD3: LicenseType.COMMERCIAL_OK,    # Apache-2.0
    PredictionBackend.CHAI1: LicenseType.COMMERCIAL_OK,        # Apache-2.0
    PredictionBackend.PROTENIX: LicenseType.COMMERCIAL_OK,     # Apache-2.0
    PredictionBackend.AF3: LicenseType.NON_COMMERCIAL,         # CC-BY-NC-SA 4.0 + non-commercial weights
    PredictionBackend.AQAFFINITY: LicenseType.COMMERCIAL_OK,   # Open
}


class JobStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETE = "complete"
    FAILED = "failed"


class PredictionInput(BaseModel):
    """Input for a structure prediction job."""

    sequences: list[str] = Field(..., min_length=1)
    backend: PredictionBackend = PredictionBackend.BOLTZ2
    commercial_use: bool = False
    af3_noncommercial_attestation: bool = False  # Required for AF3
    # Boltz-2 specific
    recycling_steps: int = 3
    sampling_steps: int = 200
    diffusion_samples: int = 1
    use_msa: bool = True
    predict_affinity: bool = False


class ReproducibilityManifest(BaseModel):
    """Ships with every prediction output. Ensures reproducibility."""

    backend: str
    backend_version: str | None = None
    weights_hash: str | None = None
    training_data_source: str | None = None
    input_sequence_hash: str
    parameters: dict[str, Any]
    seed: int | None = None
    runtime_env: dict[str, str]
    gpu_type: str | None = None
    timestamp_utc: float
    foldcopilot_version: str = "1.1.0"

    @classmethod
    def create(
        cls,
        backend: str,
        sequences: list[str],
        parameters: dict[str, Any],
        backend_version: str | None = None,
        weights_hash: str | None = None,
        training_data_source: str | None = None,
        seed: int | None = None,
        gpu_type: str | None = None,
    ) -> ReproducibilityManifest:
        seq_hash = hashlib.sha256(
            "|".join(sorted(sequences)).encode()
        ).hexdigest()

        return cls(
            backend=backend,
            backend_version=backend_version,
            weights_hash=weights_hash,
            training_data_source=training_data_source,
            input_sequence_hash=seq_hash,
            parameters=parameters,
            seed=seed,
            runtime_env={
                "python_version": platform.python_version(),
                "platform": platform.platform(),
                "architecture": platform.machine(),
            },
            gpu_type=gpu_type,
            timestamp_utc=time.time(),
        )


class PredictionResult(BaseModel):
    """Result of a structure prediction job."""

    job_id: str
    status: JobStatus
    backend: str
    sequences: list[str]

    # Output paths/URIs (not raw content — keep MCP responses compact)
    output_pdb_path: str | None = None
    output_cif_path: str | None = None
    confidence_json_path: str | None = None

    # Compact summary (always included in MCP response)
    mean_plddt: float | None = None
    per_chain_plddt: dict[str, float] | None = None
    predicted_tm_score: float | None = None
    predicted_affinity: float | None = None  # Boltz-2 affinity

    manifest: ReproducibilityManifest | None = None
    error_message: str | None = None
    elapsed_seconds: float | None = None

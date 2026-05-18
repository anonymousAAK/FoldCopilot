"""Prediction tools — run structure predictions via supported backends.

v0.3: Boltz-2 (MIT, fast, includes affinity).
v0.5: + OpenFold3 (Apache-2.0, AF3 reproduction) + Chai-1 (Apache-2.0).
"""

from __future__ import annotations

from foldcopilot.clients import boltz2_client, chai1_client, openfold3_client
from foldcopilot.models.prediction import (
    BACKEND_LICENSES,
    LicenseType,
    PredictionBackend,
    PredictionInput,
)


def check_license_compatibility(
    backend: PredictionBackend, commercial: bool
) -> str | None:
    """Check if backend is compatible with commercial use. Returns error message or None."""
    if commercial and BACKEND_LICENSES.get(backend) == LicenseType.NON_COMMERCIAL:
        return (
            f"{backend.value} is not licensed for commercial use. "
            f"Use boltz2 (MIT) or openfold3 (Apache-2.0) instead."
        )
    return None


async def predict_structure(
    sequences: list[str],
    backend: str = "boltz2",
    commercial_use: bool = False,
    recycling_steps: int = 3,
    sampling_steps: int = 200,
    diffusion_samples: int = 1,
    use_msa: bool = True,
    predict_affinity: bool = False,
) -> dict:
    """Run a structure prediction using the specified backend.

    Returns a PredictionResult with output paths, confidence scores,
    and a reproducibility manifest.
    """
    # Validate backend
    try:
        backend_enum = PredictionBackend(backend)
    except ValueError:
        available = [b.value for b in PredictionBackend]
        return {
            "error": f"Unknown backend '{backend}'. Available: {available}",
            "status": "failed",
        }

    # License check
    license_err = check_license_compatibility(backend_enum, commercial_use)
    if license_err:
        return {"error": license_err, "status": "failed"}

    prediction_input = PredictionInput(
        sequences=sequences,
        backend=backend_enum,
        commercial_use=commercial_use,
        recycling_steps=recycling_steps,
        sampling_steps=sampling_steps,
        diffusion_samples=diffusion_samples,
        use_msa=use_msa,
        predict_affinity=predict_affinity,
    )

    # Route to the right backend
    if backend_enum == PredictionBackend.BOLTZ2:
        result = await boltz2_client.predict_local(prediction_input)
    elif backend_enum == PredictionBackend.OPENFOLD3:
        result = await openfold3_client.predict_local(prediction_input)
    elif backend_enum == PredictionBackend.CHAI1:
        result = await chai1_client.predict_local(prediction_input)
    else:
        return {"error": f"Backend {backend} not yet implemented.", "status": "failed"}

    return result.model_dump()


def get_backend_status(backend: str = "boltz2") -> dict:
    """Check installation status of a prediction backend."""
    status_map = {
        "boltz2": boltz2_client.get_boltz_status,
        "openfold3": openfold3_client.get_status,
        "chai1": chai1_client.get_status,
    }
    fn = status_map.get(backend)
    if fn is None:
        return {"error": f"Unknown backend: {backend}. Available: {list(status_map)}"}
    return fn()


def list_backends() -> dict:
    """List all available prediction backends and their status."""
    status_fns = {
        PredictionBackend.BOLTZ2: boltz2_client.get_boltz_status,
        PredictionBackend.OPENFOLD3: openfold3_client.get_status,
        PredictionBackend.CHAI1: chai1_client.get_status,
    }

    backends = []
    for b in PredictionBackend:
        status = status_fns[b]()
        backends.append({
            "name": b.value,
            "license": BACKEND_LICENSES[b].value,
            "implemented": True,
            "installed": status.get("installed", False),
            "gpu_available": status.get("gpu_available", False),
        })

    return {"backends": backends}

"""Prediction tools — run structure predictions via supported backends.

v0.3: Boltz-2 (MIT, fast, includes affinity).
v0.5: + OpenFold3 (Apache-2.0, AF3 reproduction) + Chai-1 (Apache-2.0).
v0.9: + AF3 (BYO-weights, non-commercial) + AQAffinity (on top of OpenFold3).
"""

from __future__ import annotations

from foldcopilot.clients import boltz2_client, chai1_client, openfold3_client
from foldcopilot.models.prediction import (
    BACKEND_LICENSES,
    LicenseType,
    PredictionBackend,
    PredictionInput,
)
from foldcopilot.utils.validation import ValidationError, validate_sequences


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
    af3_noncommercial_attestation: bool = False,
    recycling_steps: int = 3,
    sampling_steps: int = 200,
    diffusion_samples: int = 1,
    use_msa: bool = True,
    predict_affinity: bool = False,
    ctx: object | None = None,
) -> dict:
    """Run a structure prediction using the specified backend.

    Returns a PredictionResult with output paths, confidence scores,
    and a reproducibility manifest.
    """
    # Validate sequences at system boundary
    try:
        sequences = validate_sequences(sequences)
    except ValidationError as e:
        return {"error": str(e), "status": "failed"}

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

    # AF3 requires explicit non-commercial attestation (BYO-weights gate)
    if backend_enum == PredictionBackend.AF3 and not af3_noncommercial_attestation:
        return {
            "error": (
                "AF3 requires af3_noncommercial_attestation=True. "
                "AlphaFold 3 weights are licensed under CC-BY-NC-SA 4.0 and "
                "may only be used for non-commercial purposes. You must supply "
                "your own weights (BYO-weights). Set af3_noncommercial_attestation=True "
                "to confirm non-commercial use."
            ),
            "status": "failed",
        }

    prediction_input = PredictionInput(
        sequences=sequences,
        backend=backend_enum,
        commercial_use=commercial_use,
        af3_noncommercial_attestation=af3_noncommercial_attestation,
        recycling_steps=recycling_steps,
        sampling_steps=sampling_steps,
        diffusion_samples=diffusion_samples,
        use_msa=use_msa,
        predict_affinity=predict_affinity,
    )

    # Report progress: MSA + inference stages
    async def _progress(step: int, total: int, message: str) -> None:
        if ctx and hasattr(ctx, "report_progress"):
            await ctx.report_progress(step, total, message=message)

    await _progress(1, 4, f"Starting {backend} prediction")
    await _progress(2, 4, "Running MSA alignment" if use_msa else "Skipping MSA (single-sequence mode)")

    # Route to the right backend
    if backend_enum == PredictionBackend.BOLTZ2:
        result = await boltz2_client.predict_local(prediction_input)
    elif backend_enum == PredictionBackend.OPENFOLD3:
        result = await openfold3_client.predict_local(prediction_input)
    elif backend_enum == PredictionBackend.CHAI1:
        result = await chai1_client.predict_local(prediction_input)
    elif backend_enum == PredictionBackend.AF3:
        result = await openfold3_client.predict_local(prediction_input, af3_mode=True)
    elif backend_enum == PredictionBackend.AQAFFINITY:
        result = await openfold3_client.predict_local(prediction_input, aqaffinity_mode=True)
    else:
        return {"error": f"Backend {backend} not yet implemented.", "status": "failed"}

    await _progress(3, 4, "Parsing output and building reproducibility manifest")

    return result.model_dump()


def _get_af3_status() -> dict:
    """Check AF3 BYO-weights status."""
    base = openfold3_client.get_status()
    return {
        "installed": base.get("installed", False),
        "gpu_available": base.get("gpu_available", False),
        "note": "AF3 requires BYO-weights (CC-BY-NC-SA 4.0, non-commercial only)",
        "weights_required": True,
    }


def _get_aqaffinity_status() -> dict:
    """Check AQAffinity status (runs on top of OpenFold3)."""
    base = openfold3_client.get_status()
    return {
        "installed": base.get("installed", False),
        "gpu_available": base.get("gpu_available", False),
        "note": "AQAffinity runs on top of OpenFold3 with SandboxAQ affinity head",
    }


def get_backend_status(backend: str = "boltz2") -> dict:
    """Check installation status of a prediction backend."""
    status_map = {
        "boltz2": boltz2_client.get_boltz_status,
        "openfold3": openfold3_client.get_status,
        "chai1": chai1_client.get_status,
        "alphafold3": _get_af3_status,
        "aqaffinity": _get_aqaffinity_status,
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
        PredictionBackend.AF3: _get_af3_status,
        PredictionBackend.AQAFFINITY: _get_aqaffinity_status,
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

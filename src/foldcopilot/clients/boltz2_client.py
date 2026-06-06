"""Boltz-2 prediction client.

Wraps Boltz-2 (MIT license, Passaro et al. 2025) for local or remote
structure prediction. Client-only — uses user-supplied compute.

Supports:
- Local execution via `boltz predict` CLI
- Remote execution via Modal (reference deploy)
"""

from __future__ import annotations

import asyncio
import json
import os
import shutil
import time
import uuid
from pathlib import Path
from typing import Any

import httpx

from foldcopilot.models.prediction import (
    JobStatus,
    PredictionInput,
    PredictionResult,
    ReproducibilityManifest,
)

# Output directory for predictions
PREDICTION_DIR = Path.home() / ".cache" / "foldcopilot" / "predictions"


def _check_boltz_available() -> str | None:
    """Check if boltz CLI is available. Returns path or None."""
    return shutil.which("boltz")


def _write_input_fasta(sequences: list[str], work_dir: Path) -> Path:
    """Write sequences to a FASTA file for Boltz-2 input."""
    fasta_path = work_dir / "input.fasta"
    with open(fasta_path, "w") as f:
        for i, seq in enumerate(sequences):
            chain_id = chr(65 + i)  # A, B, C, ...
            f.write(f">chain_{chain_id}\n{seq}\n")
    return fasta_path


def _write_input_yaml(
    sequences: list[str], work_dir: Path, predict_affinity: bool = False
) -> Path:
    """Write Boltz-2 YAML input format for multi-chain or affinity prediction."""
    yaml_path = work_dir / "input.yaml"
    lines = ["version: 1", "sequences:"]
    for i, seq in enumerate(sequences):
        chain_id = chr(65 + i)
        lines.append("  - protein:")
        lines.append(f"      id: {chain_id}")
        lines.append(f"      sequence: {seq}")

    if predict_affinity and len(sequences) >= 2:
        lines.append("properties:")
        lines.append("  - affinity:")
        lines.append("      binder: A")
        lines.append("      target: B")

    yaml_path.write_text("\n".join(lines))
    return yaml_path


def _parse_boltz_output(output_dir: Path) -> dict[str, Any]:
    """Parse Boltz-2 output directory for results."""
    result: dict[str, Any] = {}

    # Find output structure files
    for ext in ("pdb", "cif"):
        candidates = list(output_dir.rglob(f"*.{ext}"))
        if candidates:
            result[f"output_{ext}_path"] = str(candidates[0])

    # Find confidence scores
    confidence_files = list(output_dir.rglob("*confidence*")) + list(
        output_dir.rglob("*scores*")
    )
    if confidence_files:
        result["confidence_json_path"] = str(confidence_files[0])

    # Parse confidence JSON if available
    for cf in confidence_files:
        if cf.suffix == ".json":
            try:
                data = json.loads(cf.read_text())
                if "plddt" in data:
                    import numpy as np
                    scores = data["plddt"]
                    if isinstance(scores, list):
                        result["mean_plddt"] = round(float(np.mean(scores)), 1)
                if "ptm" in data:
                    result["predicted_tm_score"] = data["ptm"]
                if "affinity" in data:
                    result["predicted_affinity"] = data["affinity"]
            except (json.JSONDecodeError, ImportError):
                pass

    # Try parsing pLDDT from PDB B-factors if no confidence JSON
    if "mean_plddt" not in result and "output_pdb_path" in result:
        try:
            result["mean_plddt"] = _parse_plddt_from_pdb(
                Path(result["output_pdb_path"])
            )
        except Exception:
            pass

    return result


def _parse_plddt_from_pdb(pdb_path: Path) -> float | None:
    """Extract mean pLDDT from PDB B-factor column (CA atoms)."""
    import numpy as np

    scores = []
    for line in pdb_path.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                scores.append(float(line[60:66].strip()))
            except ValueError:
                pass
    if scores:
        return round(float(np.mean(scores)), 1)
    return None


async def predict_local(
    prediction_input: PredictionInput,
) -> PredictionResult:
    """Run Boltz-2 prediction locally via CLI.

    Requires `boltz` to be installed and accessible in PATH.
    """
    job_id = str(uuid.uuid4())[:8]
    start_time = time.time()

    boltz_path = _check_boltz_available()
    if not boltz_path:
        return PredictionResult(
            job_id=job_id,
            status=JobStatus.FAILED,
            backend="boltz2",
            sequences=prediction_input.sequences,
            error_message=(
                "Boltz-2 CLI not found in PATH. Install with: "
                "pip install boltz (requires GPU). "
                "See https://github.com/jwohlwend/boltz for setup."
            ),
        )

    # Create working directory
    work_dir = PREDICTION_DIR / job_id
    work_dir.mkdir(parents=True, exist_ok=True)
    output_dir = work_dir / "output"
    output_dir.mkdir(exist_ok=True)

    try:
        # Write input
        if prediction_input.predict_affinity or len(prediction_input.sequences) > 1:
            input_path = _write_input_yaml(
                prediction_input.sequences,
                work_dir,
                predict_affinity=prediction_input.predict_affinity,
            )
        else:
            input_path = _write_input_fasta(prediction_input.sequences, work_dir)

        # Build command
        cmd = [
            boltz_path,
            "predict",
            str(input_path),
            "--out_dir", str(output_dir),
            "--recycling_steps", str(prediction_input.recycling_steps),
            "--sampling_steps", str(prediction_input.sampling_steps),
            "--diffusion_samples", str(prediction_input.diffusion_samples),
        ]

        if not prediction_input.use_msa:
            cmd.append("--no_msa")

        # Run
        process = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            cwd=str(work_dir),
        )
        stdout, stderr = await process.communicate()

        elapsed = time.time() - start_time

        if process.returncode != 0:
            return PredictionResult(
                job_id=job_id,
                status=JobStatus.FAILED,
                backend="boltz2",
                sequences=prediction_input.sequences,
                error_message=f"Boltz-2 exited with code {process.returncode}: {stderr.decode()[:500]}",
                elapsed_seconds=round(elapsed, 1),
            )

        # Parse output
        parsed = _parse_boltz_output(output_dir)

        # Build reproducibility manifest
        manifest = ReproducibilityManifest.create(
            backend="boltz2",
            sequences=prediction_input.sequences,
            parameters={
                "recycling_steps": prediction_input.recycling_steps,
                "sampling_steps": prediction_input.sampling_steps,
                "diffusion_samples": prediction_input.diffusion_samples,
                "use_msa": prediction_input.use_msa,
                "predict_affinity": prediction_input.predict_affinity,
            },
            training_data_source="Boltz-2 training data (Passaro et al., bioRxiv 2025.06.14.659707)",
        )

        # Save manifest
        manifest_path = output_dir / "reproducibility_manifest.json"
        manifest_path.write_text(manifest.model_dump_json(indent=2))

        return PredictionResult(
            job_id=job_id,
            status=JobStatus.COMPLETE,
            backend="boltz2",
            sequences=prediction_input.sequences,
            output_pdb_path=parsed.get("output_pdb_path"),
            output_cif_path=parsed.get("output_cif_path"),
            confidence_json_path=parsed.get("confidence_json_path"),
            mean_plddt=parsed.get("mean_plddt"),
            predicted_tm_score=parsed.get("predicted_tm_score"),
            predicted_affinity=parsed.get("predicted_affinity"),
            manifest=manifest,
            elapsed_seconds=round(elapsed, 1),
        )

    except Exception as e:
        return PredictionResult(
            job_id=job_id,
            status=JobStatus.FAILED,
            backend="boltz2",
            sequences=prediction_input.sequences,
            error_message=str(e),
            elapsed_seconds=round(time.time() - start_time, 1),
        )


def get_boltz_status() -> dict:
    """Check Boltz-2 installation status and environment."""
    boltz_path = _check_boltz_available()
    gpu_available = _check_gpu()
    return {
        "installed": boltz_path is not None,
        "path": boltz_path,
        "gpu_available": gpu_available,
        "prediction_dir": str(PREDICTION_DIR),
        "setup_instructions": (
            "Install Boltz-2: pip install boltz\n"
            "Requires NVIDIA GPU with CUDA. See https://github.com/jwohlwend/boltz\n"
            "For cloud GPU: use Modal (modal.com) or RunPod (runpod.io)"
        ) if not boltz_path else None,
    }


def _check_gpu() -> bool:
    """Check if CUDA GPU is available."""
    try:
        import torch
        return torch.cuda.is_available()
    except ImportError:
        return False


async def predict_nim(
    prediction_input: PredictionInput,
    nim_endpoint: str = "https://health.api.nvidia.com/v1/biology/mit/boltz-2",
    api_key: str | None = None,
) -> PredictionResult:
    """Run Boltz-2 prediction via NVIDIA NIM REST API.

    Alternative to local CLI — enables GPU-less cloud inference.
    Requires NVIDIA API key (env: NVIDIA_API_KEY or passed directly).
    NIM v1.5+ returns PAE matrix in response.
    """
    key = api_key or os.environ.get("NVIDIA_API_KEY")
    if not key:
        return PredictionResult(
            job_id=f"nim-{int(time.time())}",
            status=JobStatus.FAILED,
            backend="boltz2-nim",
            sequences=prediction_input.sequences,
            error_message="NVIDIA_API_KEY not set. Get one at build.nvidia.com",
        )

    start_time = time.time()
    job_id = f"nim-{int(time.time())}"

    async with httpx.AsyncClient(timeout=120) as client:
        payload = {
            "sequences": prediction_input.sequences,
            "diffusion_samples": prediction_input.diffusion_samples,
            "recycling_steps": prediction_input.recycling_steps,
            "sampling_steps": prediction_input.sampling_steps,
        }
        resp = await client.post(
            nim_endpoint,
            json=payload,
            headers={"Authorization": f"Bearer {key}", "Content-Type": "application/json"},
        )
        if resp.status_code != 200:
            return PredictionResult(
                job_id=job_id,
                status=JobStatus.FAILED,
                backend="boltz2-nim",
                sequences=prediction_input.sequences,
                error_message=f"NIM API error {resp.status_code}: {resp.text[:500]}",
                elapsed_seconds=round(time.time() - start_time, 1),
            )

        data = resp.json()
        elapsed = time.time() - start_time

        manifest = ReproducibilityManifest.create(
            backend="boltz2-nim",
            sequences=prediction_input.sequences,
            parameters={
                "diffusion_samples": prediction_input.diffusion_samples,
                "recycling_steps": prediction_input.recycling_steps,
                "sampling_steps": prediction_input.sampling_steps,
                "nim_endpoint": nim_endpoint,
            },
            training_data_source="Boltz-2 (Passaro et al., bioRxiv 2025.06.14.659707)",
        )

        return PredictionResult(
            job_id=job_id,
            status=JobStatus.COMPLETE,
            backend="boltz2-nim",
            sequences=prediction_input.sequences,
            mean_plddt=data.get("mean_plddt"),
            predicted_tm_score=data.get("predicted_tm_score"),
            predicted_affinity=data.get("predicted_affinity"),
            manifest=manifest,
            elapsed_seconds=round(elapsed, 1),
        )


def get_nim_status() -> dict:
    """Check NVIDIA NIM Boltz-2 API availability."""
    has_key = bool(os.environ.get("NVIDIA_API_KEY"))
    return {
        "available": has_key,
        "endpoint": "https://health.api.nvidia.com/v1/biology/mit/boltz-2",
        "setup": "Set NVIDIA_API_KEY environment variable. Get key at build.nvidia.com" if not has_key else "API key configured",
    }

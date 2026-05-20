"""Protenix prediction client.

Protenix-v2 (Apache-2.0, ByteDance, Apr 2026). 464M params, superior
antibody-antigen prediction. Best open Ab-Ag predictor as of May 2026.
"""

from __future__ import annotations

import asyncio
import json
import shutil
import time
import uuid
from pathlib import Path
from typing import Any

from foldcopilot.models.prediction import (
    JobStatus,
    PredictionInput,
    PredictionResult,
    ReproducibilityManifest,
)

PREDICTION_DIR = Path.home() / ".cache" / "foldcopilot" / "predictions"


def _check_protenix_available() -> str | None:
    """Check if Protenix CLI is available."""
    return shutil.which("protenix")


def _write_input_json(sequences: list[str], work_dir: Path) -> Path:
    """Write Protenix input JSON format (AF3-style)."""
    chains = []
    for i, seq in enumerate(sequences):
        chains.append({
            "id": chr(65 + i),
            "sequence": seq,
            "type": "protein",
        })

    input_data = {
        "name": "foldcopilot_job",
        "sequences": chains,
        "modelSeeds": [42],
    }

    input_path = work_dir / "input.json"
    input_path.write_text(json.dumps(input_data, indent=2))
    return input_path


def _parse_output(output_dir: Path) -> dict[str, Any]:
    """Parse Protenix output directory."""
    result: dict[str, Any] = {}

    for ext in ("pdb", "cif"):
        candidates = list(output_dir.rglob(f"*.{ext}"))
        if candidates:
            result[f"output_{ext}_path"] = str(candidates[0])

    confidence_files = list(output_dir.rglob("*confidence*")) + list(
        output_dir.rglob("*summary*")
    )
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
                if "iptm" in data:
                    result["predicted_iptm"] = data["iptm"]
                result["confidence_json_path"] = str(cf)
            except (json.JSONDecodeError, ImportError):
                pass

    if "mean_plddt" not in result and "output_pdb_path" in result:
        try:
            from foldcopilot.clients.boltz2_client import _parse_plddt_from_pdb
            result["mean_plddt"] = _parse_plddt_from_pdb(Path(result["output_pdb_path"]))
        except Exception:
            pass

    return result


async def predict_local(
    prediction_input: PredictionInput,
) -> PredictionResult:
    """Run Protenix prediction locally.

    Protenix-v2 (Apache-2.0, ByteDance, Apr 2026). 464M params, superior
    antibody-antigen prediction. Best open Ab-Ag predictor as of May 2026.
    """
    backend_label = "protenix"
    job_id = str(uuid.uuid4())[:8]
    start_time = time.time()

    protenix_path = _check_protenix_available()
    if not protenix_path:
        return PredictionResult(
            job_id=job_id,
            status=JobStatus.FAILED,
            backend=backend_label,
            sequences=prediction_input.sequences,
            error_message=(
                "Protenix not found in PATH. "
                "Install Protenix: https://github.com/bytedance/Protenix\n"
                "Apache-2.0 license — commercial use OK.\n"
                "Protenix-v2: 464M params, best open antibody-antigen predictor."
            ),
        )

    work_dir = PREDICTION_DIR / f"protenix_{job_id}"
    work_dir.mkdir(parents=True, exist_ok=True)
    output_dir = work_dir / "output"
    output_dir.mkdir(exist_ok=True)

    try:
        input_path = _write_input_json(prediction_input.sequences, work_dir)

        cmd = [
            protenix_path,
            "--input", str(input_path),
            "--output_dir", str(output_dir),
        ]

        if not prediction_input.use_msa:
            cmd.append("--no_msa")

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
                backend=backend_label,
                sequences=prediction_input.sequences,
                error_message=f"Protenix exited with code {process.returncode}: {stderr.decode()[:500]}",
                elapsed_seconds=round(elapsed, 1),
            )

        parsed = _parse_output(output_dir)
        manifest = ReproducibilityManifest.create(
            backend=backend_label,
            sequences=prediction_input.sequences,
            parameters={
                "use_msa": prediction_input.use_msa,
                "diffusion_samples": prediction_input.diffusion_samples,
            },
        )

        manifest_path = output_dir / "reproducibility_manifest.json"
        manifest_path.write_text(manifest.model_dump_json(indent=2))

        return PredictionResult(
            job_id=job_id,
            status=JobStatus.COMPLETE,
            backend=backend_label,
            sequences=prediction_input.sequences,
            output_pdb_path=parsed.get("output_pdb_path"),
            output_cif_path=parsed.get("output_cif_path"),
            confidence_json_path=parsed.get("confidence_json_path"),
            mean_plddt=parsed.get("mean_plddt"),
            predicted_tm_score=parsed.get("predicted_tm_score"),
            manifest=manifest,
            elapsed_seconds=round(elapsed, 1),
        )

    except Exception as e:
        return PredictionResult(
            job_id=job_id,
            status=JobStatus.FAILED,
            backend=backend_label,
            sequences=prediction_input.sequences,
            error_message=str(e),
            elapsed_seconds=round(time.time() - start_time, 1),
        )


def get_status() -> dict:
    protenix_path = _check_protenix_available()
    return {
        "installed": protenix_path is not None,
        "path": protenix_path,
        "gpu_available": _check_gpu(),
        "prediction_dir": str(PREDICTION_DIR),
        "setup_instructions": (
            "Install Protenix: https://github.com/bytedance/Protenix\n"
            "Apache-2.0 license — commercial use OK.\n"
            "Protenix-v2: 464M params, best open antibody-antigen predictor."
        ) if not protenix_path else None,
    }


def _check_gpu() -> bool:
    try:
        import torch
        return torch.cuda.is_available()
    except ImportError:
        return False

"""Chai-1 prediction client.

Wraps Chai-1 (Apache-2.0, Chai Discovery, Sep 2024).
Multi-chain structure prediction with ligand support.
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


def _check_chai_available() -> str | None:
    """Check if Chai-1 CLI / Python module is available."""
    path = shutil.which("chai")
    if path:
        return path
    # Check if chai_lab module is importable
    try:
        import importlib
        importlib.import_module("chai_lab")
        return "chai_lab (python module)"
    except ImportError:
        return None


def _write_input_fasta(sequences: list[str], work_dir: Path) -> Path:
    """Write Chai-1 input FASTA."""
    fasta_path = work_dir / "input.fasta"
    with open(fasta_path, "w") as f:
        for i, seq in enumerate(sequences):
            chain_id = chr(65 + i)
            f.write(f">protein|name=chain_{chain_id}\n{seq}\n")
    return fasta_path


def _parse_output(output_dir: Path) -> dict[str, Any]:
    """Parse Chai-1 output directory."""
    result: dict[str, Any] = {}

    for ext in ("pdb", "cif"):
        candidates = list(output_dir.rglob(f"*.{ext}"))
        if candidates:
            # Chai-1 outputs ranked models; take rank 0
            ranked = sorted(candidates, key=lambda p: p.name)
            result[f"output_{ext}_path"] = str(ranked[0])

    score_files = list(output_dir.rglob("*scores*")) + list(
        output_dir.rglob("*confidence*")
    )
    for sf in score_files:
        if sf.suffix == ".json":
            try:
                data = json.loads(sf.read_text())
                if "plddt" in data:
                    import numpy as np
                    result["mean_plddt"] = round(float(np.mean(data["plddt"])), 1)
                if "ptm" in data:
                    result["predicted_tm_score"] = data["ptm"]
                result["confidence_json_path"] = str(sf)
            except (json.JSONDecodeError, ImportError):
                pass

    if "mean_plddt" not in result and "output_pdb_path" in result:
        try:
            from foldcopilot.clients.boltz2_client import _parse_plddt_from_pdb
            result["mean_plddt"] = _parse_plddt_from_pdb(Path(result["output_pdb_path"]))
        except Exception:
            pass

    return result


async def predict_local(prediction_input: PredictionInput) -> PredictionResult:
    """Run Chai-1 prediction locally."""
    job_id = str(uuid.uuid4())[:8]
    start_time = time.time()

    chai_path = _check_chai_available()
    if not chai_path:
        return PredictionResult(
            job_id=job_id,
            status=JobStatus.FAILED,
            backend="chai1",
            sequences=prediction_input.sequences,
            error_message=(
                "Chai-1 not found. Install with: pip install chai_lab\n"
                "See: https://github.com/chaidiscovery/chai-lab\n"
                "Apache-2.0 license — commercial use OK."
            ),
        )

    work_dir = PREDICTION_DIR / f"chai1_{job_id}"
    work_dir.mkdir(parents=True, exist_ok=True)
    output_dir = work_dir / "output"
    output_dir.mkdir(exist_ok=True)

    try:
        input_path = _write_input_fasta(prediction_input.sequences, work_dir)

        # Try CLI first, fall back to Python API
        cli_path = shutil.which("chai")
        if cli_path:
            cmd = [
                cli_path,
                "fold",
                str(input_path),
                "--output-dir", str(output_dir),
                "--num-diffusion-samples", str(prediction_input.diffusion_samples),
            ]

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
                    backend="chai1",
                    sequences=prediction_input.sequences,
                    error_message=f"Chai-1 exited with code {process.returncode}: {stderr.decode()[:500]}",
                    elapsed_seconds=round(elapsed, 1),
                )
        else:
            # Python API path
            cmd = [
                "python", "-c",
                f"from chai_lab.chai1 import run_inference; "
                f"run_inference(fasta_file='{input_path}', output_dir='{output_dir}', "
                f"num_trunk_recycles={prediction_input.recycling_steps}, "
                f"num_diffn_timesteps={prediction_input.sampling_steps})"
            ]

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
                    backend="chai1",
                    sequences=prediction_input.sequences,
                    error_message=f"Chai-1 Python API failed: {stderr.decode()[:500]}",
                    elapsed_seconds=round(elapsed, 1),
                )

        elapsed = time.time() - start_time
        parsed = _parse_output(output_dir)

        manifest = ReproducibilityManifest.create(
            backend="chai1",
            sequences=prediction_input.sequences,
            parameters={
                "recycling_steps": prediction_input.recycling_steps,
                "sampling_steps": prediction_input.sampling_steps,
                "diffusion_samples": prediction_input.diffusion_samples,
                "use_msa": prediction_input.use_msa,
            },
        )

        manifest_path = output_dir / "reproducibility_manifest.json"
        manifest_path.write_text(manifest.model_dump_json(indent=2))

        return PredictionResult(
            job_id=job_id,
            status=JobStatus.COMPLETE,
            backend="chai1",
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
            backend="chai1",
            sequences=prediction_input.sequences,
            error_message=str(e),
            elapsed_seconds=round(time.time() - start_time, 1),
        )


def get_status() -> dict:
    chai_path = _check_chai_available()
    return {
        "installed": chai_path is not None,
        "path": chai_path,
        "gpu_available": _check_gpu(),
        "prediction_dir": str(PREDICTION_DIR),
        "setup_instructions": (
            "Install Chai-1: pip install chai_lab\n"
            "See: https://github.com/chaidiscovery/chai-lab\n"
            "Apache-2.0 license — commercial use OK."
        ) if not chai_path else None,
    }


def _check_gpu() -> bool:
    try:
        import torch
        return torch.cuda.is_available()
    except ImportError:
        return False

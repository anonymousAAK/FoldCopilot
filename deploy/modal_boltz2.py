"""
Modal deploy script for Boltz-2 GPU inference.

Usage:
    # First time: set up Modal account and token
    pip install modal
    modal token new

    # Deploy as a persistent endpoint
    modal deploy deploy/modal_boltz2.py

    # Or run a one-off prediction
    modal run deploy/modal_boltz2.py --sequence "MKTVRQERLKS..."

    # Call the deployed endpoint from Python
    import modal
    f = modal.Function.from_name("foldcopilot-boltz2", "predict")
    result = f.remote(sequence="MKTVRQERLKS...", name="my_protein")

Cost estimate:
    A10G: ~$0.60/hr  -> ~$0.003 per 20s prediction
    A100: ~$2.78/hr  -> ~$0.015 per 20s prediction
"""

import modal
import subprocess
import tempfile
import json
from pathlib import Path

app = modal.App("foldcopilot-boltz2")

boltz_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("boltz", "biopython")
    .run_commands("python -c 'import boltz; print(boltz.__version__)'")
)


@app.function(
    image=boltz_image,
    gpu=modal.gpu.A10G(),          # swap to modal.gpu.A100() for faster inference
    timeout=600,                    # 10 min max per prediction
    memory=32768,                   # 32 GB RAM
    volumes={"/cache": modal.Volume.from_name("boltz-cache", create_if_missing=True)},
)
def predict(
    sequence: str,
    name: str = "query",
    recycling_steps: int = 3,
    sampling_steps: int = 200,
    diffusion_samples: int = 1,
) -> dict:
    """Run Boltz-2 structure prediction on a single protein sequence.

    Args:
        sequence: Amino acid sequence (one-letter codes).
        name: Label for the prediction job.
        recycling_steps: Number of recycling iterations.
        sampling_steps: Number of diffusion sampling steps.
        diffusion_samples: Number of independent samples to generate.

    Returns:
        dict with keys: name, pdb_content, confidence_scores, output_dir
    """
    work_dir = Path(tempfile.mkdtemp())
    fasta_path = work_dir / f"{name}.fasta"

    # Write input FASTA
    fasta_path.write_text(f">{name}\n{sequence}\n")

    # Run boltz predict
    cmd = [
        "boltz", "predict",
        str(fasta_path),
        "--out_dir", str(work_dir / "output"),
        "--recycling_steps", str(recycling_steps),
        "--sampling_steps", str(sampling_steps),
        "--diffusion_samples", str(diffusion_samples),
        "--output_format", "pdb",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=540)

    if result.returncode != 0:
        return {
            "error": True,
            "name": name,
            "stderr": result.stderr[-2000:] if result.stderr else "",
            "stdout": result.stdout[-2000:] if result.stdout else "",
        }

    # Collect output PDB files
    output_dir = work_dir / "output"
    pdb_files = list(output_dir.rglob("*.pdb"))

    if not pdb_files:
        return {"error": True, "name": name, "message": "No PDB output generated"}

    # Read the first (best) PDB
    best_pdb = sorted(pdb_files)[0]
    pdb_content = best_pdb.read_text()

    # Cache the result
    cache_dir = Path("/cache") / name
    cache_dir.mkdir(parents=True, exist_ok=True)
    (cache_dir / f"{name}.pdb").write_text(pdb_content)

    return {
        "error": False,
        "name": name,
        "pdb_content": pdb_content,
        "pdb_files": [str(p) for p in pdb_files],
        "num_residues": len(sequence),
        "num_samples": diffusion_samples,
    }


@app.local_entrypoint()
def main(
    sequence: str = "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG",
    name: str = "example",
):
    """CLI entrypoint for `modal run deploy/modal_boltz2.py`."""
    print(f"Predicting structure for '{name}' ({len(sequence)} residues)...")
    result = predict.remote(sequence=sequence, name=name)

    if result.get("error"):
        print(f"Prediction failed: {result}")
    else:
        print(f"Prediction complete. {len(result['pdb_files'])} PDB file(s) generated.")
        # Save locally
        out_path = Path(f"{name}_boltz2.pdb")
        out_path.write_text(result["pdb_content"])
        print(f"Best structure saved to {out_path}")

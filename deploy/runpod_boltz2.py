"""
RunPod serverless handler for Boltz-2 structure prediction.

Setup:
    1. Create a RunPod account at https://www.runpod.io
    2. Build a Docker image from this handler:

        FROM pytorch/pytorch:2.3.0-cuda12.1-cudnn8-runtime
        RUN pip install boltz runpod biopython
        COPY runpod_boltz2.py /handler.py
        CMD ["python", "/handler.py"]

    3. Push to Docker Hub or RunPod registry:
        docker build -t yourusername/boltz2-runpod:latest -f Dockerfile .
        docker push yourusername/boltz2-runpod:latest

    4. Create a Serverless Endpoint on RunPod:
        - Template: your pushed image
        - GPU: A10G (24GB) or A100 (40/80GB)
        - Max workers: 1-3 (scale to taste)

    5. Send requests:
        import runpod
        runpod.api_key = "YOUR_API_KEY"
        endpoint = runpod.Endpoint("YOUR_ENDPOINT_ID")
        result = endpoint.run_sync({
            "input": {
                "sequence": "MKTVRQERLKS...",
                "name": "my_protein"
            }
        })

Cost estimate:
    A10G: ~$0.00031/s -> ~$0.006 per 20s prediction
    A100: ~$0.00076/s -> ~$0.015 per 20s prediction
    (Only billed while actively processing)
"""

import runpod
import subprocess
import tempfile
from pathlib import Path


def handler(event: dict) -> dict:
    """RunPod serverless handler for Boltz-2 prediction.

    Input schema:
        {
            "sequence": str,         # Required. Amino acid sequence.
            "name": str,             # Optional. Job label. Default: "query".
            "recycling_steps": int,  # Optional. Default: 3.
            "sampling_steps": int,   # Optional. Default: 200.
            "diffusion_samples": int # Optional. Default: 1.
        }

    Returns:
        {
            "pdb_content": str,      # PDB file content of best prediction.
            "num_residues": int,
            "num_samples": int,
            "pdb_files": list[str]
        }
    """
    job_input = event["input"]

    sequence = job_input["sequence"]
    name = job_input.get("name", "query")
    recycling_steps = job_input.get("recycling_steps", 3)
    sampling_steps = job_input.get("sampling_steps", 200)
    diffusion_samples = job_input.get("diffusion_samples", 1)

    # Validate sequence
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    cleaned = sequence.upper().strip()
    if not cleaned or not all(c in valid_aa for c in cleaned):
        return {"error": f"Invalid amino acid sequence. Use standard one-letter codes."}

    work_dir = Path(tempfile.mkdtemp())
    fasta_path = work_dir / f"{name}.fasta"
    fasta_path.write_text(f">{name}\n{cleaned}\n")

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
            "error": "Boltz-2 prediction failed",
            "stderr": result.stderr[-2000:] if result.stderr else "",
        }

    output_dir = work_dir / "output"
    pdb_files = sorted(output_dir.rglob("*.pdb"))

    if not pdb_files:
        return {"error": "No PDB output generated"}

    best_pdb = pdb_files[0]
    pdb_content = best_pdb.read_text()

    return {
        "pdb_content": pdb_content,
        "name": name,
        "num_residues": len(cleaned),
        "num_samples": diffusion_samples,
        "pdb_files": [p.name for p in pdb_files],
    }


if __name__ == "__main__":
    runpod.serverless.start({"handler": handler})

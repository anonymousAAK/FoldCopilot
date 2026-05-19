# FoldCopilot Deploy Manifests

Reference deployment configs for running Boltz-2 inference on GPU. FoldCopilot is client-only -- it connects to these backends but never hosts them.

## Options

| Option | File | GPU | Best for |
|---|---|---|---|
| **Modal** | `modal_boltz2.py` | A10G / A100 | Quick experiments, pay-per-second, no infra to manage |
| **RunPod** | `runpod_boltz2.py` | A10G / A100 | Serverless with longer cold-start tolerance |
| **NVIDIA NIM** | `nim_boltz2.sh` | Any NVIDIA >=24GB | Self-hosted, on-prem, or persistent cloud VM |

## Cost Estimates

| Platform | GPU | Per-prediction (20s) | Hourly |
|---|---|---|---|
| Modal | A10G | ~$0.003 | ~$0.60 |
| Modal | A100 | ~$0.015 | ~$2.78 |
| RunPod | A10G | ~$0.006 | ~$1.10 |
| RunPod | A100 | ~$0.015 | ~$2.70 |
| NIM (cloud VM) | A10G | VM cost | ~$0.60-1.00 (spot) |
| NIM (on-prem) | Any | Electricity | -- |

## Prerequisites

**All options** require:
- Python 3.11+
- NVIDIA GPU with >=24 GB VRAM (A10G, L40S, A100, H100)

**Modal** (`modal_boltz2.py`):
```bash
pip install modal
modal token new
modal deploy deploy/modal_boltz2.py    # persistent endpoint
modal run deploy/modal_boltz2.py       # one-off run
```

**RunPod** (`runpod_boltz2.py`):
```bash
pip install runpod
# Build Docker image, push, create serverless endpoint on runpod.io
# See docstring in runpod_boltz2.py for full instructions
```

**NVIDIA NIM** (`nim_boltz2.sh`):
```bash
export NGC_API_KEY="your-key"   # from https://org.ngc.nvidia.com/setup/api-key
bash deploy/nim_boltz2.sh
# Serves REST API on http://localhost:8000
```

## Connecting to FoldCopilot

Once a backend is running, point FoldCopilot at it:

```python
# FoldCopilot auto-detects local boltz installation.
# For remote endpoints, set the environment variable:
export BOLTZ2_ENDPOINT="http://localhost:8000"
```

## These are reference configs

Not production infrastructure. Adjust GPU types, timeouts, scaling, and auth to your needs.

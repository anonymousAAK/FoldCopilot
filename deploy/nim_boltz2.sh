#!/usr/bin/env bash
# ------------------------------------------------------------------
# NVIDIA NIM Boltz-2 self-hosted container
#
# Prerequisites:
#   - Docker with NVIDIA Container Toolkit (nvidia-docker2)
#   - NVIDIA GPU with >=24 GB VRAM (A10G, L40S, A100, H100)
#   - NGC API key: https://org.ngc.nvidia.com/setup/api-key
#
# Usage:
#   export NGC_API_KEY="your-ngc-api-key"
#   bash deploy/nim_boltz2.sh
#
# Cost estimate (self-hosted):
#   Cloud VM with A10G: ~$0.60-1.00/hr (AWS/GCP/Azure spot)
#   Cloud VM with A100: ~$2.50-3.50/hr
#   On-prem: electricity only
# ------------------------------------------------------------------

set -euo pipefail

# ---- Configuration ----
NIM_IMAGE="nvcr.io/nim/mit/boltz2:1.6.0"
CONTAINER_NAME="boltz2-nim"
HOST_PORT="${BOLTZ2_PORT:-8000}"
MODEL_CACHE="${BOLTZ2_CACHE_DIR:-$HOME/.cache/boltz2-nim}"

# ---- Preflight checks ----
if [ -z "${NGC_API_KEY:-}" ]; then
    echo "ERROR: NGC_API_KEY is not set."
    echo "Get your key at https://org.ngc.nvidia.com/setup/api-key"
    echo "  export NGC_API_KEY=\"your-key\""
    exit 1
fi

if ! command -v docker &> /dev/null; then
    echo "ERROR: docker is not installed."
    exit 1
fi

if ! docker info 2>/dev/null | grep -q "Runtimes.*nvidia"; then
    echo "WARNING: NVIDIA container runtime may not be configured."
    echo "See: https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html"
fi

# ---- Login to NGC registry ----
echo "Logging in to NGC container registry..."
echo "${NGC_API_KEY}" | docker login nvcr.io -u '$oauthtoken' --password-stdin

# ---- Pull image ----
echo "Pulling ${NIM_IMAGE}..."
docker pull "${NIM_IMAGE}"

# ---- Create cache directory ----
mkdir -p "${MODEL_CACHE}"

# ---- Stop existing container if running ----
if docker ps -a --format '{{.Names}}' | grep -q "^${CONTAINER_NAME}$"; then
    echo "Stopping existing ${CONTAINER_NAME} container..."
    docker rm -f "${CONTAINER_NAME}" > /dev/null 2>&1
fi

# ---- Run container ----
echo "Starting Boltz-2 NIM on port ${HOST_PORT}..."
docker run -d \
    --name "${CONTAINER_NAME}" \
    --gpus all \
    -p "${HOST_PORT}:8000" \
    -v "${MODEL_CACHE}:/opt/nim/.cache" \
    -e "NGC_API_KEY=${NGC_API_KEY}" \
    "${NIM_IMAGE}"

echo ""
echo "Container '${CONTAINER_NAME}' started."
echo "Waiting for health check..."

# ---- Wait for readiness ----
for i in $(seq 1 60); do
    if curl -sf "http://localhost:${HOST_PORT}/v1/health/ready" > /dev/null 2>&1; then
        echo "Boltz-2 NIM is ready on http://localhost:${HOST_PORT}"
        break
    fi
    if [ "$i" -eq 60 ]; then
        echo "WARNING: Health check timed out after 60s. Check logs:"
        echo "  docker logs ${CONTAINER_NAME}"
        exit 1
    fi
    sleep 1
done

# ---- Usage example ----
cat << 'EOF'

--- Example: single-sequence prediction ---

curl -s http://localhost:8000/v1/structure/predict \
  -H "Content-Type: application/json" \
  -d '{
    "sequences": [
      {
        "protein": {
          "id": "A",
          "sequence": "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"
        }
      }
    ],
    "recycling_steps": 3,
    "sampling_steps": 200,
    "diffusion_samples": 1
  }' \
  -o prediction.zip

unzip prediction.zip -d prediction_output/
echo "Structure saved to prediction_output/"

--- Example: multi-chain complex ---

curl -s http://localhost:8000/v1/structure/predict \
  -H "Content-Type: application/json" \
  -d '{
    "sequences": [
      {"protein": {"id": "A", "sequence": "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQ"}},
      {"protein": {"id": "B", "sequence": "DIAYLRSLGYNIVATPRGYVLAGG"}}
    ]
  }' \
  -o complex.zip

--- Logs & cleanup ---

docker logs -f boltz2-nim          # stream logs
docker stop boltz2-nim             # stop
docker rm boltz2-nim               # remove

EOF

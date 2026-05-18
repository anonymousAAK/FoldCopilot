"""AlphaFold Database API client.

Wraps the AFDB REST API (alphafold.ebi.ac.uk/api) for structure lookup,
pLDDT scores, and PAE matrices.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import httpx
import numpy as np

AFDB_BASE = "https://alphafold.ebi.ac.uk/api"
AFDB_FILES = "https://alphafold.ebi.ac.uk/files"

# Cache directory — content-addressed by (uniprot_id, data_type)
_CACHE_DIR = Path.home() / ".cache" / "foldcopilot" / "afdb"


def _cache_key(uniprot_id: str, data_type: str) -> Path:
    h = hashlib.sha256(f"{uniprot_id}:{data_type}".encode()).hexdigest()[:16]
    return _CACHE_DIR / f"{uniprot_id}_{data_type}_{h}.json"


def _read_cache(key: Path) -> Any | None:
    if key.exists():
        return json.loads(key.read_text())
    return None


def _write_cache(key: Path, data: Any) -> None:
    key.parent.mkdir(parents=True, exist_ok=True)
    key.write_text(json.dumps(data))


async def get_prediction_metadata(
    uniprot_id: str, *, client: httpx.AsyncClient | None = None
) -> dict:
    """Fetch AFDB prediction metadata for a UniProt accession."""
    cache_key = _cache_key(uniprot_id, "metadata")
    cached = _read_cache(cache_key)
    if cached is not None:
        return cached

    own_client = client is None
    if own_client:
        client = httpx.AsyncClient(timeout=30)
    try:
        url = f"{AFDB_BASE}/prediction/{uniprot_id}"
        resp = await client.get(url)
        resp.raise_for_status()
        data = resp.json()
        # AFDB returns a list; take first entry
        if isinstance(data, list):
            data = data[0]
        _write_cache(cache_key, data)
        return data
    finally:
        if own_client:
            await client.aclose()


async def get_plddt_scores(
    uniprot_id: str, *, client: httpx.AsyncClient | None = None
) -> list[float]:
    """Fetch per-residue pLDDT scores from AFDB CIF confidence data.

    Uses the summary JSON endpoint which contains per-residue confidence.
    Falls back to parsing the CIF B-factor column.
    """
    cache_key = _cache_key(uniprot_id, "plddt")
    cached = _read_cache(cache_key)
    if cached is not None:
        return cached

    own_client = client is None
    if own_client:
        client = httpx.AsyncClient(timeout=60)
    try:
        # Try the confidence summary endpoint first
        meta = await get_prediction_metadata(uniprot_id, client=client)
        # Fetch the PDB file — pLDDT is stored in the B-factor column
        pdb_url = meta.get("pdbUrl")
        if not pdb_url:
            entry_id = meta.get("entryId", f"AF-{uniprot_id}-F1")
            pdb_url = f"{AFDB_FILES}/{entry_id}-model_v4.pdb"

        resp = await client.get(pdb_url)
        resp.raise_for_status()
        pdb_text = resp.text

        # Parse B-factor from ATOM records (pLDDT is in columns 61-66)
        plddt_by_residue: dict[int, float] = {}
        for line in pdb_text.splitlines():
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                res_seq = int(line[22:26].strip())
                b_factor = float(line[60:66].strip())
                plddt_by_residue[res_seq] = b_factor

        scores = [plddt_by_residue[k] for k in sorted(plddt_by_residue)]
        _write_cache(cache_key, scores)
        return scores
    finally:
        if own_client:
            await client.aclose()


async def get_pae_matrix(
    uniprot_id: str, *, client: httpx.AsyncClient | None = None
) -> np.ndarray | None:
    """Fetch the PAE matrix from AFDB. Returns NxN numpy array or None."""
    cache_key = _cache_key(uniprot_id, "pae")
    cached = _read_cache(cache_key)
    if cached is not None:
        return np.array(cached)

    own_client = client is None
    if own_client:
        client = httpx.AsyncClient(timeout=60)
    try:
        meta = await get_prediction_metadata(uniprot_id, client=client)
        pae_url = meta.get("paeImageUrl", "").replace("png", "json")
        if not pae_url:
            entry_id = meta.get("entryId", f"AF-{uniprot_id}-F1")
            pae_url = (
                f"{AFDB_FILES}/{entry_id}-predicted_aligned_error_v4.json"
            )

        resp = await client.get(pae_url)
        if resp.status_code == 404:
            return None
        resp.raise_for_status()

        pae_data = resp.json()
        # AFDB PAE JSON format: list with one dict containing
        # "predicted_aligned_error" as a 2D list
        if isinstance(pae_data, list) and len(pae_data) > 0:
            matrix = pae_data[0].get(
                "predicted_aligned_error",
                pae_data[0].get("pae"),
            )
        elif isinstance(pae_data, dict):
            matrix = pae_data.get(
                "predicted_aligned_error", pae_data.get("pae")
            )
        else:
            return None

        if matrix is None:
            return None

        arr = np.array(matrix, dtype=np.float32)
        _write_cache(cache_key, arr.tolist())
        return arr
    finally:
        if own_client:
            await client.aclose()

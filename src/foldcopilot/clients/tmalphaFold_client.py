"""TMalphaFold / OPM client — membrane orientation data for transmembrane proteins.

TMalphaFold provides predicted transmembrane topology and membrane positioning
for AlphaFold structures. OPM (Orientations of Proteins in Membranes) provides
experimentally validated membrane boundaries.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import httpx

TMAF_BASE = "https://tmalphafold.ttk.hu/api"
OPM_BASE = "https://opm.phar.umich.edu/api"
CACHE_DIR = Path.home() / ".cache" / "foldcopilot" / "tmalphaFold"


def _cache_path(uniprot_id: str, data_type: str) -> Path:
    key = hashlib.sha256(f"{uniprot_id}:{data_type}".encode()).hexdigest()[:16]
    return CACHE_DIR / f"{key}.json"


async def get_membrane_topology(uniprot_id: str) -> dict[str, Any]:
    """Get predicted transmembrane topology from TMalphaFold.

    Returns:
        Transmembrane segments, membrane insertion angle, and topology type.
    """
    cache = _cache_path(uniprot_id, "topology")
    if cache.exists():
        import json
        return json.loads(cache.read_text())

    async with httpx.AsyncClient(timeout=30) as client:
        try:
            resp = await client.get(f"{TMAF_BASE}/protein/{uniprot_id}")
            resp.raise_for_status()
            data = resp.json()
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                return {
                    "uniprot_id": uniprot_id,
                    "found": False,
                    "message": f"No TMalphaFold data for {uniprot_id}. May not be a transmembrane protein.",
                }
            return {"error": f"TMalphaFold API error: {e.response.status_code}"}
        except httpx.RequestError as e:
            return {"error": f"TMalphaFold request failed: {str(e)}"}

    tm_segments = data.get("tm_segments", [])
    result = {
        "uniprot_id": uniprot_id,
        "found": True,
        "topology_type": data.get("type", "unknown"),
        "tm_segments": [
            {
                "start": seg.get("start"),
                "end": seg.get("end"),
                "type": seg.get("type", "transmembrane"),
            }
            for seg in tm_segments
        ],
        "n_tm_segments": len(tm_segments),
        "membrane_insertion_angle": data.get("tilt_angle"),
        "source": "TMalphaFold",
    }

    # Cache result
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    import json
    cache.write_text(json.dumps(result, indent=2))

    return result


async def get_opm_orientation(uniprot_id: str) -> dict[str, Any]:
    """Get experimentally validated membrane orientation from OPM.

    Falls back gracefully if protein not in OPM database.
    """
    async with httpx.AsyncClient(timeout=30) as client:
        try:
            resp = await client.get(f"{OPM_BASE}/search", params={"query": uniprot_id})
            resp.raise_for_status()
            data = resp.json()
        except (httpx.HTTPStatusError, httpx.RequestError):
            return {
                "uniprot_id": uniprot_id,
                "found": False,
                "message": "OPM data not available.",
            }

    if not data or (isinstance(data, list) and len(data) == 0):
        return {
            "uniprot_id": uniprot_id,
            "found": False,
            "message": f"No OPM entry for {uniprot_id}.",
        }

    entry = data[0] if isinstance(data, list) else data
    return {
        "uniprot_id": uniprot_id,
        "found": True,
        "pdb_id": entry.get("pdbid"),
        "type": entry.get("type_subunit", "unknown"),
        "thickness": entry.get("thickness"),
        "tilt_angle": entry.get("tilt"),
        "source": "OPM",
    }


async def get_membrane_context(uniprot_id: str) -> dict[str, Any]:
    """Combined membrane context: TMalphaFold topology + OPM orientation.

    Use this for GPCR and other transmembrane protein analysis.
    """
    topology = await get_membrane_topology(uniprot_id)
    opm = await get_opm_orientation(uniprot_id)

    is_tm = topology.get("found", False) or opm.get("found", False)

    return {
        "uniprot_id": uniprot_id,
        "is_transmembrane": is_tm,
        "tmalphaFold": topology,
        "opm": opm,
        "interpretation": (
            f"Predicted {topology.get('n_tm_segments', 0)} transmembrane segments. "
            f"Topology type: {topology.get('topology_type', 'unknown')}."
            if topology.get("found") else
            "No transmembrane topology predicted. May be a soluble protein."
        ),
    }

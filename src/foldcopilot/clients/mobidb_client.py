"""MobiDB API client — fetches disorder/mobility annotations.

MobiDB (https://mobidb.bio.unipd.it) aggregates disorder predictions and
curated annotations from multiple sources.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import httpx

MOBIDB_API = "https://mobidb.bio.unipd.it/api/download"

_CACHE_DIR = Path.home() / ".cache" / "foldcopilot" / "mobidb"


def _cache_path(uniprot_id: str) -> Path:
    h = hashlib.sha256(uniprot_id.encode()).hexdigest()[:12]
    return _CACHE_DIR / f"{uniprot_id}_{h}.json"


async def get_mobidb_regions(
    uniprot_id: str, *, client: httpx.AsyncClient | None = None
) -> list[dict]:
    """Return MobiDB disorder annotations for a UniProt accession.

    Returns list of dicts with keys: start, end, annotation, source.
    Empty list if no MobiDB entry exists.
    """
    cp = _cache_path(uniprot_id)
    if cp.exists():
        return json.loads(cp.read_text())

    own_client = client is None
    if own_client:
        client = httpx.AsyncClient(timeout=30)
    try:
        resp = await client.get(
            f"{MOBIDB_API}?acc={uniprot_id}&format=json",
        )
        if resp.status_code == 404:
            regions: list[dict] = []
            _cache(cp, regions)
            return regions

        resp.raise_for_status()
        data = resp.json()

        regions = []
        # MobiDB returns a dict with accession as key, or a list
        entries = data if isinstance(data, list) else [data]

        for entry in entries:
            # Consensus disorder
            for feature in entry.get("mobidb_consensus", {}).get("disorder", {}).get("predictors", []):
                for region in feature.get("regions", []):
                    regions.append({
                        "start": region[0],
                        "end": region[1],
                        "annotation": "disorder",
                        "source": f"mobidb_{feature.get('method', 'consensus')}",
                    })

            # Also check curated data
            for feature in entry.get("mobidb_consensus", {}).get("disorder", {}).get("curated", []):
                for region in feature.get("regions", []):
                    regions.append({
                        "start": region[0],
                        "end": region[1],
                        "annotation": "disorder_curated",
                        "source": f"mobidb_{feature.get('source', 'curated')}",
                    })

        # Deduplicate
        seen = set()
        unique = []
        for r in regions:
            key = (r["start"], r["end"])
            if key not in seen:
                seen.add(key)
                unique.append(r)

        _cache(cp, unique)
        return unique
    finally:
        if own_client:
            await client.aclose()


def _cache(path: Path, data: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))

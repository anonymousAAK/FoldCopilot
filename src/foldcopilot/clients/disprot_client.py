"""DisProt API client — fetches intrinsically disordered region annotations.

DisProt (https://disprot.org) is the curated ground-truth database for IDRs.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import httpx

DISPROT_API = "https://disprot.org/api"

_CACHE_DIR = Path.home() / ".cache" / "foldcopilot" / "disprot"


def _cache_path(uniprot_id: str) -> Path:
    h = hashlib.sha256(uniprot_id.encode()).hexdigest()[:12]
    return _CACHE_DIR / f"{uniprot_id}_{h}.json"


async def get_disprot_regions(
    uniprot_id: str, *, client: httpx.AsyncClient | None = None
) -> list[dict]:
    """Return DisProt disorder annotations for a UniProt accession.

    Returns list of dicts with keys: start, end, annotation, evidence_type.
    When available (DisProt 2026 / IDPO), also includes optional keys:
    - term_namespace: IDPO ontology namespace (e.g. "IDPO:00076")
    - method: MIADE experimental/computational method
    - conditions: MIADE experimental conditions
    Empty list if no DisProt entry exists.
    """
    cp = _cache_path(uniprot_id)
    if cp.exists():
        return json.loads(cp.read_text())

    own_client = client is None
    if own_client:
        client = httpx.AsyncClient(timeout=30)
    try:
        # DisProt search by UniProt accession
        resp = await client.get(
            f"{DISPROT_API}/{uniprot_id}",
            headers={"Accept": "application/json"},
        )
        if resp.status_code == 404:
            regions: list[dict] = []
            _cache(cp, regions)
            return regions

        resp.raise_for_status()
        data = resp.json()

        regions = []
        for region in data.get("disprot_consensus", {}).get("structural_state", []):
            if region.get("type") in ("disorder", "D"):
                entry: dict = {
                    "start": region["start"],
                    "end": region["end"],
                    "annotation": "disorder",
                    "evidence_type": "disprot_consensus",
                }
                if region.get("term_namespace"):
                    entry["term_namespace"] = region["term_namespace"]
                regions.append(entry)

        # Also check individual annotations
        for ann in data.get("regions", []):
            term = ann.get("term_name", "").lower()
            if "disorder" in term or "flexible" in term:
                entry = {
                    "start": ann["start"],
                    "end": ann["end"],
                    "annotation": ann.get("term_name", "disorder"),
                    "evidence_type": ann.get("evidence_code", "unknown"),
                }
                # DisProt 2026: IDPO ontology namespace
                if ann.get("term_namespace"):
                    entry["term_namespace"] = ann["term_namespace"]
                # DisProt 2026: MIADE evidence details
                if ann.get("method"):
                    entry["method"] = ann["method"]
                if ann.get("conditions"):
                    entry["conditions"] = ann["conditions"]
                regions.append(entry)

        # Deduplicate by (start, end)
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

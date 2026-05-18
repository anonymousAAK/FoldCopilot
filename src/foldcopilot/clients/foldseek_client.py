"""Foldseek API client — structural similarity search.

Wraps the Foldseek webserver REST API (search.foldseek.com/api).
No existing MCP wrapper for Foldseek exists as of May 2026.
"""

from __future__ import annotations

import asyncio
import hashlib
import json
from pathlib import Path
from typing import Any

import httpx

FOLDSEEK_API = "https://search.foldseek.com/api"

_CACHE_DIR = Path.home() / ".cache" / "foldcopilot" / "foldseek"


def _cache_key(query_hash: str) -> Path:
    return _CACHE_DIR / f"{query_hash}.json"


def _read_cache(key: Path) -> Any | None:
    if key.exists():
        return json.loads(key.read_text())
    return None


def _write_cache(key: Path, data: Any) -> None:
    key.parent.mkdir(parents=True, exist_ok=True)
    key.write_text(json.dumps(data))


async def submit_search(
    pdb_content: str,
    *,
    databases: list[str] | None = None,
    mode: str = "3diaa",
    client: httpx.AsyncClient | None = None,
) -> str:
    """Submit a structure search to Foldseek. Returns a ticket ID for polling.

    Args:
        pdb_content: PDB file content as string.
        databases: List of databases to search. Defaults to ["afdb50", "afdb-swissprot", "pdb100"].
        mode: Search mode — "3diaa" (default, fast) or "tmalign" (slower, more sensitive).
    """
    if databases is None:
        databases = ["afdb50", "afdb-swissprot", "pdb100"]

    own_client = client is None
    if own_client:
        client = httpx.AsyncClient(timeout=60)
    try:
        # Foldseek API expects multipart form data
        files = {"q": ("query.pdb", pdb_content, "text/plain")}
        data = {
            "mode": mode,
            "database[]": databases,
        }
        resp = await client.post(
            f"{FOLDSEEK_API}/ticket",
            files=files,
            data=data,
        )
        resp.raise_for_status()
        result = resp.json()
        ticket_id = result.get("id")
        if not ticket_id:
            raise ValueError(f"Foldseek did not return a ticket ID: {result}")
        return ticket_id
    finally:
        if own_client:
            await client.aclose()


async def poll_result(
    ticket_id: str,
    *,
    max_wait: int = 300,
    poll_interval: int = 5,
    client: httpx.AsyncClient | None = None,
) -> dict:
    """Poll Foldseek for search results. Blocks until complete or timeout.

    Returns the full result dict with alignments.
    """
    # Check cache
    ck = _cache_key(ticket_id)
    cached = _read_cache(ck)
    if cached is not None:
        return cached

    own_client = client is None
    if own_client:
        client = httpx.AsyncClient(timeout=30)
    try:
        elapsed = 0
        while elapsed < max_wait:
            resp = await client.get(f"{FOLDSEEK_API}/ticket/{ticket_id}")
            resp.raise_for_status()
            data = resp.json()

            status = data.get("status")
            if status == "COMPLETE":
                # Fetch actual results
                result_resp = await client.get(
                    f"{FOLDSEEK_API}/result/{ticket_id}/0",
                )
                result_resp.raise_for_status()
                results = result_resp.json()
                _write_cache(ck, results)
                return results
            elif status == "ERROR":
                raise RuntimeError(f"Foldseek search failed: {data.get('error', 'unknown')}")
            elif status in ("PENDING", "RUNNING"):
                await asyncio.sleep(poll_interval)
                elapsed += poll_interval
            else:
                raise RuntimeError(f"Unknown Foldseek status: {status}")

        raise TimeoutError(f"Foldseek search timed out after {max_wait}s (ticket: {ticket_id})")
    finally:
        if own_client:
            await client.aclose()


async def search_pdb(
    pdb_content: str,
    *,
    databases: list[str] | None = None,
    mode: str = "3diaa",
    client: httpx.AsyncClient | None = None,
) -> dict:
    """Submit a Foldseek search and wait for results. Convenience wrapper.

    Returns parsed results with alignments.
    """
    own_client = client is None
    if own_client:
        client = httpx.AsyncClient(timeout=60)
    try:
        ticket = await submit_search(
            pdb_content, databases=databases, mode=mode, client=client
        )
        return await poll_result(ticket, client=client)
    finally:
        if own_client:
            await client.aclose()


def parse_alignments(raw_results: dict, max_hits: int = 20) -> list[dict]:
    """Parse Foldseek raw results into a clean list of hits.

    Returns top hits with: target, probability, evalue, score, TM-score,
    aligned length, sequence identity, description.
    """
    hits = []
    # Foldseek returns results per database, each with alignments
    results_list = raw_results.get("results", [])
    if not isinstance(results_list, list):
        results_list = [results_list]

    for db_result in results_list:
        db_name = db_result.get("db", "unknown")
        for alignment in db_result.get("alignments", [])[:max_hits]:
            hit = {
                "database": db_name,
                "target": alignment.get("target", ""),
                "target_description": alignment.get("tDescription", ""),
                "probability": alignment.get("prob"),
                "evalue": alignment.get("eval"),
                "score": alignment.get("score"),
                "tm_score": alignment.get("tmScore"),
                "aligned_length": alignment.get("alnLength"),
                "sequence_identity": alignment.get("seqId"),
                "query_start": alignment.get("qStartPos"),
                "query_end": alignment.get("qEndPos"),
                "target_start": alignment.get("tStartPos"),
                "target_end": alignment.get("tEndPos"),
            }
            hits.append(hit)

    # Sort by evalue (lower better), then TM-score (higher better)
    hits.sort(key=lambda h: (h.get("evalue") or 999, -(h.get("tm_score") or 0)))
    return hits[:max_hits]

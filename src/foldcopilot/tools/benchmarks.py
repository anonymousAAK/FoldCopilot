"""Benchmarking harness — evaluate prediction quality against known structures.

v0.9: Public eval suite supporting CASP16 monomers, DisProt hallucination set,
and custom benchmark datasets. Ships results as structured reports.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

import numpy as np

from foldcopilot.tools.ensemble import (
    compute_ca_distances,
    compute_rmsd,
    parse_ca_coords,
    parse_plddt_from_pdb,
)


# Built-in benchmark datasets
BENCHMARK_DATASETS = {
    "disprot_hallucination": {
        "name": "DisProt AF3 Hallucination Set",
        "description": (
            "72 DisProt proteins where AF3 hallucinations are documented. "
            "Based on arXiv 2510.15939: 22% of IDR residues falsely predicted as ordered."
        ),
        "source": "arXiv 2510.15939",
        "size": 72,
        "type": "hallucination_detection",
    },
    "casp16_monomers": {
        "name": "CASP16 Monomer Targets",
        "description": (
            "Free modeling targets from CASP16 with experimental structures. "
            "CASP16 monomer assessment published (PMC 12157625, 2026). "
            "Single-domain fold prediction largely solved — no target folds "
            "incorrectly predicted. AF3 key for confidence estimation."
        ),
        "source": "predictioncenter.org",
        "reference": "PMC 12157625",
        "size": None,  # varies
        "type": "structure_accuracy",
    },
    "casp16_multimers": {
        "name": "CASP16 Multimer Targets",
        "description": (
            "Multimer/complex targets from CASP16 with experimental structures. "
            "CASP16 monomer assessment published (PMC 12157625, 2026). "
            "Multimer prediction remains challenging: <25% high quality. "
            "AF3 key for confidence estimation."
        ),
        "source": "predictioncenter.org",
        "reference": "PMC 12157625",
        "size": None,  # varies
        "type": "structure_accuracy",
    },
    "custom": {
        "name": "Custom Dataset",
        "description": "User-provided PDB pairs for benchmarking.",
        "source": "user",
        "size": None,
        "type": "structure_accuracy",
    },
}


def list_benchmark_datasets() -> dict:
    """List available benchmark datasets."""
    return {"datasets": BENCHMARK_DATASETS}


def evaluate_structure_pair(
    predicted_pdb: str,
    reference_pdb: str,
    target_name: str = "target",
) -> dict:
    """Evaluate a single predicted structure against a reference (experimental).

    Computes:
    - Global CA-RMSD
    - Per-residue CA distance distribution
    - GDT-TS (Global Distance Test - Total Score)
    - pLDDT-accuracy correlation
    - Residue-level accuracy breakdown
    """
    pred_coords = parse_ca_coords(predicted_pdb)
    ref_coords = parse_ca_coords(reference_pdb)
    pred_plddt = parse_plddt_from_pdb(predicted_pdb)

    if not pred_coords or not ref_coords:
        return {"error": "Could not parse CA atoms from one or both PDB inputs.", "target": target_name}

    ca_distances = compute_ca_distances(pred_coords, ref_coords)
    if not ca_distances:
        return {"error": "No overlapping residues between predicted and reference.", "target": target_name}

    rmsd = compute_rmsd(ca_distances)
    distances = np.array(list(ca_distances.values()))

    # GDT-TS: fraction of CA atoms within 1, 2, 4, 8 Angstroms
    gdt_1 = float(np.mean(distances < 1.0))
    gdt_2 = float(np.mean(distances < 2.0))
    gdt_4 = float(np.mean(distances < 4.0))
    gdt_8 = float(np.mean(distances < 8.0))
    gdt_ts = (gdt_1 + gdt_2 + gdt_4 + gdt_8) / 4.0

    # pLDDT vs accuracy correlation
    plddt_accuracy_corr = None
    if pred_plddt:
        matched_plddt = []
        matched_dist = []
        for res_idx, dist in ca_distances.items():
            if res_idx in pred_plddt:
                matched_plddt.append(pred_plddt[res_idx])
                matched_dist.append(dist)

        if len(matched_plddt) >= 3:
            # Higher pLDDT should correlate with lower distance (better accuracy)
            plddt_arr = np.array(matched_plddt)
            dist_arr = np.array(matched_dist)
            plddt_accuracy_corr = float(np.corrcoef(plddt_arr, -dist_arr)[0, 1])

    # Accuracy breakdown by pLDDT bucket
    bucket_accuracy = _accuracy_by_plddt_bucket(ca_distances, pred_plddt)

    return {
        "target": target_name,
        "aligned_residues": len(ca_distances),
        "ca_rmsd": round(rmsd, 2),
        "gdt_ts": round(gdt_ts, 4),
        "gdt_1A": round(gdt_1, 4),
        "gdt_2A": round(gdt_2, 4),
        "gdt_4A": round(gdt_4, 4),
        "gdt_8A": round(gdt_8, 4),
        "mean_ca_distance": round(float(np.mean(distances)), 2),
        "median_ca_distance": round(float(np.median(distances)), 2),
        "max_ca_distance": round(float(np.max(distances)), 2),
        "plddt_accuracy_correlation": (
            round(plddt_accuracy_corr, 3) if plddt_accuracy_corr is not None else None
        ),
        "accuracy_by_plddt_bucket": bucket_accuracy,
    }


def evaluate_batch(
    pairs: list[dict],
) -> dict:
    """Evaluate a batch of predicted vs reference structure pairs.

    Args:
        pairs: List of dicts with keys:
            - predicted_pdb: str (PDB content)
            - reference_pdb: str (PDB content)
            - target_name: str (optional)
    """
    start = time.time()
    results = []
    for i, pair in enumerate(pairs):
        result = evaluate_structure_pair(
            pair["predicted_pdb"],
            pair["reference_pdb"],
            target_name=pair.get("target_name", f"target_{i+1}"),
        )
        results.append(result)

    # Aggregate statistics
    valid = [r for r in results if "error" not in r]
    if valid:
        rmsds = [r["ca_rmsd"] for r in valid]
        gdt_scores = [r["gdt_ts"] for r in valid]
        summary = {
            "total_targets": len(pairs),
            "successful": len(valid),
            "failed": len(pairs) - len(valid),
            "mean_rmsd": round(float(np.mean(rmsds)), 2),
            "median_rmsd": round(float(np.median(rmsds)), 2),
            "mean_gdt_ts": round(float(np.mean(gdt_scores)), 4),
            "median_gdt_ts": round(float(np.median(gdt_scores)), 4),
            "best_target": valid[int(np.argmin(rmsds))]["target"],
            "worst_target": valid[int(np.argmax(rmsds))]["target"],
        }
    else:
        summary = {
            "total_targets": len(pairs),
            "successful": 0,
            "failed": len(pairs),
        }

    return {
        "summary": summary,
        "per_target": results,
        "elapsed_seconds": round(time.time() - start, 2),
    }


def _accuracy_by_plddt_bucket(
    ca_distances: dict[int, float],
    plddt: dict[int, float],
) -> dict:
    """Break down prediction accuracy by pLDDT confidence bucket."""
    buckets: dict[str, list[float]] = {
        "very_high_gt90": [],
        "high_70_90": [],
        "low_50_70": [],
        "very_low_lt50": [],
    }

    for res_idx, dist in ca_distances.items():
        score = plddt.get(res_idx)
        if score is None:
            continue
        if score > 90:
            buckets["very_high_gt90"].append(dist)
        elif score > 70:
            buckets["high_70_90"].append(dist)
        elif score > 50:
            buckets["low_50_70"].append(dist)
        else:
            buckets["very_low_lt50"].append(dist)

    result = {}
    for bucket_name, distances in buckets.items():
        if distances:
            arr = np.array(distances)
            result[bucket_name] = {
                "count": len(distances),
                "mean_ca_distance": round(float(np.mean(arr)), 2),
                "fraction_under_2A": round(float(np.mean(arr < 2.0)), 3),
            }
        else:
            result[bucket_name] = {"count": 0}

    return result


def generate_benchmark_report(
    batch_results: dict,
    dataset_name: str = "custom",
    backend_name: str = "unknown",
) -> dict:
    """Generate a formatted benchmark report suitable for publication.

    Returns structured data that can be included in a JOSS paper,
    bioRxiv preprint, or Zenodo dataset.
    """
    summary = batch_results.get("summary", {})
    per_target = batch_results.get("per_target", [])

    # Leaderboard entry format
    leaderboard_entry = {
        "backend": backend_name,
        "dataset": dataset_name,
        "n_targets": summary.get("total_targets", 0),
        "mean_rmsd": summary.get("mean_rmsd"),
        "median_rmsd": summary.get("median_rmsd"),
        "mean_gdt_ts": summary.get("mean_gdt_ts"),
        "median_gdt_ts": summary.get("median_gdt_ts"),
        "timestamp": time.time(),
    }

    # Per-bucket accuracy across all targets
    all_bucket_data: dict[str, list[float]] = {}
    for result in per_target:
        if "error" in result:
            continue
        for bucket, data in result.get("accuracy_by_plddt_bucket", {}).items():
            if data.get("mean_ca_distance") is not None:
                all_bucket_data.setdefault(bucket, []).append(data["mean_ca_distance"])

    aggregated_buckets = {}
    for bucket, values in all_bucket_data.items():
        arr = np.array(values)
        aggregated_buckets[bucket] = {
            "mean_ca_distance_across_targets": round(float(np.mean(arr)), 2),
            "n_targets_with_data": len(values),
        }

    return {
        "report_type": "benchmark",
        "dataset": BENCHMARK_DATASETS.get(dataset_name, {"name": dataset_name}),
        "leaderboard_entry": leaderboard_entry,
        "plddt_calibration": aggregated_buckets,
        "full_results": batch_results,
        "citation": (
            "If you use this benchmark in your research, please cite: "
            "FoldCopilot (https://github.com/adarsh/FoldCopilot)"
        ),
    }

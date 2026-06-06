"""Notebook export — emit Jupyter/Colab notebooks that reproduce analyses.

Generates a self-contained notebook from a FoldCopilot analysis result,
enabling end-to-end reproducibility and citation multiplier.
"""

from __future__ import annotations

import json
import time


def _make_cell(cell_type: str, source: str | list[str]) -> dict:
    """Create a Jupyter notebook cell."""
    if isinstance(source, str):
        source = source.splitlines(keepends=True)
    return {
        "cell_type": cell_type,
        "metadata": {},
        "source": source,
        **({"outputs": [], "execution_count": None} if cell_type == "code" else {}),
    }


def export_confidence_notebook(
    uniprot_id: str,
    confidence_report: dict,
) -> dict:
    """Export a confidence assessment as a reproducible Jupyter notebook.

    Returns notebook content as a dict (nbformat v4) and suggested filename.
    """
    cells = [
        _make_cell("markdown", [
            f"# FoldCopilot Confidence Assessment: {uniprot_id}\n",
            "\n",
            f"Generated: {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}\n",
            "\n",
            "This notebook reproduces the FoldCopilot confidence analysis.\n",
            "Run all cells to regenerate results.\n",
        ]),
        _make_cell("code", [
            "# Install dependencies\n",
            "# !pip install foldcopilot httpx numpy matplotlib\n",
        ]),
        _make_cell("code", [
            "import json\n",
            "import httpx\n",
            "import numpy as np\n",
            "\n",
            f'UNIPROT_ID = "{uniprot_id}"\n',
        ]),
        _make_cell("markdown", "## 1. Fetch pLDDT scores from AlphaFold DB\n"),
        _make_cell("code", [
            'url = f"https://alphafold.ebi.ac.uk/api/prediction/{UNIPROT_ID}"\n',
            "resp = httpx.get(url)\n",
            "resp.raise_for_status()\n",
            "metadata = resp.json()[0]\n",
            "print(f\"Model: {metadata.get(\\\"modelCreatedDate\\\", \\\"unknown\\\")}\")\n",
            "print(f\"Gene: {metadata.get(\\\"gene\\\", \\\"unknown\\\")}\")\n",
        ]),
        _make_cell("code", [
            "# Fetch PDB and extract pLDDT from B-factors\n",
            "pdb_url = metadata['pdbUrl']\n",
            "pdb_text = httpx.get(pdb_url).text\n",
            "\n",
            "plddt_scores = []\n",
            "for line in pdb_text.splitlines():\n",
            "    if line.startswith('ATOM') and line[12:16].strip() == 'CA':\n",
            "        plddt_scores.append(float(line[60:66]))\n",
            "\n",
            "print(f'Residues: {len(plddt_scores)}')\n",
            "print(f'Mean pLDDT: {np.mean(plddt_scores):.1f}')\n",
        ]),
        _make_cell("markdown", "## 2. pLDDT Distribution\n"),
        _make_cell("code", [
            "import matplotlib.pyplot as plt\n",
            "\n",
            "fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))\n",
            "\n",
            "# Per-residue plot\n",
            "ax1.plot(range(1, len(plddt_scores)+1), plddt_scores, linewidth=0.5)\n",
            "ax1.axhline(y=90, color='green', linestyle='--', alpha=0.5, label='>90: very high')\n",
            "ax1.axhline(y=70, color='orange', linestyle='--', alpha=0.5, label='>70: high')\n",
            "ax1.axhline(y=50, color='red', linestyle='--', alpha=0.5, label='>50: low')\n",
            "ax1.set_xlabel('Residue')\n",
            "ax1.set_ylabel('pLDDT')\n",
            "ax1.set_title(f'Per-residue pLDDT — {UNIPROT_ID}')\n",
            "ax1.legend()\n",
            "\n",
            "# Histogram\n",
            "ax2.hist(plddt_scores, bins=50, edgecolor='black', alpha=0.7)\n",
            "ax2.set_xlabel('pLDDT')\n",
            "ax2.set_ylabel('Count')\n",
            "ax2.set_title('pLDDT Distribution')\n",
            "plt.tight_layout()\n",
            "plt.show()\n",
        ]),
        _make_cell("markdown", "## 3. Original FoldCopilot Report\n"),
        _make_cell("code", [
            "# FoldCopilot confidence report (cached from original analysis)\n",
            f"report = {json.dumps(confidence_report, indent=2)}\n",
            "print(json.dumps(report, indent=2))\n",
        ]),
        _make_cell("markdown", [
            "## Citation\n",
            "\n",
            "If you use this analysis, please cite:\n",
            "- FoldCopilot (https://github.com/adarsh/FoldCopilot)\n",
            "- Jumper et al., Nature 596, 583–589 (2021) — AlphaFold\n",
        ]),
    ]

    notebook = {
        "nbformat": 4,
        "nbformat_minor": 5,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {"name": "python", "version": "3.11.0"},
            "foldcopilot": {
                "version": "0.1.0",
                "uniprot_id": uniprot_id,
                "generated_utc": time.time(),
            },
        },
        "cells": cells,
    }

    return {
        "notebook": notebook,
        "filename": f"foldcopilot_{uniprot_id}_confidence.ipynb",
        "format": "nbformat_v4",
    }


def export_benchmark_notebook(
    batch_results: dict,
    dataset_name: str = "custom",
    backend_name: str = "unknown",
) -> dict:
    """Export benchmark results as a reproducible notebook."""
    cells = [
        _make_cell("markdown", [
            f"# FoldCopilot Benchmark Report: {dataset_name}\n",
            "\n",
            f"Backend: {backend_name}\n",
            f"Generated: {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}\n",
        ]),
        _make_cell("code", [
            "import json\n",
            "import numpy as np\n",
            "import matplotlib.pyplot as plt\n",
            "\n",
            f"results = {json.dumps(batch_results, indent=2)}\n",
        ]),
        _make_cell("markdown", "## Summary Statistics\n"),
        _make_cell("code", [
            "summary = results['summary']\n",
            "for k, v in summary.items():\n",
            "    print(f'{k}: {v}')\n",
        ]),
        _make_cell("markdown", "## Per-Target Results\n"),
        _make_cell("code", [
            "per_target = results.get('per_target', [])\n",
            "valid = [t for t in per_target if 'error' not in t]\n",
            "if valid:\n",
            "    rmsds = [t['ca_rmsd'] for t in valid]\n",
            "    gdts = [t['gdt_ts'] for t in valid]\n",
            "    names = [t['target'] for t in valid]\n",
            "\n",
            "    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))\n",
            "    ax1.bar(names, rmsds)\n",
            "    ax1.set_ylabel('CA-RMSD (Å)')\n",
            "    ax1.set_title('Per-target RMSD')\n",
            "    ax1.tick_params(axis='x', rotation=45)\n",
            "\n",
            "    ax2.bar(names, gdts)\n",
            "    ax2.set_ylabel('GDT-TS')\n",
            "    ax2.set_title('Per-target GDT-TS')\n",
            "    ax2.tick_params(axis='x', rotation=45)\n",
            "    plt.tight_layout()\n",
            "    plt.show()\n",
        ]),
    ]

    notebook = {
        "nbformat": 4,
        "nbformat_minor": 5,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {"name": "python", "version": "3.11.0"},
        },
        "cells": cells,
    }

    return {
        "notebook": notebook,
        "filename": f"foldcopilot_benchmark_{dataset_name}_{backend_name}.ipynb",
        "format": "nbformat_v4",
    }

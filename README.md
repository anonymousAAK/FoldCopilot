<p align="center">
  <img src="https://img.shields.io/badge/MCP-Native-5A67D8?style=for-the-badge" alt="MCP Native" />
  <img src="https://img.shields.io/badge/License-MIT-22C55E?style=for-the-badge" alt="MIT License" />
  <img src="https://img.shields.io/badge/Python-3.11+-3776AB?style=for-the-badge&logo=python&logoColor=white" alt="Python 3.11+" />
  <img src="https://img.shields.io/badge/Backends-Boltz--2_|_OpenFold3_|_Chai--1-F59E0B?style=for-the-badge" alt="Backends" />
</p>

<h1 align="center">FoldCopilot</h1>

<p align="center">
  <strong>The confidence-aware interpretation copilot for protein structure predictions.</strong><br/>
  Not another wrapper. The interpretation layer the field is missing.
</p>

<p align="center">
  <a href="#quickstart">Quickstart</a> &middot;
  <a href="#why-foldcopilot">Why FoldCopilot</a> &middot;
  <a href="#tools">Tools</a> &middot;
  <a href="#architecture">Architecture</a> &middot;
  <a href="#benchmarks">Benchmarks</a> &middot;
  <a href="#citation">Citation</a>
</p>

---

## The Problem

AlphaFold changed biology. But it also introduced a new failure mode: **confidently wrong predictions that researchers trust without question.**

- **22% of intrinsically disordered residues** are falsely predicted as ordered by AlphaFold 3 ([arXiv 2510.15939](https://arxiv.org/abs/2510.15939))
- **No existing tool** cross-checks predictions against DisProt/MobiDB ground truth
- **No existing tool** compares outputs across Boltz-2, OpenFold3, and Chai-1 to surface disagreement
- **No existing MCP server** wraps Foldseek for structural similarity search
- Researchers copy-paste pLDDT scores without understanding what they mean

FoldCopilot fixes this. It sits between prediction backends and the researcher, adding the interpretation layer that turns raw predictions into trustworthy structural insights.

---

## Why FoldCopilot

| What exists today | What FoldCopilot adds |
|---|---|
| AFDB lookup servers return raw data | **Confidence interpretation** with pLDDT bucketing, PAE analysis, hallucination detection |
| Predictions come with no context | **IDR cross-checking** against DisProt 2026 + MobiDB ground truth |
| Each backend is a silo | **Cross-model disagreement detection** — when Boltz-2 and OpenFold3 disagree, you need to know |
| No Foldseek MCP exists | **First MCP wrapper for Foldseek** with confidence-filtered structural search |
| Wrappers host inference (expensive) | **Client-only architecture** — your GPU, your cloud, $0 from us |
| Commercial use is a license minefield | **Automatic license routing** — commercial queries never touch non-commercial weights |

<details>
<summary><strong>Competitive landscape (May 2026)</strong></summary>

| Project | Stars | What it does | What it doesn't do |
|---|---|---|---|
| ProteinMCP | ~4 | AF2-era protein engineering, 38 tools | No AF3, no confidence interpretation, no ensembling |
| AlphaFold-MCP-Server | 33 | AFDB REST lookup | Cannot run predictions, no interpretation |
| ChatMol/molecule-mcp | ~85 | PyMOL/ChimeraX visualization | No prediction, no confidence |
| BioinfoMCP | 38 tools | Classical NGS pipelines | Zero structure predictors |
| **FoldCopilot** | **New** | **Confidence interpretation + ensemble disagreement + 3 backends + Foldseek** | **This is the gap** |

</details>

---

## Quickstart

### Install

```bash
pip install foldcopilot
```

### Run as MCP Server

```bash
# Stdio transport (Claude Desktop, Cursor, etc.)
foldcopilot

# Or with Python
python -m foldcopilot.server
```

### Claude Desktop Configuration

Add to your `claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "foldcopilot": {
      "command": "foldcopilot",
      "env": {}
    }
  }
}
```

### First Query

Once connected, ask Claude:

> "Assess the confidence of AlphaFold's prediction for human p53 (P04637). Flag any regions where AlphaFold might be hallucinating structure in disordered regions."

FoldCopilot will:
1. Fetch the AFDB structure and per-residue pLDDT scores
2. Analyze the PAE matrix for domain positioning reliability
3. Cross-check against DisProt and MobiDB for known disordered regions
4. Flag hallucination warnings where AF predicts order in known IDRs
5. Return a structured `ConfidenceReport` with severity ratings and caveats

---

## Tools

### Confidence Interpretation (Moat 1)

| Tool | Description |
|---|---|
| `assess_confidence` | **Core tool.** pLDDT + PAE + DisProt/MobiDB IDR cross-check + hallucination warnings. Returns a structured `ConfidenceReport`. |
| `lookup_structure` | AFDB metadata by UniProt accession |
| `get_plddt_scores` | Per-residue confidence scores with bucket distribution |
| `get_pae_summary` | Predicted Aligned Error matrix summary |

### Structural Search (First Foldseek MCP)

| Tool | Description |
|---|---|
| `foldseek_search` | Search Foldseek with raw PDB content across AFDB, PDB, and more |
| `foldseek_search_uniprot` | Fetch AFDB structure by UniProt ID, then search Foldseek |
| `find_confident_homologs` | **AFDB-Foldseek loop**: find structural homologs where *both* query and target have confident structures |

### Structure Prediction

| Tool | Description |
|---|---|
| `predict_structure` | Run predictions via Boltz-2 (MIT), OpenFold3 (Apache-2.0), or Chai-1 (Apache-2.0). BYO-compute. |
| `check_backend_status` | Verify backend installation and GPU availability |
| `list_prediction_backends` | All backends with license type and installation status |

### Ensemble Comparison (Moat 2)

| Tool | Description |
|---|---|
| `compare_predictions` | **Cross-model disagreement detection.** Feed two PDB outputs, get per-residue agreement classification, disagreement spans, RMSD, and pLDDT correlation. |
| `compare_prediction_files` | Same as above, from file paths |

---

## Architecture

```
                          +------------------+
                          |   Claude / LLM   |
                          +--------+---------+
                                   |
                              MCP Protocol
                                   |
                          +--------v---------+
                          |   FoldCopilot    |
                          |   MCP Server     |
                          +--------+---------+
                                   |
                 +-----------------+-----------------+
                 |                 |                 |
        +--------v------+ +------v-------+ +-------v--------+
        | Confidence    | | Foldseek     | | Prediction     |
        | Interpreter   | | Search       | | Engine         |
        +--------+------+ +------+-------+ +-------+--------+
                 |                |                  |
        +--------v------+ +------v-------+ +-------v--------+
        | AFDB + DisProt| | Foldseek    | | Boltz-2        |
        | + MobiDB      | | Web API     | | OpenFold3      |
        +---------------+ +-------------+ | Chai-1         |
                                           +----------------+
                                                  |
                                           +------v-------+
                                           | Ensemble     |
                                           | Comparator   |
                                           +--------------+
```

### Design Principles

- **Interpretation over wrapping.** The value is in the confidence layer, not the API calls.
- **Client-only.** We never host inference. Your GPU, your cloud, your cost.
- **License-aware.** Commercial queries are automatically routed to MIT/Apache-2.0 backends. AF3 weights require explicit non-commercial attestation.
- **Reproducible.** Every prediction ships with a `ReproducibilityManifest`: model version, weights hash, input hash, parameters, runtime environment, GPU type, timestamp.
- **Compact responses.** Raw PDB/CIF files are saved locally. MCP responses contain summaries, paths, and URIs — never multi-megabyte payloads.
- **Content-addressed caching.** Predictions are deterministic and expensive. Cache key: `(sequence_hash, model_version, params_hash)`. Every cache hit is pure win.

---

## Prediction Backends

| Backend | License | Speed | Affinity | Status |
|---|---|---|---|---|
| **Boltz-2** | MIT | ~20s/GPU | Yes (Pearson r=0.66) | Default |
| **OpenFold3** | Apache-2.0 | ~minutes | No | Commercial-safe AF3 |
| **Chai-1** | Apache-2.0 | ~minutes | No | Multi-chain + ligands |

### License Routing

```
commercial_use=True  --> Boltz-2 (MIT) or OpenFold3 (Apache-2.0) or Chai-1 (Apache-2.0)
commercial_use=False --> All backends available
AF3 weights          --> NEVER auto-selected. BYO-weights + non-commercial attestation only.
Chai-2               --> NOT supported. Closed API, ToS prohibits relay.
```

### BYO Compute

FoldCopilot is a client. Bring your own GPU:

```bash
# Local GPU
pip install boltz

# Cloud GPU (Modal)
modal run deploy/modal_boltz2.py

# Cloud GPU (RunPod)
runpodctl start --gpu A100 --image boltz2:latest
```

---

## Hallucination Detection

AlphaFold 3 predicts ordered structure for ~22% of residues that are actually intrinsically disordered ([arXiv 2510.15939](https://arxiv.org/abs/2510.15939)). These hallucinations can mislead drug design and functional annotation.

FoldCopilot's `assess_confidence` tool cross-references every prediction against:

- **DisProt 2026** — curated ground truth for intrinsically disordered regions
- **MobiDB** — aggregated disorder predictions from multiple sources

When AlphaFold reports pLDDT > 70 for a residue that DisProt says is disordered, FoldCopilot raises a `HallucinationWarning` with severity rating:

| Severity | Condition | Action |
|---|---|---|
| **High** | AF pLDDT > 70 in known IDR | Do not trust this region. Validate experimentally. |
| **Moderate** | AF pLDDT 50-70 in known IDR | Treat with caution. Likely disordered. |
| None | AF pLDDT < 50 in known IDR | AF agrees with IDR databases. Region is disordered. |

---

## Cross-Model Disagreement

When two models agree, trust goes up. When they disagree, you need to know.

```
compare_predictions(boltz2_pdb, openfold3_pdb)
```

Returns per-residue classification:

| Agreement Level | Meaning | What to do |
|---|---|---|
| `strong_agree` | Both confident, structures match (<3A) | High structural confidence |
| `moderate_agree` | One confident, structures match | Reasonable confidence |
| `disagree` | Both confident, structures DIFFER (>3A) | **Validate experimentally** |
| `both_uncertain` | Neither model confident | Likely disordered or flexible |

Contiguous disagreement spans are flagged with interpretations and mean CA distances.

---

## Benchmarks

### Confidence Assessment

FoldCopilot's hallucination detection was developed against the dataset from [arXiv 2510.15939](https://arxiv.org/abs/2510.15939):

- 72 DisProt proteins with curated disorder annotations
- 22% of IDR residues hallucinated as ordered by AF3
- 18% of biological-process residues hallucinated

### Test Coverage

```
92 tests | 5 test modules | All passing
```

```
tests/test_confidence.py    15 tests  (pLDDT bucketing, span detection, hallucination)
tests/test_foldseek.py      13 tests  (alignment parsing, UniProt extraction)
tests/test_predict.py       17 tests  (license routing, manifests, I/O parsing)
tests/test_ensemble.py      25 tests  (RMSD, agreement classification, span detection)
tests/test_validation.py    22 tests  (sequence, UniProt, PDB input validation)
```

---

## Roadmap

- [x] **v0.1** Confidence Copilot MVP — AFDB + pLDDT + PAE + DisProt + MobiDB + hallucination detection
- [x] **v0.2** Foldseek MCP — first MCP wrapper for Foldseek + AFDB-Foldseek agent loop
- [x] **v0.3** Boltz-2 backend — MIT, 20s/GPU, affinity prediction, reproducibility manifests
- [x] **v0.4** Ensemble disagreement — cross-model comparison, per-residue agreement, span detection
- [x] **v0.5** OpenFold3 + Chai-1 backends, input validation, JOSS readiness
- [ ] **v0.6** AlphaMissense + AlphaFill cofactor transplantation
- [ ] **v0.7** Therapeutic vertical packs (Antibody, Kinase, GPCR)
- [ ] **v0.8** Education mode (`--explain` flag for plain-language confidence reading)
- [ ] **v0.9** Benchmarking harness (CASP16, Polaris-ASAP, DisProt hallucination set)
- [ ] **v1.0** JOSS submission + Zenodo dataset DOI + public leaderboard

---

## Research Use

> **FoldCopilot is for research use only.** Every `ConfidenceReport` and `EnsembleReport` includes standard caveats. Do not use for clinical decisions.

Key caveats embedded in every response:

1. pLDDT > 70 does not guarantee correctness — it indicates the model's own confidence.
2. AlphaFold can hallucinate ordered structure in intrinsically disordered regions.
3. High PAE (>10A) at interfaces suggests unreliable domain/chain positioning.
4. Cross-model agreement increases confidence but does not eliminate shared biases.

---

## Citation

If FoldCopilot is useful in your research, please cite:

```bibtex
@software{foldcopilot2026,
  title     = {FoldCopilot: Confidence-Aware Interpretation Copilot for Protein Structure Predictions},
  author    = {Adarsh},
  year      = {2026},
  url       = {https://github.com/adarsh/FoldCopilot},
  license   = {MIT}
}
```

See also: [`CITATION.cff`](CITATION.cff) for machine-readable citation metadata.

---

## Contributing

FoldCopilot is MIT-licensed and welcomes contributions. Areas of high impact:

- **New backends** — Protenix, AQAffinity, ESMFold
- **Database integrations** — AlphaMissense, AlphaFill, CATH, SCOP, SAbDab, KLIFS
- **Therapeutic verticals** — Antibody, Kinase, GPCR domain packs
- **Benchmarks** — CASP16, Polaris-ASAP evaluation sets
- **Wet-lab validation** — partner with us to validate predictions experimentally

```bash
# Development setup
git clone https://github.com/adarsh/FoldCopilot.git
cd FoldCopilot
pip install -e ".[dev]"
pytest
```

---

## Governance

This project includes a continuity plan. If the primary maintainer becomes unavailable:

1. All code is MIT-licensed and fully open
2. No hosted infrastructure to maintain — client-only architecture
3. All external API dependencies are public and documented
4. Test suite is comprehensive and self-contained
5. Reach out via GitHub Issues to volunteer as co-maintainer

---

<p align="center">
  <sub>Built for the researchers who need to know when to trust a fold and when to reach for the pipette.</sub>
</p>

---
title: 'FoldCopilot: A Confidence-Aware Interpretation Copilot for Protein Structure Predictions via the Model Context Protocol'
tags:
  - Python
  - protein structure prediction
  - AlphaFold
  - confidence interpretation
  - Model Context Protocol
  - bioinformatics
  - structural biology
authors:
  - name: Adarsh
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: Independent Researcher
    index: 1
date: June 2026
bibliography: paper.bib
---

# Summary

FoldCopilot is an open-source Model Context Protocol (MCP) server that provides confidence-aware interpretation of protein structure predictions. It integrates six prediction backends --- Boltz-2, OpenFold3, Chai-1, Protenix-v2, AlphaFold 3 (BYO-weights), and AQAffinity --- with structured confidence analysis, intrinsically disordered region (IDR) cross-validation against DisProt and MobiDB, cross-model ensemble disagreement detection, and therapeutic-area vertical analysis packs for antibodies, kinases, and GPCRs. Implemented in Python atop FastMCP 3.x, FoldCopilot exposes 33 MCP tools that enable large language model (LLM) agents to autonomously assess, compare, and interpret protein structures with scientific rigor. Every prediction is accompanied by a reproducibility manifest, and every confidence assessment embeds uncertainty caveats to prevent over-interpretation of computational models.

# Statement of Need

The release of AlphaFold 2 [@jumper2021highly] and subsequently AlphaFold 3 [@abramson2024accurate] transformed structural biology by making high-accuracy protein structure prediction widely accessible. However, these advances introduced a subtle but consequential failure mode: confidently wrong predictions, particularly in intrinsically disordered regions (IDRs). A systematic evaluation by Alderson et al. [@alderson2025idr] found that AlphaFold 3 falsely predicts ordered structure for 22% of residues annotated as disordered in DisProt [@quaglia2026disprot], with 18% of biological-process-associated residues similarly hallucinated. These hallucinations can mislead downstream analyses including drug target assessment, functional annotation, and protein engineering campaigns. Complementary work has shown that AlphaFold 2 is preferable to AlphaFold 3 for disorder detection because it avoids structural hallucinations in disordered regions [@af2_idr_preferred], suggesting a dual-model strategy where AF2 flags disorder and AF3 or Boltz-2 provides structural predictions.

Despite the critical importance of confidence interpretation, no existing tool in the MCP ecosystem combines structure prediction with systematic confidence assessment and cross-model ensemble analysis in an agentic interface. The MCP ecosystem has grown rapidly --- Smithery alone indexes over 7,000 servers --- yet none provides confidence-aware structure prediction interpretation. Wet-lab biologists frequently encounter raw pLDDT scores and PAE matrices without the domain expertise to interpret them correctly. A researcher seeing pLDDT > 70 for a given residue may assume structural reliability without checking whether that residue falls in a known IDR, where the prediction is likely a hallucination. FoldCopilot addresses this gap by providing structured, citation-backed confidence reports that translate numerical scores into actionable scientific guidance, complete with severity ratings and explicit caveats.

Furthermore, the protein structure prediction landscape has fragmented across multiple backends with different licenses, accuracy profiles, and capabilities. Boltz-2 [@passaro2025boltz2] offers MIT-licensed predictions with affinity estimation; OpenFold3 [@openfold3] provides an Apache-2.0-licensed reproduction of AlphaFold 3; Protenix-v2 [@protenix2026] achieves state-of-the-art antibody-antigen prediction under Apache-2.0. When two independently trained models agree on a structural feature, confidence increases; when they disagree, the researcher needs to know. No prior tool surfaces this cross-model disagreement information in an agentic workflow.

# State of the Field

Several MCP servers address adjacent problems in computational biology, but none occupies the confidence-aware structure prediction interpretation niche. ProteinMCP [@xu2026proteinmcp] is the closest competitor: a peer-reviewed 38-tool framework for protein engineering that wraps AlphaFold 2-era tools including MSA generation, ESM embeddings, and BindCraft binder design. However, it lacks AlphaFold 3 or Boltz-2 backends, provides no confidence interpretation or IDR cross-checking, and does not support ensemble disagreement detection. AlphaFold-MCP-Server provides REST-based lookup against the AlphaFold Database but cannot run predictions and offers no interpretation beyond raw score retrieval. BioinfoMCP automates classical next-generation sequencing pipelines (FastQC, BWA, GATK) with zero structure predictors among its 38 tools. FoldRun MCP is a mock demonstration server for Gemini Enterprise that exposes submission and status endpoints but performs no real prediction or interpretation. BioMCP and its OncoMCP extension focus on clinical genomics --- ClinicalTrials.gov, PubMed, variant databases --- with no structural biology capabilities. ChatMol/molecule-mcp provides PyMOL and ChimeraX visualization via MCP but no prediction. OmniFold containerizes AlphaFold 3, Chai-1, and Boltz-2 but is not MCP-native and lacks any interpretation layer.

FoldCopilot is, to our knowledge, the first MCP server that combines multiple prediction backends with structured confidence interpretation, IDR hallucination detection, cross-model ensemble disagreement analysis, and domain-specific therapeutic vertical packs in a single agentic interface.

# Key Features

FoldCopilot's feature set is organized around confidence interpretation as the primary value layer, with prediction backends as configurable plug-ins:

- **Structured ConfidenceReport**: Every assessment returns per-residue pLDDT bucketing (very high / confident / low / very low), PAE entropy analysis for domain and interface reliability, and IDR cross-validation against DisProt 2026 and MobiDB ground truth.
- **Hallucination detection**: Residues where a predictor reports pLDDT > 70 but DisProt annotates as disordered are flagged with severity ratings (high, moderate, none) and actionable guidance.
- **Dual-model IDR strategy**: AF2 pLDDT < 50 is used as a hallucination-free disorder signal, while AF3 or Boltz-2 provides structural predictions, following evidence that AF2 avoids IDR hallucinations [@af2_idr_preferred].
- **Cross-model ensemble disagreement detection**: Per-residue agreement classification (strong agree, moderate agree, disagree, both uncertain) with contiguous disagreement span identification and mean C-alpha distances.
- **License-aware routing**: Commercial queries are automatically routed to MIT (Boltz-2) or Apache-2.0 (OpenFold3, Chai-1, Protenix-v2) backends. AlphaFold 3 requires explicit non-commercial attestation and user-supplied weights.
- **Reproducibility manifests**: Every prediction records model version, weights hash, input sequence hash, parameters, runtime environment, GPU type, and timestamp.
- **Therapeutic vertical packs**: Antibody analysis with Kabat CDR identification and CDR-H3 confidence warnings; kinase analysis with ATP-site annotation, DFG motif detection, and KLIFS cross-reference; GPCR analysis with transmembrane helix confidence and TMalphaFold membrane topology [@tmalphaFold].
- **Education mode**: Plain-language explanations of pLDDT, PAE, and hallucination concepts with analogies and citations, designed for wet-lab biologists without computational expertise.
- **Benchmarking harness**: Evaluation against CASP16 [@casp16] targets and DisProt hallucination datasets with CA-RMSD, GDT-TS, and pLDDT calibration metrics.

# Implementation

FoldCopilot is implemented in Python using FastMCP 3.x [@fastmcp] atop the Model Context Protocol specification [@mcp_spec]. Data models are defined with Pydantic for runtime validation and serialization. Long-running structure predictions use the `@mcp.tool(task=True)` decorator implementing the MCP background-task protocol, allowing clients to receive a task identifier immediately and poll for progress updates (e.g., "MSA complete", "diffusion step 12/200", "ranking models"). A content-addressed cache keyed on `(sequence_hash, model_version, params_hash)` avoids redundant GPU-intensive computations. Raw PDB/CIF files are stored locally; MCP responses contain compact JSON summaries with file paths, keeping payloads under 50 KB. OpenTelemetry instrumentation provides per-tool-call traces and timing metrics. The server operates as a client-only architecture: it orchestrates prediction backends installed on user-supplied compute (local GPU, Modal, RunPod, or NVIDIA NIM) and hosts no inference itself. The test suite comprises 205 tests across 15 modules covering confidence analysis, ensemble comparison, prediction routing, input validation, annotations, therapeutic verticals, education mode, benchmarking, notebook export, fold-drift detection, and observability.

# Acknowledgements

FoldCopilot builds upon the foundational work of the AlphaFold team at DeepMind, the Boltz team at MIT, the AlQuraishi Lab (OpenFold3), Chai Discovery, the ByteDance Research team (Protenix), and the Foldseek developers. The DisProt and MobiDB databases provide essential ground truth for disorder annotation. The FastMCP framework and the Anthropic MCP specification provide the protocol infrastructure. We thank the open-source structural biology community for making this work possible.

# References

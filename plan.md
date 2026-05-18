# Critical Review and Improvement Plan for FoldCopilot

## TL;DR
- **Build it, but reframe and resequence.** The multi-backend prediction copilot niche is genuinely unoccupied as of May 2026 — ProteinMCP (Xu et al., Protein Science 2026; ~4 GitHub stars) is AF2-era protein-engineering, Augmented-Nature/AlphaFold-MCP-Server (33 stars) is AFDB lookup-only, ChatMol/molecule-mcp (~85 stars) is visualization-only, and BioinfoMCP's 38 auto-converted tools contain zero structure predictors. But the **defensible** product is not "wrap six folders behind MCP" — it's **confidence-aware ensemble interpretation**, which the current plan buries in v0.4 at week 7+.
- **The riskiest assumption is sequencing the killer feature last.** Pull confidence/hallucination interpretation (pLDDT + PAE + IDR flagging against DisProt 2026 / MobiDB + cross-model disagreement) into v0.1. Cut Chai-2 and hosted-Tamarind/Neurosnap wrappers from the roadmap entirely (ToS landmines). Treat AlphaFold 3 as BYO-weights only — its CC-BY-NC-SA 4.0 source license and non-commercial weights ToU make a permissive wrapper that hosts AF3 inference for arbitrary users legally untenable.
- **Optimize distribution for citations, not karma.** JOSS + bioRxiv + a Zenodo-DOI'd "AF3-hallucinated DisProt set" + CZI EOSS Cycle 7 LOI form the moat. A coordinated 48-hour HN/Reddit/X launch (Tue–Thu, mid-morning PT) drives the spike — arXiv 2511.04453 finds AI-tool HN launches average 121 stars in 24 h, 189 in 48 h, 289 in week 1, and the "Show HN" tag shows no statistical advantage after controls — but a citable scientific artifact is what compounds.

## Key Findings

### 1. Competitive landscape — narrow but real gap
- **ProteinMCP** (github.com/charlesxu90/ProteinMCP, MIT, ~4 stars, ~1 contributor, 125 commits, no test suite; bioRxiv 10.64898/2026.03.11.711149, Protein Science 35(4):e70547): 38-tool agentic framework whose shipped catalogue is `msa_mcp, alphafold2_mcp, mmseqs2_mcp, esm_mcp, prottrans_mcp, plmc_mcp, ev_onehot_mcp, bindcraft_mcp, boltzgen_mcp`. **Does not wrap AF3, Boltz-2, Chai-1, Protenix, Foldseek, or AFDB.** Workflow skills cover fitness modeling, BindCraft binder design, BoltzGen nanobody design. **No confidence interpretation, no IDR flagging, no cross-model ensembling.** Frames the agentic-engineering niche, not the prediction-interpretation niche.
- **Augmented-Nature/AlphaFold-MCP-Server** (MIT, 33 stars, 6 commits, 1 contributor, released Jul 24 2025): pure REST wrapper over `alphafold.ebi.ac.uk/api/`. Exposes `get_confidence_scores` and `analyze_confidence_regions` for cached AFDB AF2 structures; **cannot run any prediction.** Sister PDB-MCP-Server is the same pattern. Both are thin and stale.
- **ChatMol/molecule-mcp** (~85 stars, Mar 2025): PyMOL + ChimeraX + GROMACS Copilot command execution. Visualization-only; no prediction.
- **BioinfoMCP** (arXiv 2510.02139, Oct 2025): auto-converter + 38 tools, all classical NGS/sequence pipelines (FastQC, BWA, samtools, GATK, MACS3, etc.). The widely-quoted "94.7% successfully executed complex workflows" is on RNA-seq/ChIP-seq/ATAC-seq/WGS — **no structure predictor is in the 38**, so the metric does not transfer.
- **No Foldseek MCP wrapper exists.** No DisProt/MobiDB MCP. No tool combines AFDB lookup + Foldseek search + prediction in one agent loop.
- **Closest non-MCP analogue: OmniFold** (PMC12330748), a containerized wrapper around AF3 + Chai-1 + Boltz-2; not agentic, no Protenix, no MCP, no interpretation layer.
- **Commercial SaaS** (Tamarind Bio, Neurosnap, 310.ai, Latch Bio, Rowan, Recursion-hosted Boltz): closed, not MCP-native, paid. They own the enterprise path. The developer/agentic path is open.

**Verdict**: "Multi-backend agentic structure-prediction copilot with confidence-aware interpretation" is unoccupied.

### 2. Model landscape — update the backend roadmap

| Model | Status May 2026 | License | Implication |
|---|---|---|---|
| **AlphaFold 3** | Code + gated weights since Nov 2024, fully public Feb 2025 | Source: CC-BY-NC-SA 4.0; weights: non-commercial, no redistribution, no training of derivative folders | **License landmine.** BYO-weights pattern only; cannot host for commercial users. |
| **Boltz-2** | Released Jun 2025 (Passaro et al., bioRxiv 2025.06.14.659707); 20 s/GPU, predicts affinity with a Pearson correlation of 0.66 on the 4-target FEP+ subset ("Boltz-2 achieves an average Pearson correlation of 0.66, outperforming all available inexpensive physical methods and ML baselines") | **MIT** | **Default MVP backend.** |
| **Chai-1** | Open Sep 2024 | **Apache-2.0** (commercial OK) | Add at v0.4. |
| **Chai-2** | Jun 2025; closed, partner-only via Responsible Deployment policy; 16% antibody hit rate across 52 targets | Closed | **Drop entirely** — wrapping a closed API violates ToS. |
| **Protenix** | ByteDance | Apache-2.0 | Defer; lower marginal value. |
| **OpenFold3 (preview)** | Released Oct 28 2025 (AlQuraishi Lab + LLNL + Steinegger Lab); bitwise reproduction of AF3 | **Apache-2.0** | **Add as MVP backend** — commercial-safe AF3 substitute. |
| **AQAffinity** | Jan 2026 (SandboxAQ on top of OpenFold3); structure-free affinity prediction | Open | Add as affinity option. |
| **ESM3 / ESM C** | ESM3-open small (1.4B) non-commercial; ESM C 300M/600M open-weight | Mixed | Embedding pre-screening only. |
| **ESMFold** | Single-sequence baseline | MIT | Keep as fast path. |
| **Foldseek** | Webserver REST API documented; **no MCP wrapper exists** | GPL-3.0 CLI / Apache app | **First-mover opportunity.** |
| **AFDB 2025** | Redesigned interface (Bertoni et al., NAR Nov 2025, gkaf1226); Foldseek integrated, 3D-Beacons linked, AlphaMissense + AlphaFill | Open | Use new endpoints. |
| **DisProt 2026 / MobiDB** | IDR ground truth | Open | First-class resources behind confidence layer. |

**Net changes vs the plan**: **add** OpenFold3, AQAffinity, DisProt, MobiDB, Foldseek MCP, AlphaMissense, AlphaFill; **defer** Protenix; **drop** Chai-2 and any hosted-SaaS wrapper; **gate** AF3 behind BYO-weights + non-commercial banner.

### 3. Technical architecture

- **Framework**: FastMCP 3.x (Python) over the official SDK. FastMCP 3.0 wraps Docket (originally built by Prefect, "processes millions of concurrent tasks every day") and gives `@mcp.tool(task=True)` decorators implementing MCP SEP-1686 (background-task protocol) — exactly the primitive for 2–5 minute folds. Clients receive a task ID immediately, then poll/stream for progress and result.
- **Transport**: Streamable HTTP per MCP spec 2025-11-25, with stdio fallback for Claude Desktop local. Google ADK 1.25.1 (March 2026) added native streaming progress events; do not roll your own polling.
- **Long-running jobs**: `await ctx.report_progress` + structured stage messages ("MSA done", "diffusion step 12/200", "ranking models", "physics relax"). Use Redis-backed Docket for persistence — in-memory storage drops jobs on restart and is unfit for production.
- **Large binaries**: never return raw PDB/CIF/PAE in tool responses. Stash to object store (S3/MinIO/local-disk presigned URLs); return URIs plus a compact JSON summary (mean pLDDT per chain, top-N low-confidence spans, RMSD vs reference, PAE matrix hash). Keep MCP responses < ~50 kB; let the model fetch detail via a `read_resource` follow-up.
- **Auth**: OAuth 2.1 per the Nov 2025 spec; encrypted env-var vault (HasMCP pattern) for backend API keys. Security baseline matters: per BlueRock Security 2026 (cited by Apigene/Security Boulevard), 36.7% of 7,000+ public MCP servers have SSRF vulnerabilities, 43% have unsafe command-execution paths, and an April 2026 audit found 41% of official-registry servers have zero auth — do not be in those buckets.
- **Caching**: deterministic on `(sequence_hash, model_version, params_hash)`. Predictions are expensive and reproducible — every cache hit is pure win. Content-addressed local disk + optional S3 backend. BioMCP's caching layer is a usable reference.
- **Observability**: OpenTelemetry traces per tool call, per-stage timings, `/health` and `/metrics`. Surface model versions, weights hash, and seed in every output.

### 4. Scope and sequencing critique
The plan (v0.1 AFDB read → v0.2 Foldseek+search → v0.3 prediction backends → v0.4 compare+interpret → v0.5 polish → v0.6+ post-launch) **buries the killer feature at week 7+ and front-loads commodity work** that Augmented-Nature already ships.

**Proposed resequencing**:
- **v0.1 (week 1–2) — "Confidence Copilot" MVP, zero GPU.** AFDB lookup + **confidence interpretation**: pLDDT bucketing, PAE matrix summarization, IDR flagging vs DisProt 2026 + MobiDB, AlphaFold hallucination warnings. The interpretation layer ships *before* the prediction layer. The arXiv 2510.15939 IDR paper provides a ready-made evaluation set: "32 percent of residues are misaligned with DisProt, with 22 percent representing hallucinations where AlphaFold3 incorrectly predicts order in disordered" regions, and 18% of biological-process residues hallucinated — that is your concrete, demonstrable, citable problem.
- **v0.2 (week 3–4) — Foldseek MCP + AFDB↔Foldseek agent loop.** Nobody has built a Foldseek MCP. Lowest effort, highest visibility, zero GPU cost. Use the public Foldseek webserver REST API.
- **v0.3 (week 5–6) — one runnable backend (Boltz-2 only).** MIT, fast (20 s/GPU), includes affinity, GPU-deployable via Modal/RunPod/NIM.
- **v0.4 (week 7–8) — ensemble + cross-model disagreement.** Second moat: when AF3 and Boltz-2 agree at a residue/interface, surface high confidence; when they disagree, flag it. No published MCP does this.
- **v0.5 (week 9–10) — second runnable backend (OpenFold3 or Chai-1) + polish + JOSS submission.**
- **v0.6+ — Protenix, AQAffinity, AlphaMissense, AlphaFill, vertical packs.**

**Cuts from MVP**: Chai-2 (closed), Tamarind/Neurosnap/310.ai wrappers (SaaS ToS risk), Protenix (low marginal value at MVP), any non-AFDB structure DB integration beyond Foldseek, fancy UI work. MVP is CLI + Claude Desktop + a 90-second demo video.

**Riskiest assumption**: that wrapping six backends is the value. It isn't — that's commodity work a BioinfoMCP-style auto-converter replicates in a weekend. The value is **interpretation + disagreement-surfacing + packaging**. Reframe the plan around the interpretation layer with backends as plug-ins.

### 5. Differentiation and moat
- **Versus a clone of the MCP scaffold**: durable assets are (a) **curated benchmarks/datasets** — a public, versioned "AF3-hallucinated DisProt set," "AF3-vs-Boltz2 disagreement maps for the human proteome," "AFDB-Foldseek-confirmed remote-homologue set"; (b) **prompts/skills** — well-tested Claude Skills/MCP prompts for "is this structure trustworthy?" workflows; (c) **citations** — once 5–10 papers cite FoldCopilot, switching costs accumulate.
- **Versus Anthropic/OpenAI/Google building first-party**: frontier labs build general MCP scaffolding, not domain-specific interpretation. Stay close to the science — integrate DisProt 2026, MobiDB, AlphaMissense, AlphaFill cofactor transplantation, TMalphaFold, CATH/SCOP, SAbDab, KLIFS, PDBe-KB. Frontier labs will not.
- **Community moats**: Discord/Slack, biweekly office hours, <24-hour PR review (AFFiNE playbook), an AFDB cross-link, partnership with an academic group for wet-lab validation. Get listed in the official MCP registry, Smithery (~6,000 distinct servers per Apigene's March 2026 Smithery CLI guide; Smithery's quality score directly determines search-result placement), Glama, PulseMCP, mcp.so (combined raw cross-registry total is ~20,000 per toolradar.com March 2026, "most are duplicates or weekend experiments"), mcpservers.org. Use mcp-submit to fan out.

### 6. Distribution and recognition
- **Launch mechanics**: per arXiv 2511.04453 (138 AI/LLM HN launches 2024–25), HN exposure yields **121 stars at 24 h, 189 at 48 h, 289 at one week**; the "Show HN" tag itself shows no statistical advantage after controlling for other factors — **timing and content matter more than the tag**. Tuesday–Thursday, mid-morning Pacific. Concentrate Reddit (r/bioinformatics, r/MachineLearning, r/mcp), HN, X/Twitter, LinkedIn, and the bioRxiv post into a single 48-hour window. Have 100–200 baseline stars from your network before launch — readers convert at ~5% from a zero-star README.
- **Paper venues**: **JOSS** is the right fit — free, diamond OA, ~1,000-word paper, peer reviewers actually test the software; "bioinformatics" is an established tag. Pair with a **bioRxiv** preprint that documents a *new* finding (e.g., AF3-vs-Boltz2 disagreement maps for the human proteome). The software-only JOSS paper alone won't drive citations. **Bioinformatics** "Application Notes" is a strong stretch target if you have wet-lab validation; **SoftwareX** has lower visibility.
- **Grants currently open or near-open (May 2026)**:
  - **Anthropic Fellows** — next cohort starts July 20 2026; applications closed April 26 for July, rolling for later cohorts. $3,850/wk + ~$15k/month compute; US/UK/Canada only. AI-safety focus — **not a natural fit** unless reframed as agent safety / scientific correctness.
  - **CZI EOSS** — Cycle 6 closed; Cycle 7 typically opens annually. EOSS supports *maintenance of already-widely-used tools* — apply after FoldCopilot has demonstrable adoption (target: ≥1,000 stars and ≥3 citing groups before LOI).
  - **Schmidt Sciences / Virtual Institute for Scientific Software (VISS)** — funds open-source scientific software at scale; lower adoption bar than CZI.
  - **Wellcome Trust** — co-funds EOSS and the Open Research Fund.
  - **Anthropic AI for Science** — compute credits; check anthropic.com for current call.
  - **India-specific**: ANRF (Anusandhan National Research Foundation) PMECRG; Wellcome/DBT India Alliance Early Career Fellowship; iHub-Anubhuti / IIT TIH grants. Most require institutional affiliation — partner with an IISc/IIT/NCBS PI if solo.
- **Listings to hit**: official MCP registry, Smithery (work to a 100/100 quality score — affects ranking), Glama, PulseMCP, mcp.so, mcpservers.org, awesome-mcp-servers, bioconda recipe, PyPI, NVIDIA NIM blueprint.

### 7. Risks and failure modes
1. **AlphaFold 3 license**: source CC-BY-NC-SA 4.0; weights non-commercial, no redistribution, no training-derivative-folders. Hosting AF3 inference for arbitrary users from an MIT/Apache-licensed wrapper is **not legal**. Mitigation: BYO-weights pattern, route commercial users to OpenFold3 (Apache-2.0) or Boltz-2 (MIT), keep AF3 access optional and gated behind explicit user attestation.
2. **Hosted-API ToS**: wrapping Chai-2, Tamarind, Neurosnap, or 310.ai's hosted APIs as backends almost certainly violates standard SaaS ToS (no automated relay, no reseller). Mitigation: Chai-1 is Apache-2.0 — run it yourself; do not relay Chai-2; do not wrap commercial SaaS without a written agreement.
3. **Scientific-correctness liability**: a copilot that tells a biologist a hallucinated active site is "trustworthy" is a real harm vector. AF3 hallucinations in disordered regions are well-documented (arXiv 2510.15939). Mitigation: hard-coded uncertainty disclaimers in every tool response; never return "trustworthy" as a binary; emit a structured `ConfidenceReport` with caveats; cite EBI's "What AlphaFold 3 struggles with" guidance and arXiv 2510.15939 in responses; add a `--research-use-only` license banner.
4. **Upstream agentic interfaces**: Anthropic could ship a first-party "AlphaFold Skill"; Google could expose AlphaFold Server via Gemini extensions. Mitigation: own the **interpretation + ensembling** layer with deep domain integration (DisProt 2026, MobiDB, AlphaMissense, CATH) — labs will not match this depth.
5. **Maintainer burnout**: solo dev, six backends, ten weeks is aggressive. Mitigation: cut scope as above, automate tests, ≤3 backends at v1.0, write a "if I disappear" governance doc.
6. **Reproducibility / "fold drift"**: backends update; same input gives different output across versions. Mitigation: pin model versions, embed model commit hash and weights checksum in every output manifest, ship `reproducibility_manifest.json` per run.
7. **Hosted-compute cost trap**: GPU bills scale faster than donations. Mitigation: ship a **client** that uses user-supplied compute (Modal, RunPod, Lightning, NVIDIA NIM) — host nothing yourself in v1.
8. **MCP security baseline**: do not be in the 36.7% SSRF-vulnerable / 43% command-injection / 41% no-auth buckets. Bake in input validation, rate limiting, and OAuth 2.1 from v0.1.

### 8. Concrete additive features
1. **`assess_confidence(structure)` resource**: returns a structured `ConfidenceReport` with per-residue pLDDT bucket; PAE matrix entropy + interface-PAE summary; IDR cross-check against DisProt 2026 + MobiDB; explicit hallucination flag where AF3 predicts order but DisProt says disorder; cross-model agreement score when ≥2 backends run. **This is the killer feature — build it first.**
2. **AFDB ↔ Foldseek agent loop**: "find structurally similar proteins to UniProt X with confident interfaces" — no existing MCP does this.
3. **Reproducibility manifest**: every output ships with model version, weights hash, seed, MSA database snapshot date, runtime env, GPU type.
4. **Notebook export**: emit a Jupyter/Colab notebook that re-runs the analysis end-to-end — major citation multiplier.
5. **AlphaMissense + AlphaFill + TMalphaFold integration**: sequence → structure → cofactor transplantation → membrane orientation → missense pathogenicity in one prompt.
6. **Therapeutic-area vertical packs (v0.6+)**: Antibody Pack (Chai-1 + DockQ + SAbDab), Kinase Pack (Boltz-2 affinity + KLIFS), GPCR Pack (membrane orientation + ligand co-folding).
7. **Education mode**: `--explain` flag that emits plain-language reading of pLDDT/PAE/ipTM for a specific result with citations. Plays to LLM strengths and to wet-lab biologist needs.
8. **Benchmarking harness**: public eval suite (CASP16 monomers, Polaris-ASAP set, Chai-2 antibody-target set, DisProt hallucination set). Ship results as a leaderboard.
9. **Wet-lab linkage**: optional Benchling / LabArchives ELN integration so an "in silico → wet-lab" handoff is one prompt.
10. **Model-agnostic fold-drift tracker**: alert users when a backend updates and a stored prediction would now differ.

## Details — Prioritized Improvement Recommendations

### Priority 1 (do now, before writing more code)
- **Resequence v0.1 around confidence interpretation**, not AFDB CRUD. `assess_confidence` works on AFDB pre-computed structures — no GPU needed.
- **Switch the skeleton to FastMCP 3.x with `@mcp.tool(task=True)` + Streamable HTTP + Redis-backed Docket.**
- **License-routing policy at the agent layer**: if `commercial=true`, exclude AF3 weights and Chai-2; route to Boltz-2 / OpenFold3 / Protenix / Chai-1.
- **Drop Chai-2 and Tamarind/Neurosnap wrappers entirely.** ToS risk, zero upside.
- **Add DisProt 2026 and MobiDB as first-class resources.** This is the IDR/hallucination backbone.

### Priority 2 (before v1.0)
- **Ship the Foldseek MCP wrapper** — no one has, the Foldseek web API is documented, it's a week of work, and AFDB↔Foldseek is your killer demo.
- **Add OpenFold3 + AQAffinity** as default commercial-safe backends.
- **Ship `reproducibility_manifest.json`** with every prediction.
- **Implement cross-model disagreement detection** (move from v0.4 to v0.3) as the second moat.

### Priority 3 (launch and post-launch)
- **Coordinate a 48-hour launch** (Tue–Thu, mid-morning PT): bioRxiv preprint + HN + r/bioinformatics + r/mcp + X/Twitter + LinkedIn + Smithery/Glama/PulseMCP/mcp.so listings via mcp-submit. Stockpile 100–200 baseline stars from your network first.
- **Submit to JOSS** once test coverage ≥90% and CI is green.
- **Publish a benchmark dataset on Zenodo with a DOI** (DisProt hallucination set, AF3-vs-Boltz2 disagreement maps) — a citable artifact independent of the software.
- **Apply for Schmidt Sciences VISS** first (lower adoption bar). Apply for **CZI EOSS Cycle 7** only once you have ≥1,000 stars and ≥3 citing groups. Apply for **Anthropic AI for Science** compute credits early.

### Priority 4 (sustainability)
- **Do not host inference.** Ship a client + reference Modal / RunPod / NVIDIA-NIM deploy manifests. Target GPU bill: ~$0.
- **Write a "if I disappear" governance doc** in the README. Solo OSS projects that survive have explicit handoff plans.
- **Partner with an academic group** for wet-lab validation and the first 3 citations.

## Recommendations
1. **Reframe the product**: "FoldCopilot is a confidence-aware interpretation copilot for protein structure predictions" — *not* "a wrapper for six folders." Wrappers are commodity; interpretation is the moat.
2. **Resequence**: confidence interpretation in v0.1, Foldseek in v0.2, Boltz-2 (only) in v0.3, ensemble-disagreement in v0.4, JOSS submission in v0.5.
3. **Cut** Chai-2, Tamarind/Neurosnap wrappers, hosted inference, AF3-as-default. **Add** OpenFold3, AQAffinity, DisProt 2026, MobiDB, Foldseek MCP, reproducibility manifests.
4. **Distribution**: JOSS + bioRxiv + Zenodo dataset > Show HN. Coordinated 48-hour launch across all MCP directories.
5. **Benchmarks that change recommendations**: if a polished competitor reaches 500+ stars in the next 60 days, pivot from generalist to a therapeutic vertical (antibodies or kinases) where domain depth matters more. If Anthropic ships a first-party AlphaFold Skill, double down on the DisProt/MobiDB/AlphaMissense integration layer — that is where you outflank them.

## Caveats
- The "94.7% success rate" of BioinfoMCP refers to classical NGS pipelines, not structure prediction; do not over-extrapolate.
- Chai-2's antibody hit rates (16–50% depending on metric) and Boltz-2's 0.66 FEP+ correlation are developer-reported; real-world benchmarks (deepmirror.ai on DHX9 and cGAS, Bugrova et al. on SMILES vs CCD sensitivity) show mixed results — frame all reported metrics as upper-bound estimates in user-facing copy.
- The HN-launch arXiv paper (2511.04453) is a preprint based on 138 AI/LLM repos; bioinformatics-specific HN dynamics may differ.
- Anthropic Fellows / CZI EOSS / Schmidt Sciences program details are accurate as of May 2026 but cycle dates change — verify before applying.
- ProteinMCP's "4 stars" snapshot is May 2026 and may change. The competitive landscape will look materially different in 6 months if a frontier lab ships a first-party fold copilot.
- The 22% AF3-IDR-hallucination rate is from a single arXiv preprint (2510.15939) on 72 DisProt proteins; the EBI's own "What AlphaFold 3 struggles with" page confirms the qualitative phenomenon but the quantitative figure may not generalize beyond the curated set.
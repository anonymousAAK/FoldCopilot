# AI Usage Disclosure

This document satisfies the JOSS mandatory AI usage disclosure requirement.

## AI Tools Used
- **Claude (Anthropic)**: Used for code generation, test writing, documentation drafting, and research synthesis during development.

## Scope of AI Assistance
- **Code generation**: Initial implementation of MCP tool wiring, client modules, test scaffolding, and utility functions.
- **Documentation**: README sections, docstrings, plan.md research synthesis, and this disclosure.
- **Research**: Web searches for competitive landscape analysis, model benchmarks, and protocol updates.

## Human Validation
- All AI-generated code was reviewed, tested, and modified by the maintainer.
- Architectural decisions (backend selection, license routing, confidence interpretation logic) were made by the maintainer.
- Scientific accuracy of confidence interpretation, IDR flagging, and hallucination detection was validated against published literature.
- All 205+ tests pass and were verified to test meaningful behavior.

## What Was NOT AI-Generated
- Project concept and framing (confidence-aware interpretation copilot)
- Backend selection criteria and license routing policy
- Scientific domain knowledge (pLDDT interpretation, PAE semantics, IDR biology)
- Benchmark dataset curation decisions

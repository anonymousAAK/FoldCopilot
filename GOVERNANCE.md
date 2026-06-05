# FoldCopilot Governance

## Project Status
FoldCopilot is maintained as an open-source project under the MIT License.

## If the Primary Maintainer Becomes Unavailable

This section ensures project continuity if the primary maintainer cannot continue.

### Succession Plan
1. **Repository access**: The repository is public on GitHub. Any fork can continue development independently.
2. **No vendor lock-in**: All backends are open-source (MIT/Apache-2.0). No proprietary dependencies.
3. **No hosted infrastructure**: FoldCopilot is a client — it hosts nothing. No servers to maintain, no GPU bills to pay.
4. **Self-contained**: All code, tests, and documentation are in the repository. No external secrets required for development.
5. **Test suite**: 205+ tests ensure correctness. CI runs on every push.

### Handoff Checklist
If transferring maintainership:
- [ ] Transfer GitHub repository ownership
- [ ] Update CITATION.cff with new maintainer
- [ ] Update PyPI package ownership
- [ ] Update Smithery/Glama/registry listings
- [ ] Notify citing groups (if any)

### Community Continuation
- Fork freely under MIT License
- Published benchmark datasets (Zenodo DOI) persist independently
- JOSS paper (if published) provides permanent citation record
- No API keys or secrets are needed to develop or test

## Decision Making
- Single maintainer model for now
- Major architectural decisions documented in `plan.md`
- Community input welcomed via GitHub Issues and Discussions

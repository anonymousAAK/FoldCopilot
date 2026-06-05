# Contributing to FoldCopilot

Thank you for your interest in contributing to FoldCopilot — a confidence-aware interpretation copilot for protein structure predictions.

## Getting Started

1. Fork and clone the repository
2. Install in development mode:
   ```bash
   pip install -e ".[dev]"
   ```
3. Run tests:
   ```bash
   pytest tests/ -q
   ```

## Development Setup

- Python 3.11+
- Dependencies: `pip install -e ".[dev,tasks,observability]"`
- Linter: `ruff check src/ tests/`
- Formatter: `ruff format src/ tests/`

## How to Contribute

### Bug Reports
Open an issue with:
- Steps to reproduce
- Expected vs actual behavior
- Python version and OS
- Backend being used (if prediction-related)

### Feature Requests
Open an issue describing:
- The problem your feature solves
- Proposed implementation (if any)
- Which plan.md section it relates to (if applicable)

### Pull Requests
1. Create a feature branch from `main`
2. Write tests for new functionality
3. Ensure all tests pass: `pytest tests/ -q`
4. Run linter: `ruff check src/ tests/`
5. Keep PRs focused — one feature or fix per PR
6. Update README.md if adding new tools

### Adding a New Prediction Backend
1. Create a client in `src/foldcopilot/clients/` following the `boltz2_client.py` pattern
2. Add the backend enum to `src/foldcopilot/models/prediction.py`
3. Add license routing in `BACKEND_LICENSES`
4. Wire into `src/foldcopilot/tools/predict.py`
5. Register the tool in `src/foldcopilot/server.py`
6. Add tests in `tests/`
7. Update `README.md` backend table

### Adding a New MCP Tool
1. Implement logic in the appropriate `src/foldcopilot/tools/` module
2. Register with `@mcp.tool()` in `server.py`
3. Write a clear multi-sentence docstring (Smithery scores tool descriptions)
4. Add tests
5. Update tool count in `.well-known/mcp.json`

## Code Style
- Follow existing patterns in the codebase
- Use type annotations
- Validate inputs at system boundaries (user input, external APIs)
- Use `httpx.AsyncClient` for HTTP calls with explicit timeouts
- Return structured dicts from tools, not raw strings

## Testing
- All tests in `tests/` directory
- Use `unittest.mock` for external API calls
- Mock at the HTTP level, not at the function level
- Test both success and error paths

## Code of Conduct
Be respectful, constructive, and inclusive. We follow the [Contributor Covenant](https://www.contributor-covenant.org/version/2/1/code_of_conduct/).

## License
By contributing, you agree that your contributions will be licensed under the MIT License.

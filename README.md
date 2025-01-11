# SR Assistant Prototype

[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

Systematic Review Assistant Prototype.

tl;dr:

```sh
❯ make
bootstrap                      How to install the `uv` project managent tool
python/insll                   Install project local Python managed by uv
install                        Install runtime dependencies
dev-install                    Install all dependencies for development
format                         Format and fix code with ruff
lint                           Lint code with ruff
typecheck                      Run MyPy
test/unit                      Run unit tests only
test/integration               Run integration tests only
test/all                       Run all tests (unit and integration)
test/cov                       Run unit tests with coverage
test/cov-all                   Run all tests with coverage
pre-commit                     Run pre-commit on all files
clean                          Clean Python build artifacts
dev-clean                      Clean everything including development artifacts
supabase/install               Install Supabase CLI
help                           Show this help
```

## Getting Started

1. [Install uv](https://docs.astral.sh/uv/getting-started/installation/),
   [uv](https://github.com/astral-sh/uv) is a modern Python package manager written in Rust. It has a
   `pip` compatible interface under `uv pip` command and resolves dependencies up to 10x faster than
   [Poetry](https://python-poetry.org/). `uv` is to Python package management what
   [ruff](https://docs.astral.sh/ruff/) is to formatting and linting. It can also manage project
   Python versions.

   ```sh
   ❯ make bootstrap
   Pick your poison:

       brew install uv
       pipx install uv
       pip install uv
       curl -LsSf https://astral.sh/uv/install.sh | sh
       cargo install --git https://github.com/astral-sh/uv uv
       winget install --id=astral-sh.uv  -e
       scoop install main/uv
       powershell -ExecutionPolicy ByPass -c 'irm https://astral.sh/uv/install.ps1 | iex'

   Or go container:  docker run ghcr.io/astral-sh/uv ...
   ```

2. Install for development: `make dev-install` (or `uv python install && uv sync`)

3. Run the app stages:

   ```sh
   uv run streamlit run src/step1/app_step1.py
   uv run streamlit run src/step2/app_step2.py
   ```

   This will change soon when we move to a [multi-page app](https://docs.streamlit.io/develop/concepts/multipage-apps/page-and-navigation) with one entrypoint (`make run`) for the MVP prototype.

## Testing

The project uses pytest for testing, with tests separated into two categories:

- **Unit Tests** (`tests/unit/`): Fast tests that don't require external services

  - Run with: `make test-unit`
  - Run with coverage: `make test-cov`
  - These run on every PR

- **Integration Tests** (`tests/integration/`): Slower tests that may require external services
  - Run with: `make test-integration`
  - These run on PRs but failures don't block merges

Additional commands:

- `make test-all` - Run all tests without coverage
- `make test-cov-all` - Run all tests with coverage

## Repo Overview and Conventions

- [uv](https://docs.astral.sh/uv/) is used for Python project management.
- [CICD](.github/) is with GitHub Actions, the usual stuff, no packaging or releases yet. Checks are allowed to fail.
- [ruff](https://docs.astral.sh/ruff/) is used for linting and formatting.
- [pre-commit](https://pre-commit.com)] is used for linting and formatting. Tests run on push.
- [mypy](https://mypy.readthedocs.io/en/stable/index.html) is used for type checking.
- [pytest](https://docs.pytest.org/en/latest/) is used for unit tests with coverage.

**TODO**:

- Local dev (docker compose)
- SupaBase and Qdrant free-tier cloud services.
- Install [Renovate GitHub App](https://github.com/apps/renovate) for automatic dependency updates. Like Dependabot, but better. Config in `.github/renovate.json`.

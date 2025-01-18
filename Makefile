.PHONY: help bootstrap python python.list install install.prod format lint ruff.fix typecheck test.unit test.integration test.all clean clean.lean security supabase.cli supabase.dbdev submodules

.DEFAULT_GOAL := help

bootstrap:  ## How to install the `uv` project managent tool
	@echo "Pick your poison:"
	@echo
	@echo "    brew install uv"
	@echo "    pipx install uv"
	@echo "    pip install uv"
	@echo "    curl -LsSf https://astral.sh/uv/install.sh | sh"
	@echo "    cargo install --git https://github.com/astral-sh/uv uv"
	@echo "    winget install --id=astral-sh.uv  -e"
	@echo "    scoop install main/uv"
	@echo "    powershell -ExecutionPolicy ByPass -c 'irm https://astral.sh/uv/install.ps1 | iex'"
	@echo
	@echo "Or go container:  docker run ghcr.io/astral-sh/uv ..."

python.insall: ## Install project local Python managed by uv
	uv python install

python.list: ## Show Python versions available with uv
	uv python list

install.prod: python  ## Install runtime dependencies
	uv sync --no-dev

deno.cli: ## Install Deno CLI, required by the VSCode extension used to develop edge function
	# Goes in $DXG_DATA_HOME/deno/bin, add to path, let the installer do it, or
	# symlink to $HOME/.local/bin or whatever it's on macoS.
	curl -fsSL https://deno.land/install.sh | sh
	deno --version
	@echo "Install the deno vscode extension (canary), use the egde-functions workspace"

install: python ## Install all dependencies for development
	uv sync
	uv run pre-commit install
	uv run pre-commit install --hook-type pre-push
	@echo "If needed, also run 'make deno.cli' if developing edge functions"

format:  ## Format and fix code with ruff
	uv run ruff format
	uv run ruff check --fix

ruff.fix:
	uv run check --fix

lint:  ## Lint code with ruff
	uv run ruff check

typecheck:  ## Run MyPy
	uv run mypy

security:
	uv sync --group security
	uv run pip-audit
	uv run ruff --select="S"

test.unit: ## Run unit tests with coverage
	uv run pytest tests/unit --cov=src --cov-report=term --cov-report=html
	@echo "Coverage report available at htmlcov/index.html"

test.integration: ## Run integration tests only
	uv run pytest tests/integration -m integration

test.all: ## Run all tests with coverage
	uv run pytest tests --cov=src --cov-report=term --cov-report=html
	@echo "Coverage report available at htmlcov/index.html"

pre-commit:  ## Run pre-commit on all files
	uv run pre-commit run --all-files

clean: clean.dev  ## Clean Python build artifacts
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
	find . -type f -name "*.pyd" -delete
	find . -type f -name ".coverage" -delete
	find . -type d -name "*.egg" -exec rm -rf {} +
	find . -type d -name ".eggs" -exec rm -rf {} +
	find . -type d -name ".pytest_cache" -exec rm -rf {} +
	find . -type d -name ".ruff_cache" -exec rm -rf {} +
	find . -type d -name ".mypy_cache" -exec rm -rf {} +


clean.dev: clean  ## Clean everything including development artifacts
	rm -rf .venv/
	rm -rf htmlcov/
	rm -rf .coverage*
	rm -rf .pytest_cache/
	rm -rf .mypy_cache/
	rm -rf .ruff_cache/
	rm -rf .pre-commit-cache/

supabase.cli:  ## Install Supabase CLI
	brew install supabase/tap/supabase

supabase.dbdev: ## Install dbdev DB package manager
	brew install supabase/tap/dbdev

submodules:  ## Sync git submodules (database-build, etc.)
	git submodule update --init --recursive
	@echo ">> tools/database-build synced: in-browser WASM AI assisted Postgres IDE"

help:  ## Show this help
	@grep -E '^[a-z/.:A-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

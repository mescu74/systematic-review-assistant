.PHONY: help bootstrap python show-python install dev-install format lint type-check test-unit test-integration test-all test-cov test-cov-all clean dev-clean

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

python/insll: ## Install project local Python managed by uv
	uv python install

python/liw5: ## Show Python versions available with uv
	uv python list

install: python  ## Install runtime dependencies
	uv sync --no-group dev

dev-install: python  ## Install all dependencies for development
	uv sync
	pre-commit install
	pre-commit install --hook-type pre-push

format:  ## Format and fix code with ruff
	uv run ruff format
	uv run ruff check --fix

lint:  ## Lint code with ruff
	uv run ruff check

typecheck:  ## Run MyPy
	uv run mypy

test/unit: ## Run unit tests only
	uv run pytest tests/unit

test/integration: ## Run integration tests only
	uv run pytest tests/integration -m integration

test/all: ## Run all tests (unit and integration)
	uv run pytest tests

test/cov: ## Run unit tests with coverage
	uv run pytest tests/unit --cov=src --cov-report=term --cov-report=html
	@echo "Coverage report available at htmlcov/index.html"

test/cov-all: ## Run all tests with coverage
	uv run pytest tests --cov=src --cov-report=term --cov-report=html
	@echo "Coverage report available at htmlcov/index.html"

pre-commit:  ## Run pre-commit on all files
	pre-commit run --all-files

clean:  ## Clean Python build artifacts
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


dev-clean: clean  ## Clean everything including development artifacts
	rm -rf .venv/
	rm -rf htmlcov/
	rm -rf .coverage*
	rm -rf .pytest_cache/
	rm -rf .mypy_cache/
	rm -rf .ruff_cache/
	rm -rf .pre-commit-cache/

supabase/install:  ## Install Supabase CLI
	brew install supabase/tap/supabase

help:  ## Show this help
	@grep -E '^[a-z/:A-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

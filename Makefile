# Neurofuzzy task runner. Run `make help` for the list of targets.
.DEFAULT_GOAL := help
PYTHON ?= python
RESULTS ?= benchmarks/results
FIGURES ?= benchmarks/figures
SEEDS ?= 10

.PHONY: help install install-dev test lint format typecheck cov benchmark figures report repro docker clean

help: ## Show this help.
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-14s\033[0m %s\n", $$1, $$2}'

install: ## Install the package (runtime only).
	$(PYTHON) -m pip install -e .

install-dev: ## Install with dev + viz extras.
	$(PYTHON) -m pip install -e ".[dev]"

test: ## Run the test suite.
	$(PYTHON) -m pytest

cov: ## Run tests with coverage.
	$(PYTHON) -m pytest --cov=neurofuzzy --cov-report=term-missing

lint: ## Lint with ruff.
	ruff check neurofuzzy tests experiments

format: ## Auto-format with ruff.
	ruff format neurofuzzy tests experiments

typecheck: ## Static type check with mypy.
	mypy

benchmark: ## Run the full benchmark (neurofuzzy + baselines) and aggregate.
	$(PYTHON) -m neurofuzzy.benchmark --all --seeds $(SEEDS) --out $(RESULTS)

figures: ## Generate all benchmark figures.
	$(PYTHON) -m neurofuzzy.visualize --results $(RESULTS) --out $(FIGURES)

report: benchmark figures ## Reproduce results + figures end to end.

repro: install-dev report ## One-command reproduction from a clean checkout.

clean: ## Remove caches and build artifacts.
	rm -rf build dist *.egg-info .pytest_cache .ruff_cache .mypy_cache
	find . -type d -name __pycache__ -prune -exec rm -rf {} +

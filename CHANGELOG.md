# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this project adheres
to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2026-06-29

### Added
- **Unified benchmark harness** (`neurofuzzy.benchmark`): resumable multi-seed
  runner with per-seed metric logging, bootstrap 95% confidence intervals,
  paired Wilcoxon signed-rank significance tests, and auto-generated Markdown
  and LaTeX reports plus a reproducibility manifest (package versions, git
  commit, configuration).
- **Reference baselines** (`neurofuzzy.baselines`): mean-of-features,
  best-single-feature, and linear-regression baselines with a common API, so the
  repository is usable as a standalone benchmark.
- **Publication-quality figures** (`neurofuzzy.visualize`): metric comparison
  bars with CIs, train/test generalization-gap chart, prediction scatter plots,
  and learned-fuzzy-set visualizations.
- Console entry points `neurofuzzy-benchmark` and `neurofuzzy-figures`.
- Tooling: `Makefile`, `Dockerfile`, pinned `requirements-lock.txt`,
  `.pre-commit-config.yaml`, and ruff/mypy/coverage configuration.
- Project metadata and community health files: `CITATION.cff`, `SECURITY.md`,
  `CODE_OF_CONDUCT.md`, `ROADMAP.md`, issue/PR templates.
- Committed, regenerable benchmark artifacts under `benchmarks/`.

### Changed
- **Vectorized `FuzzyController`** (~50x faster batch prediction) while keeping
  results numerically identical to the previous per-sample implementation
  (verified bit-for-bit over hundreds of random controllers). Output fuzzy sets
  are now precomputed once at construction time.
- Replaced the placeholder README results table with **real, reproducible**
  numbers (mean with 95% CIs over 10 seeds) and an honest discussion of the
  observed train/test generalization gap.
- Trimmed unused dependencies (`scikit-learn`, `tqdm`); made `matplotlib` an
  optional `viz` extra.
- Fixed `environment.yml` to install the local package instead of a non-existent
  PyPI distribution.
- `NeurofuzzyModel.fit` now uses the standard `logging` module instead of `print`.

### Notes
- All reported metrics are produced by the benchmark harness and depend on the
  documented optimizer configuration and the small benchmark sample sizes.
  No results are hard-coded or fabricated.

## [1.0.0] - 2023
- Initial Python port of the legacy Java prototype accompanying the
  *Data & Knowledge Engineering* (2023) article.

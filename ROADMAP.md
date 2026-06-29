# Roadmap

This roadmap captures planned directions. Contributions toward any item are
welcome — see [CONTRIBUTING.md](CONTRIBUTING.md) and the "good first issue" label.

## Near term
- Additional word/sentence-similarity benchmark datasets with automatic
  downloaders (e.g., RG-65, WordSim-353, SimLex-999) and license metadata.
- More baselines (ridge/lasso, gradient boosting, simple MLP) under the same
  evaluation protocol.
- Cross-validation and nested-CV protocols to complement the repeated holdout.

## Medium term
- Pluggable defuzzification and membership-function families with a registry.
- Regularization / model-selection to reduce the observed train/test gap on
  small datasets.
- A hosted documentation site (MkDocs Material) on GitHub Pages.

## Longer term
- A public leaderboard generated from the benchmark artifacts.
- Optional integration points for modern neural sentence encoders as feature
  providers.

Have an idea? Open a feature request using the issue template.

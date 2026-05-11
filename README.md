# 🧠 Neurofuzzy Semantic Similarity Measurement

> Martinez-Gil, J., Mokadem, R., Küng, J., & Hameurlain, A. (2023). **Neurofuzzy semantic similarity measurement.** *Data & Knowledge Engineering*, 145, 102155. https://doi.org/10.1016/j.datak.2023.102155

[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![CI](https://github.com/jorge-martinez-gil/dke2023/actions/workflows/ci.yml/badge.svg)](https://github.com/jorge-martinez-gil/dke2023/actions/workflows/ci.yml)
[![Paper](https://img.shields.io/badge/Paper-DKE%202023-green.svg)](https://doi.org/10.1016/j.datak.2023.102155)
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.datak.2023.102155-blue)](https://doi.org/10.1016/j.datak.2023.102155)

## Abstract

This repository provides a fully reproducible Python implementation of neurofuzzy semantic similarity measurement. The approach combines neural similarity features with a Mamdani fuzzy inference controller optimized by Differential Evolution. The resulting system maps heterogeneous neural signals to robust semantic similarity estimates in the range `[0, 1]`.

## Key Contributions

- Clean Python package (`neurofuzzy`) replacing the legacy Java prototype
- Pure NumPy fuzzy controller with explicit parameterization and centroid defuzzification
- Differential Evolution optimization using SciPy for reproducible parameter learning
- Reproducible experiment CLI for MC and Geresid datasets with multi-run statistics
- Unit tests and CI-ready workflow for reliability

## Architecture

```text
Word Pairs → [Neural Encoder] → Feature Vectors → [Fuzzy Controller] → Similarity Score
                                      ↑
                           [Differential Evolution Optimizer]
```

## Quick Start

```bash
pip install -e .
```

```python
from neurofuzzy.data_loader import load_dataset
from neurofuzzy.model import NeurofuzzyModel

X, y = load_dataset("datasets/mc.txt")
model = NeurofuzzyModel(maxiter=50, popsize=10, verbose=False).fit(X, y)
pred = model.predict(X[:5])
```

```bash
python experiments/run_experiment.py --dataset mc --runs 10 --output results/
```

## Datasets

| Dataset | Samples | Features used | Description |
|---|---:|---:|---|
| MC (`datasets/mc.txt`) | 30 | 4 | Miller & Charles-style benchmark pairs with neural features |
| Geresid (`datasets/geresid.txt`) | 50 | 4 | Geresid benchmark pairs with neural features |

## Reproducibility

```bash
conda env create -f environment.yml
conda activate neurofuzzy
python experiments/run_experiment.py --dataset mc --runs 30
python experiments/run_experiment.py --dataset geresid --runs 30
```

## Project Structure

```text
dke2023/
├── datasets/
│   ├── mc.txt
│   └── geresid.txt
├── neurofuzzy/
│   ├── __init__.py
│   ├── data_loader.py
│   ├── fuzzy_controller.py
│   ├── metrics.py
│   ├── model.py
│   └── optimizer.py
├── experiments/
│   └── run_experiment.py
├── tests/
│   ├── test_data_loader.py
│   ├── test_metrics.py
│   └── test_model.py
├── src/dke2023/
│   └── dke2023.java
├── .github/workflows/ci.yml
├── CONTRIBUTING.md
├── environment.yml
├── pyproject.toml
└── requirements.txt
```

## Results

| Dataset | Pearson r | Spearman ρ |
|---|---:|---:|
| MC | ~0.89 | ~0.87 |
| Geresid | ~0.78 | ~0.76 |

> Values are representative placeholders aligned with the publication context; run the experiment script to generate current reproducible values in your environment.

## Citation

```bibtex
@article{martinez2023neurofuzzy,
  author       = {Jorge Martinez-Gil and Riad Mokadem and Josef K{"u}ng and Abdelkader Hameurlain},
  title        = {Neurofuzzy semantic similarity measurement},
  journal      = {Data \& Knowledge Engineering},
  volume       = {145},
  pages        = {102155},
  year         = {2023},
  doi          = {10.1016/j.datak.2023.102155},
  url          = {https://doi.org/10.1016/j.datak.2023.102155}
}
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

## License

MIT License. See [LICENSE](LICENSE).

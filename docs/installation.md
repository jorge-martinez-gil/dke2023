# Installation

## With pip (editable)

```bash
git clone https://github.com/jorge-martinez-gil/dke2023.git
cd dke2023
pip install -e .            # runtime
pip install -e ".[dev]"     # + tests, linting, type-checking, figures
```

## With conda

```bash
conda env create -f environment.yml
conda activate neurofuzzy
```

## Exact pins (for reproduction)

```bash
pip install -r requirements-lock.txt
pip install -e .
```

## Docker

```bash
docker build -t neurofuzzy .
docker run --rm -v "$PWD/benchmarks:/app/benchmarks" neurofuzzy
```

## Requirements

Python 3.9+ and NumPy, SciPy, pandas. Figures additionally require matplotlib
(`pip install -e ".[viz]"`).

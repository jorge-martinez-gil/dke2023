# Quickstart

## Fit and predict

```python
from neurofuzzy import NeurofuzzyModel
from neurofuzzy.data_loader import load_dataset

X, y = load_dataset("datasets/mc.txt")
model = NeurofuzzyModel(maxiter=40, popsize=8, random_state=0, verbose=False).fit(X, y)
print(model.predict(X[:5]))
print(model.evaluate(X, y))
```

## Run the benchmark

```bash
python -m neurofuzzy.benchmark --all --seeds 10 --out benchmarks/results
python -m neurofuzzy.visualize --results benchmarks/results --out benchmarks/figures
```

The benchmark is **resumable** and respects `--time-budget <seconds>`, so you can
complete a long run in several short sessions. See the
[benchmark protocol](benchmark.md) for details.

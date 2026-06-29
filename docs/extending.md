# Extending

## Add a dataset

1. Place a text file in `datasets/` using the [dataset format](datasets.md).
2. Register it in `neurofuzzy/benchmark.py`:

   ```python
   DATASET_FILES = {"mc": "mc.txt", "geresid": "geresid.txt", "mine": "mine.txt"}
   ```
3. Run `python -m neurofuzzy.benchmark --all --datasets mine`.

See [`examples/custom_dataset.py`](https://github.com/jorge-martinez-gil/dke2023/blob/main/examples/custom_dataset.py).

## Add a method / baseline

Implement the minimal API and register it:

```python
import numpy as np
from neurofuzzy.baselines import BaseSimilarityBaseline, BASELINES

class MyBaseline(BaseSimilarityBaseline):
    name = "my-baseline"
    def fit(self, X, y):
        # learn parameters from the training split
        return self
    def predict(self, X):
        return np.clip(..., 0.0, 1.0)

BASELINES[MyBaseline.name] = MyBaseline
```

It is then automatically included in the benchmark, confidence intervals, and
significance tests. See [`examples/add_baseline.py`](https://github.com/jorge-martinez-gil/dke2023/blob/main/examples/add_baseline.py).

## Add a metric

Add a function to `neurofuzzy/metrics.py` and include it in `evaluate_all`. Add a
unit test asserting its value on a known input.

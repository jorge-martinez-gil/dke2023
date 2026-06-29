# Datasets

## Format

Each dataset is a plain-text file with one word/text pair per line:

```text
ground_truth, feature_1, feature_2, feature_3, feature_4[, ...]
```

- The first column is the human-annotated similarity (the target).
- The loader uses the next four columns as input features (`n_features=4`).
- Extra columns are ignored.

## Bundled datasets

| Dataset | File | Pairs | Notes |
|---|---|---:|---|
| MC | `datasets/mc.txt` | 30 | Miller–Charles-style pairs; rows carry 7 features, the first 4 are used. |
| Geresid | `datasets/geresid.txt` | 50 | Benchmark pairs with 4 features. |

## Loading

```python
from neurofuzzy.data_loader import load_dataset
X, y = load_dataset("datasets/mc.txt", n_features=4)
```

The loader validates row width and raises a clear `ValueError` on malformed rows.

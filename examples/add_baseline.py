"""Implement a custom baseline and benchmark it against the built-ins.

Any object with ``fit(X, y)`` and ``predict(X)`` works. Registering it in
``neurofuzzy.baselines.BASELINES`` makes it available to the benchmark CLI.

Run::

    python examples/add_baseline.py
"""

from __future__ import annotations

import numpy as np

from neurofuzzy.baselines import BaseSimilarityBaseline
from neurofuzzy.data_loader import load_dataset, train_test_split
from neurofuzzy.metrics import evaluate_all


class MedianFeatureBaseline(BaseSimilarityBaseline):
    """Predict the per-pair median of the input features (robust to outliers)."""

    name = "median-of-features"

    def predict(self, X: np.ndarray) -> np.ndarray:
        return np.clip(np.median(np.asarray(X, dtype=float), axis=1), 0.0, 1.0)


def main() -> None:
    X, y = load_dataset("datasets/mc.txt")
    X_tr, X_te, y_tr, y_te = train_test_split(X, y, train_ratio=0.6, random_state=0)

    model = MedianFeatureBaseline().fit(X_tr, y_tr)
    metrics = evaluate_all(y_te, model.predict(X_te))
    print(f"{model.name} test metrics:", {k: round(v, 4) for k, v in metrics.items()})

    print(
        "\nTo include it in the full benchmark, add this to "
        "neurofuzzy/baselines.py:\n"
        "    BASELINES[MedianFeatureBaseline.name] = MedianFeatureBaseline\n"
        "then run:  python -m neurofuzzy.benchmark --all"
    )


if __name__ == "__main__":
    main()

"""Evaluate the model on your own dataset file.

Dataset format: one pair per line, comma-separated, as

    ground_truth, feature_1, feature_2, feature_3, feature_4

Run::

    python examples/custom_dataset.py path/to/your_dataset.txt
"""

from __future__ import annotations

import sys

from neurofuzzy import NeurofuzzyModel
from neurofuzzy.data_loader import load_dataset, train_test_split


def main(path: str) -> None:
    X, y = load_dataset(path, n_features=4)
    X_tr, X_te, y_tr, y_te = train_test_split(X, y, train_ratio=0.6, random_state=0)
    model = NeurofuzzyModel(
        maxiter=40, popsize=8, train_ratio=1.0, random_state=0, verbose=False
    ).fit(X_tr, y_tr)
    print(f"Loaded {X.shape[0]} pairs from {path}")
    print("Test metrics:", {k: round(v, 4) for k, v in model.evaluate(X_te, y_te).items()})


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(1)
    main(sys.argv[1])

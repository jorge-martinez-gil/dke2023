"""Minimal end-to-end example: load data, fit, predict, evaluate.

Run from the repository root::

    python examples/minimal_example.py
"""

from __future__ import annotations

from neurofuzzy import NeurofuzzyModel
from neurofuzzy.data_loader import load_dataset, train_test_split


def main() -> None:
    X, y = load_dataset("datasets/mc.txt")
    X_train, X_test, y_train, y_test = train_test_split(X, y, train_ratio=0.6, random_state=0)

    # Small, fast configuration so the example finishes quickly.
    model = NeurofuzzyModel(maxiter=40, popsize=8, train_ratio=1.0, random_state=0, verbose=False)
    model.fit(X_train, y_train)

    print("Predictions (first 5 test pairs):", model.predict(X_test[:5]).round(3))
    metrics = model.evaluate(X_test, y_test)
    for name, value in metrics.items():
        print(f"  {name:24s} {value:.4f}")


if __name__ == "__main__":
    main()

"""Reference baseline methods for semantic-similarity benchmarking.

These deliberately simple, dependency-light baselines turn the repository into a
usable *benchmark*: any new method can be compared against the same train/test
protocol and metrics. None of the numbers here are hard-coded -- every score is
computed from the data at run time.

All baselines follow a minimal scikit-learn-style API: ``fit(X, y)`` returns
``self`` and ``predict(X)`` returns an array of similarity scores in ``[0, 1]``.
"""

from __future__ import annotations

from typing import Callable

import numpy as np

from .metrics import pearson


class BaseSimilarityBaseline:
    """Common interface for similarity baselines."""

    #: Human-readable, citable name used in reports and figures.
    name: str = "base"

    def fit(self, X: np.ndarray, y: np.ndarray) -> BaseSimilarityBaseline:
        """Fit the baseline (no-op for parameter-free baselines)."""
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:  # pragma: no cover - abstract
        raise NotImplementedError


class MeanFeatureBaseline(BaseSimilarityBaseline):
    """Predict the arithmetic mean of the input similarity features.

    A parameter-free ensemble: it assumes every neural feature is an equally
    valid similarity estimate and averages them.
    """

    name = "mean-of-features"

    def predict(self, X: np.ndarray) -> np.ndarray:
        arr = np.asarray(X, dtype=float)
        return np.clip(arr.mean(axis=1), 0.0, 1.0)


class BestSingleFeatureBaseline(BaseSimilarityBaseline):
    """Select, on the training split, the single most correlated feature.

    The chosen feature index is then used unchanged at prediction time. This is
    the classic "pick the best individual measure" baseline.
    """

    name = "best-single-feature"

    def __init__(self) -> None:
        self.feature_index_: int | None = None

    def fit(self, X: np.ndarray, y: np.ndarray) -> BestSingleFeatureBaseline:
        arr = np.asarray(X, dtype=float)
        target = np.asarray(y, dtype=float)
        correlations = [abs(pearson(target, arr[:, j])) for j in range(arr.shape[1])]
        self.feature_index_ = int(np.argmax(correlations))
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        if self.feature_index_ is None:
            raise RuntimeError("BestSingleFeatureBaseline must be fitted before predict().")
        arr = np.asarray(X, dtype=float)
        return np.clip(arr[:, self.feature_index_], 0.0, 1.0)


class LinearRegressionBaseline(BaseSimilarityBaseline):
    """Ordinary least-squares regression of the target on the features.

    Predictions are clipped to ``[0, 1]``. Implemented with ``numpy.linalg``
    only, so it adds no heavy dependencies.
    """

    name = "linear-regression"

    def __init__(self) -> None:
        self.coef_: np.ndarray | None = None

    def fit(self, X: np.ndarray, y: np.ndarray) -> LinearRegressionBaseline:
        arr = np.asarray(X, dtype=float)
        target = np.asarray(y, dtype=float)
        design = np.hstack([np.ones((arr.shape[0], 1)), arr])
        self.coef_, *_ = np.linalg.lstsq(design, target, rcond=None)
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        if self.coef_ is None:
            raise RuntimeError("LinearRegressionBaseline must be fitted before predict().")
        arr = np.asarray(X, dtype=float)
        design = np.hstack([np.ones((arr.shape[0], 1)), arr])
        return np.clip(design @ self.coef_, 0.0, 1.0)


#: Registry mapping baseline names to zero-argument factories.
BASELINES: dict[str, Callable[[], BaseSimilarityBaseline]] = {
    MeanFeatureBaseline.name: MeanFeatureBaseline,
    BestSingleFeatureBaseline.name: BestSingleFeatureBaseline,
    LinearRegressionBaseline.name: LinearRegressionBaseline,
}


def get_baseline(name: str) -> BaseSimilarityBaseline:
    """Instantiate a baseline by name.

    Parameters
    ----------
    name : str
        One of the keys in :data:`BASELINES`.

    Raises
    ------
    KeyError
        If ``name`` is not a registered baseline.
    """
    if name not in BASELINES:
        raise KeyError(f"Unknown baseline {name!r}. Available: {sorted(BASELINES)}")
    return BASELINES[name]()

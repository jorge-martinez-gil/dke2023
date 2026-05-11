"""Differential Evolution optimizer for the neurofuzzy controller."""

from __future__ import annotations

from typing import Callable, Optional

import numpy as np
from scipy.optimize import OptimizeResult, differential_evolution

from .fuzzy_controller import FuzzyController
from .metrics import pearson


def optimize(
    X_train: np.ndarray,
    y_train: np.ndarray,
    n_params: int = 40,
    maxiter: int = 200,
    popsize: int = 15,
    seed: int = 42,
    callback: Optional[Callable[..., bool]] = None,
) -> OptimizeResult:
    """Optimize fuzzy parameters with Differential Evolution.

    Parameters
    ----------
    X_train : numpy.ndarray
        Training feature matrix with shape ``(n_samples, 4)``.
    y_train : numpy.ndarray
        Training targets with shape ``(n_samples,)``.
    n_params : int, default=40
        Number of parameters to optimize.
    maxiter : int, default=200
        Maximum number of DE iterations.
    popsize : int, default=15
        Population size multiplier.
    seed : int, default=42
        Random seed.
    callback : callable or None, default=None
        Optional callback passed to SciPy DE.

    Returns
    -------
    scipy.optimize.OptimizeResult
        Optimization result with ``.x`` containing best parameters.
    """

    X_arr = np.asarray(X_train, dtype=float)
    y_arr = np.asarray(y_train, dtype=float)

    def objective(params: np.ndarray) -> float:
        controller = FuzzyController(params)
        preds = controller.predict(X_arr)
        return -pearson(y_arr, preds)

    bounds = [(0.0, 1.0)] * n_params
    return differential_evolution(
        objective,
        bounds=bounds,
        maxiter=maxiter,
        popsize=popsize,
        seed=seed,
        callback=callback,
        polish=True,
        updating='deferred',
    )

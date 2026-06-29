"""Evaluation metrics for semantic similarity predictions."""

from __future__ import annotations

import warnings

import numpy as np
from scipy.stats import pearsonr, spearmanr


def _to_arrays(y_true: np.ndarray, y_pred: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    a = np.asarray(y_true, dtype=float).ravel()
    b = np.asarray(y_pred, dtype=float).ravel()
    if a.shape != b.shape:
        raise ValueError(f"Shape mismatch: y_true {a.shape} vs y_pred {b.shape}")
    return a, b


def pearson(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Compute Pearson correlation coefficient.

    Parameters
    ----------
    y_true : numpy.ndarray
        Ground-truth values.
    y_pred : numpy.ndarray
        Predicted values.

    Returns
    -------
    float
        Pearson correlation in ``[-1, 1]``.
    """
    y_t, y_p = _to_arrays(y_true, y_pred)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        value = pearsonr(y_t, y_p).statistic
    return float(0.0 if np.isnan(value) else value)


def spearman(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Compute Spearman rank correlation coefficient.

    Parameters
    ----------
    y_true : numpy.ndarray
        Ground-truth values.
    y_pred : numpy.ndarray
        Predicted values.

    Returns
    -------
    float
        Spearman correlation in ``[-1, 1]``.
    """
    y_t, y_p = _to_arrays(y_true, y_pred)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        value = spearmanr(y_t, y_p).statistic
    return float(0.0 if np.isnan(value) else value)


def cosine_similarity(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Compute cosine similarity between two vectors.

    Parameters
    ----------
    y_true : numpy.ndarray
        Ground-truth values.
    y_pred : numpy.ndarray
        Predicted values.

    Returns
    -------
    float
        Cosine similarity in ``[-1, 1]``.
    """
    y_t, y_p = _to_arrays(y_true, y_pred)
    denom = float(np.linalg.norm(y_t) * np.linalg.norm(y_p))
    if denom <= 1e-12:
        return 0.0
    return float(np.dot(y_t, y_p) / denom)


def mean_absolute_error(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Compute mean absolute error.

    Parameters
    ----------
    y_true : numpy.ndarray
        Ground-truth values.
    y_pred : numpy.ndarray
        Predicted values.

    Returns
    -------
    float
        MAE value.
    """
    y_t, y_p = _to_arrays(y_true, y_pred)
    return float(np.mean(np.abs(y_t - y_p)))


def root_mean_squared_error(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Compute root mean squared error.

    Parameters
    ----------
    y_true : numpy.ndarray
        Ground-truth values.
    y_pred : numpy.ndarray
        Predicted values.

    Returns
    -------
    float
        RMSE value.
    """
    y_t, y_p = _to_arrays(y_true, y_pred)
    return float(np.sqrt(np.mean((y_t - y_p) ** 2)))


def evaluate_all(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    """Compute all evaluation metrics.

    Parameters
    ----------
    y_true : numpy.ndarray
        Ground-truth values.
    y_pred : numpy.ndarray
        Predicted values.

    Returns
    -------
    dict
        Metrics dictionary with Pearson, Spearman, cosine similarity,
        MAE, and RMSE.
    """
    return {
        "pearson": pearson(y_true, y_pred),
        "spearman": spearman(y_true, y_pred),
        "cosine_similarity": cosine_similarity(y_true, y_pred),
        "mean_absolute_error": mean_absolute_error(y_true, y_pred),
        "root_mean_squared_error": root_mean_squared_error(y_true, y_pred),
    }

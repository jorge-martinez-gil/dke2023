"""Data loading and splitting utilities for neurofuzzy experiments."""

from __future__ import annotations

from pathlib import Path

import numpy as np


def load_dataset(filepath: str | Path, n_features: int = 4) -> tuple[np.ndarray, np.ndarray]:
    """Load a similarity dataset from a text file.

    Parameters
    ----------
    filepath : str or pathlib.Path
        Path to a dataset where each non-empty line contains comma-separated
        numeric values: ``ground_truth, feature_1, ..., feature_n``.
    n_features : int, default=4
        Number of feature columns to load after the ground-truth column.

    Returns
    -------
    tuple of numpy.ndarray
        ``(X, y)`` where ``X`` has shape ``(n_samples, n_features)`` and
        ``y`` has shape ``(n_samples,)``.

    Raises
    ------
    ValueError
        If a row has fewer columns than required.
    """
    rows: list[list[float]] = []
    path = Path(filepath)
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            parts = [token.strip() for token in line.split(",") if token.strip()]
            if len(parts) < n_features + 1:
                raise ValueError(
                    f"Invalid row in {filepath!s}: expected at least "
                    f"{n_features + 1} columns, got {len(parts)}"
                )
            rows.append([float(value) for value in parts[: n_features + 1]])

    data = np.asarray(rows, dtype=float)
    if data.ndim != 2 or data.shape[1] != n_features + 1:
        raise ValueError(
            f"Invalid dataset shape for {filepath!s}: expected "
            f"(?, {n_features + 1}), got {data.shape}"
        )

    y = data[:, 0]
    X = data[:, 1 : n_features + 1]
    return X, y


def train_test_split(
    X: np.ndarray,
    y: np.ndarray,
    train_ratio: float = 0.6,
    random_state: int | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Shuffle and split features/targets into train and test subsets.

    Parameters
    ----------
    X : numpy.ndarray
        Feature matrix with shape ``(n_samples, n_features)``.
    y : numpy.ndarray
        Target vector with shape ``(n_samples,)``.
    train_ratio : float, default=0.6
        Fraction of samples assigned to the training split.
    random_state : int or None, default=None
        Seed for reproducible shuffling.

    Returns
    -------
    tuple of numpy.ndarray
        ``(X_train, X_test, y_train, y_test)``.

    Raises
    ------
    ValueError
        If shapes are inconsistent or ``train_ratio`` is invalid.
    """
    if X.shape[0] != y.shape[0]:
        raise ValueError("X and y must contain the same number of samples")
    if not 0.0 < train_ratio < 1.0:
        raise ValueError("train_ratio must be in the open interval (0, 1)")

    rng = np.random.default_rng(random_state)
    indices = np.arange(X.shape[0])
    rng.shuffle(indices)

    split_idx = int(X.shape[0] * train_ratio)
    train_idx = indices[:split_idx]
    test_idx = indices[split_idx:]

    return X[train_idx], X[test_idx], y[train_idx], y[test_idx]

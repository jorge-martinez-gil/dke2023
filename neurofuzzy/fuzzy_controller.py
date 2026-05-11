"""Parameterized Mamdani fuzzy inference controller in pure NumPy."""

from __future__ import annotations

import numpy as np


def trimf(x: np.ndarray | float, a: float, b: float, c: float) -> np.ndarray:
    """Compute triangular membership values.

    Parameters
    ----------
    x : numpy.ndarray or float
        Input values.
    a : float
        Left foot.
    b : float
        Peak.
    c : float
        Right foot.

    Returns
    -------
    numpy.ndarray
        Membership values in ``[0, 1]``.
    """
    x_arr = np.asarray(x, dtype=float)
    denom_left = b - a
    denom_right = c - b
    left = np.zeros_like(x_arr, dtype=float)
    right = np.zeros_like(x_arr, dtype=float)
    if abs(denom_left) > 1e-12:
        left = (x_arr - a) / denom_left
    if abs(denom_right) > 1e-12:
        right = (c - x_arr) / denom_right
    return np.clip(np.minimum(left, right), 0.0, 1.0)


def trapmf(x: np.ndarray | float, a: float, b: float, c: float, d: float) -> np.ndarray:
    """Compute trapezoidal membership values.

    Parameters
    ----------
    x : numpy.ndarray or float
        Input values.
    a : float
        Left foot.
    b : float
        Left shoulder.
    c : float
        Right shoulder.
    d : float
        Right foot.

    Returns
    -------
    numpy.ndarray
        Membership values in ``[0, 1]``.
    """
    x_arr = np.asarray(x, dtype=float)
    rise = np.ones_like(x_arr, dtype=float)
    fall = np.ones_like(x_arr, dtype=float)
    if abs(b - a) > 1e-12:
        rise = (x_arr - a) / (b - a)
    if abs(d - c) > 1e-12:
        fall = (d - x_arr) / (d - c)
    top = np.minimum(np.minimum(rise, 1.0), np.minimum(fall, 1.0))
    return np.clip(top, 0.0, 1.0)


class FuzzyController:
    """Mamdani fuzzy controller with 40 bounded parameters."""

    def __init__(self, params: np.ndarray):
        """Initialize fuzzy controller.

        Parameters
        ----------
        params : numpy.ndarray
            Flat parameter vector of length 40.
        """
        vector = np.asarray(params, dtype=float).ravel()
        if vector.shape[0] != 40:
            raise ValueError(f'Expected params of length 40, got {vector.shape[0]}')
        self.params = np.clip(vector, 0.0, 1.0)
        self._y_grid = np.linspace(0.0, 1.0, 201)

    def _input_memberships(self, x: np.ndarray) -> np.ndarray:
        centers = self.params[8:20].reshape(4, 3)
        widths = np.maximum(self.params[20:32].reshape(4, 3), 1e-3)
        x_clip = np.clip(x, 0.0, 1.0)

        memberships = np.zeros((4, 3), dtype=float)
        for i in range(4):
            for j in range(3):
                c = centers[i, j]
                w = widths[i, j]
                memberships[i, j] = float(trimf(x_clip[i], c - w, c, c + w))
        return memberships

    def _output_sets(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        b = np.sort(self.params[0:8])
        low = trapmf(self._y_grid, 0.0, b[0], b[1], b[2])
        mid = trapmf(self._y_grid, b[1], b[2], b[3], b[4])
        high = trapmf(self._y_grid, b[4], b[5], b[6], 1.0)
        return low, mid, high

    def evaluate(self, features: np.ndarray) -> float:
        """Map a 4-element feature vector to similarity in [0, 1].

        Parameters
        ----------
        features : numpy.ndarray
            Input vector of shape ``(4,)``.

        Returns
        -------
        float
            Defuzzified similarity score.
        """
        x = np.asarray(features, dtype=float).ravel()
        if x.shape[0] != 4:
            raise ValueError(f'Expected 4 input features, got {x.shape[0]}')

        memberships = self._input_memberships(x)
        dominance = memberships.max(axis=1)
        top_inputs = np.argsort(dominance)[-2:]

        m_a = memberships[top_inputs[0]]
        m_b = memberships[top_inputs[1]]

        rule_scales = np.ones(9, dtype=float)
        rule_scales[:8] = np.clip(self.params[32:40], 0.0, 1.0)

        consequent_map = np.array(
            [
                [0, 0, 1],
                [0, 1, 2],
                [1, 2, 2],
            ]
        )

        term_activation = np.zeros(3, dtype=float)
        idx = 0
        for ia in range(3):
            for ib in range(3):
                firing = m_a[ia] * m_b[ib] * rule_scales[idx]
                term = consequent_map[ia, ib]
                term_activation[term] = max(term_activation[term], firing)
                idx += 1

        out_low, out_mid, out_high = self._output_sets()
        aggregated = np.maximum.reduce(
            [
                np.minimum(out_low, term_activation[0]),
                np.minimum(out_mid, term_activation[1]),
                np.minimum(out_high, term_activation[2]),
            ]
        )

        denom = float(np.sum(aggregated))
        if denom <= 1e-12:
            return 0.5
        score = float(np.sum(self._y_grid * aggregated) / denom)
        return float(np.clip(score, 0.0, 1.0))

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Vectorized prediction over an ``(N, 4)`` feature matrix.

        Parameters
        ----------
        X : numpy.ndarray
            Input matrix of shape ``(n_samples, 4)``.

        Returns
        -------
        numpy.ndarray
            Predicted similarity vector of shape ``(n_samples,)``.
        """
        array = np.asarray(X, dtype=float)
        if array.ndim != 2 or array.shape[1] != 4:
            raise ValueError(f'Expected X shape (n_samples, 4), got {array.shape}')
        return np.array([self.evaluate(row) for row in array], dtype=float)

"""Parameterized Mamdani fuzzy inference controller in pure NumPy.

The controller maps a 4-dimensional vector of neural similarity features to a
single semantic-similarity score in ``[0, 1]``. It implements a compact Mamdani
fuzzy inference system with three triangular membership functions per input,
three trapezoidal output sets, nine fuzzy rules whose firing strengths are
modulated by eight learnable scales (the ninth rule keeps a neutral weight of
``1.0``), and centroid defuzzification over a fixed 201-point output grid.

All 40 parameters are bounded to ``[0, 1]`` and learned by Differential
Evolution (see :mod:`neurofuzzy.optimizer`).

Notes
-----
This module is a clean NumPy re-implementation inspired by the Java/jFuzzyLogic
prototype that accompanied the original publication. It is *not* a bit-for-bit
port of the legacy FCL controller; instead it provides an equivalent, fully
reproducible, dependency-light controller with an explicit parameterization.
"""

from __future__ import annotations

import numpy as np

#: Number of input features consumed by the controller.
N_INPUTS = 4
#: Number of triangular terms per input.
N_TERMS = 3
#: Total number of learnable parameters.
N_PARAMS = 40
#: Number of points in the defuzzification grid.
GRID_SIZE = 201


def trimf(x, a: float, b: float, c: float) -> np.ndarray:
    """Compute triangular membership values in ``[0, 1]``.

    Parameters
    ----------
    x : numpy.ndarray or float
        Input value(s).
    a, b, c : float
        Left foot, peak, and right foot of the triangle.
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


def trapmf(x, a: float, b: float, c: float, d: float) -> np.ndarray:
    """Compute trapezoidal membership values in ``[0, 1]``.

    Parameters
    ----------
    x : numpy.ndarray or float
        Input value(s).
    a, b, c, d : float
        Left foot, left shoulder, right shoulder, and right foot.
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
    """Mamdani fuzzy controller with 40 bounded parameters.

    Parameters
    ----------
    params : numpy.ndarray
        Flat parameter vector of length 40. Values are clipped to ``[0, 1]``.

    Attributes
    ----------
    params : numpy.ndarray
        The clipped parameter vector actually used by the controller.

    Notes
    -----
    Output sets and the defuzzification grid are precomputed once at
    construction time. Batch :meth:`predict` is fully vectorized and produces
    results numerically identical to looping :meth:`evaluate`.
    """

    OUTPUT_BOUNDARY_SLICE = slice(0, 8)
    INPUT_CENTERS_SLICE = slice(8, 20)
    INPUT_WIDTHS_SLICE = slice(20, 32)
    RULE_SCALES_SLICE = slice(32, 40)  # 8 learned scales; 9th rule keeps scale 1.0

    # Mamdani rule base: for the two most-activated inputs, the (term_a, term_b)
    # pair maps to an output term index (0=low, 1=mid, 2=high).
    CONSEQUENT_MAP = np.array(
        [
            [0, 0, 1],
            [0, 1, 2],
            [1, 2, 2],
        ]
    )

    def __init__(self, params: np.ndarray):
        vector = np.asarray(params, dtype=float).ravel()
        if vector.shape[0] != N_PARAMS:
            raise ValueError(f"Expected params of length {N_PARAMS}, got {vector.shape[0]}")
        self.params = np.clip(vector, 0.0, 1.0)
        self._y_grid = np.linspace(0.0, 1.0, GRID_SIZE)

        self._centers = self.params[self.INPUT_CENTERS_SLICE].reshape(N_INPUTS, N_TERMS)
        self._widths = np.maximum(
            self.params[self.INPUT_WIDTHS_SLICE].reshape(N_INPUTS, N_TERMS), 1e-3
        )
        self._rule_scales = np.ones(9, dtype=float)
        self._rule_scales[:8] = np.clip(self.params[self.RULE_SCALES_SLICE], 0.0, 1.0)

        low, mid, high = self._build_output_sets()
        self._out_low, self._out_mid, self._out_high = low, mid, high
        self._out_stack = np.stack([low, mid, high])  # (3, GRID_SIZE)

    def _build_output_sets(self):
        b = np.sort(self.params[self.OUTPUT_BOUNDARY_SLICE])
        low = trapmf(self._y_grid, 0.0, b[0], b[1], b[2])
        mid = trapmf(self._y_grid, b[1], b[2], b[3], b[4])
        high = trapmf(self._y_grid, b[4], b[5], b[6], 1.0)
        return low, mid, high

    def _output_sets(self):
        """Return the precomputed ``(low, mid, high)`` output membership arrays."""
        return self._out_low, self._out_mid, self._out_high

    def _memberships(self, x: np.ndarray) -> np.ndarray:
        """Triangular membership values for one sample, shape ``(4, 3)``."""
        x_clip = np.clip(x, 0.0, 1.0).reshape(N_INPUTS, 1)
        left = (x_clip - (self._centers - self._widths)) / self._widths
        right = ((self._centers + self._widths) - x_clip) / self._widths
        return np.clip(np.minimum(left, right), 0.0, 1.0)

    def _input_memberships(self, x: np.ndarray) -> np.ndarray:
        """Backward-compatible alias for :meth:`_memberships`."""
        return self._memberships(np.asarray(x, dtype=float).ravel())

    def _all_memberships(self, X: np.ndarray) -> np.ndarray:
        """Vectorized memberships for a batch, shape ``(N, 4, 3)``."""
        x_clip = np.clip(X, 0.0, 1.0)[:, :, None]
        centers = self._centers[None]
        widths = self._widths[None]
        left = (x_clip - (centers - widths)) / widths
        right = ((centers + widths) - x_clip) / widths
        return np.clip(np.minimum(left, right), 0.0, 1.0)

    def _term_activations(self, memberships: np.ndarray) -> np.ndarray:
        """Aggregate rule firings into output-term activations, shape ``(3,)``."""
        dominance = memberships.max(axis=1)
        top_inputs = np.argsort(dominance)[-2:]
        m_a = memberships[top_inputs[0]]
        m_b = memberships[top_inputs[1]]

        term_activation = np.zeros(N_TERMS, dtype=float)
        idx = 0
        for ia in range(N_TERMS):
            for ib in range(N_TERMS):
                firing = m_a[ia] * m_b[ib] * self._rule_scales[idx]
                term = self.CONSEQUENT_MAP[ia, ib]
                term_activation[term] = max(term_activation[term], firing)
                idx += 1
        return term_activation

    def _term_activations_batch(self, all_mem: np.ndarray) -> np.ndarray:
        """Vectorized rule aggregation for a batch, shape ``(N, 3)``.

        Identical results to applying :meth:`_term_activations` to each row.
        """
        n = all_mem.shape[0]
        dominance = all_mem.max(axis=2)
        top = np.argsort(dominance, axis=1)[:, -2:]
        rows = np.arange(n)
        m_a = all_mem[rows, top[:, 0]]
        m_b = all_mem[rows, top[:, 1]]

        scales = self._rule_scales.reshape(N_TERMS, N_TERMS)  # idx = ia*3 + ib
        firings = m_a[:, :, None] * m_b[:, None, :] * scales[None]  # (N, 3, 3)

        term_acts = np.zeros((n, N_TERMS), dtype=float)
        for term in range(N_TERMS):
            mask = term == self.CONSEQUENT_MAP
            term_acts[:, term] = firings[:, mask].max(axis=1)
        return term_acts

    def _defuzzify(self, term_activation: np.ndarray) -> float:
        """Centroid defuzzification for a single ``(3,)`` activation vector."""
        aggregated = np.maximum.reduce(
            [
                np.minimum(self._out_low, term_activation[0]),
                np.minimum(self._out_mid, term_activation[1]),
                np.minimum(self._out_high, term_activation[2]),
            ]
        )
        denom = float(np.sum(aggregated))
        if denom <= 1e-12:
            return 0.5
        score = float(np.sum(self._y_grid * aggregated) / denom)
        return float(np.clip(score, 0.0, 1.0))

    def evaluate(self, features: np.ndarray) -> float:
        """Map a 4-element feature vector to similarity in ``[0, 1]``."""
        x = np.asarray(features, dtype=float).ravel()
        if x.shape[0] != N_INPUTS:
            raise ValueError(f"Expected {N_INPUTS} input features, got {x.shape[0]}")
        return self._defuzzify(self._term_activations(self._memberships(x)))

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Vectorized prediction over an ``(N, 4)`` feature matrix.

        Produces results numerically identical to calling :meth:`evaluate` on
        each row, but precomputes static fuzzy sets once and vectorizes the
        membership, rule-aggregation, and centroid steps across samples.
        """
        array = np.asarray(X, dtype=float)
        if array.ndim != 2 or array.shape[1] != N_INPUTS:
            raise ValueError(f"Expected X shape (n_samples, {N_INPUTS}), got {array.shape}")
        if array.shape[0] == 0:
            return np.empty(0, dtype=float)

        all_mem = self._all_memberships(array)
        term_acts = self._term_activations_batch(all_mem)

        clipped = np.minimum(self._out_stack[None, :, :], term_acts[:, :, None])
        aggregated = clipped.max(axis=1)
        denom = np.sum(aggregated, axis=1)
        numer = np.sum(self._y_grid * aggregated, axis=1)
        scores = np.full(array.shape[0], 0.5, dtype=float)
        valid = denom > 1e-12
        scores[valid] = numer[valid] / denom[valid]
        return np.clip(scores, 0.0, 1.0)

"""Tests for the vectorized fuzzy controller."""

import numpy as np
import pytest

from neurofuzzy.fuzzy_controller import N_PARAMS, FuzzyController, trapmf, trimf


def test_param_length_validation():
    with pytest.raises(ValueError):
        FuzzyController(np.zeros(N_PARAMS - 1))


def test_output_range():
    rng = np.random.default_rng(0)
    for _ in range(20):
        c = FuzzyController(rng.random(N_PARAMS))
        out = c.predict(rng.random((10, 4)))
        assert out.shape == (10,)
        assert np.all(out >= 0.0) and np.all(out <= 1.0)


def test_batch_matches_scalar_bit_exact():
    rng = np.random.default_rng(1)
    for _ in range(50):
        c = FuzzyController(rng.random(N_PARAMS))
        X = rng.random((rng.integers(1, 15), 4))
        batch = c.predict(X)
        scalar = np.array([c.evaluate(row) for row in X])
        assert np.array_equal(batch, scalar)


def test_empty_input_returns_empty():
    c = FuzzyController(np.full(N_PARAMS, 0.5))
    assert c.predict(np.empty((0, 4))).shape == (0,)


def test_evaluate_wrong_shape_raises():
    c = FuzzyController(np.full(N_PARAMS, 0.5))
    with pytest.raises(ValueError):
        c.evaluate(np.array([0.1, 0.2, 0.3]))


def test_trimf_and_trapmf_bounds():
    grid = np.linspace(-1, 2, 50)
    assert np.all((trimf(grid, 0, 0.5, 1) >= 0) & (trimf(grid, 0, 0.5, 1) <= 1))
    assert np.all((trapmf(grid, 0, 0.2, 0.8, 1) >= 0) & (trapmf(grid, 0, 0.2, 0.8, 1) <= 1))


def test_params_are_clipped():
    c = FuzzyController(np.full(N_PARAMS, 5.0))
    assert np.all(c.params <= 1.0) and np.all(c.params >= 0.0)

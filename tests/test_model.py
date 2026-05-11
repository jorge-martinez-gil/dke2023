"""Unit tests for neurofuzzy model and controller."""

import numpy as np

from neurofuzzy.fuzzy_controller import FuzzyController
from neurofuzzy.model import NeurofuzzyModel


def test_model_instantiation() -> None:
    model = NeurofuzzyModel()
    assert model.n_features == 4


def test_fuzzy_controller_output_range() -> None:
    controller = FuzzyController(np.full(40, 0.5))
    value = controller.evaluate(np.array([0.2, 0.4, 0.6, 0.8]))
    assert 0.0 <= value <= 1.0


def test_fuzzy_controller_predict_matches_scalar() -> None:
    controller = FuzzyController(np.full(40, 0.5))
    X = np.array([[0.1, 0.2, 0.3, 0.4], [0.9, 0.8, 0.7, 0.6]])
    y_vec = controller.predict(X)
    y_scalar = np.array([controller.evaluate(row) for row in X])
    assert np.allclose(y_vec, y_scalar)

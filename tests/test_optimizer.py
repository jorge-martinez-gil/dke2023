"""Tests for the Differential Evolution optimizer and model plumbing."""

import numpy as np
import pytest

from neurofuzzy.data_loader import load_dataset
from neurofuzzy.model import NeurofuzzyModel
from neurofuzzy.optimizer import optimize

ROOT = __import__("pathlib").Path(__file__).resolve().parents[1]


def _small_data():
    X, y = load_dataset(ROOT / "datasets" / "mc.txt", n_features=4)
    return X[:12], y[:12]


def test_optimize_returns_40_params_and_is_deterministic():
    X, y = _small_data()
    r1 = optimize(X, y, maxiter=3, popsize=3, seed=0)
    r2 = optimize(X, y, maxiter=3, popsize=3, seed=0)
    assert r1.x.shape == (40,)
    assert np.array_equal(r1.x, r2.x)  # same seed -> identical result


def test_model_fit_predict_roundtrip(tmp_path):
    X, y = _small_data()
    model = NeurofuzzyModel(maxiter=2, popsize=3, train_ratio=1.0, random_state=0, verbose=False)
    model.fit(X, y)
    pred = model.predict(X)
    assert pred.shape == (X.shape[0],)
    assert np.all((pred >= 0) & (pred <= 1))

    path = tmp_path / "params.npy"
    model.save(path)
    reloaded = NeurofuzzyModel.load(path)
    assert np.array_equal(reloaded.predict(X), pred)


def test_predict_before_fit_raises():
    with pytest.raises(RuntimeError):
        NeurofuzzyModel().predict(np.zeros((2, 4)))


def test_fit_wrong_feature_count_raises():
    model = NeurofuzzyModel()
    with pytest.raises(ValueError):
        model.fit(np.zeros((5, 3)), np.zeros(5))

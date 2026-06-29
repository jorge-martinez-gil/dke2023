"""Tests for reference baselines."""

import numpy as np
import pytest

from neurofuzzy.baselines import (
    BASELINES,
    BestSingleFeatureBaseline,
    LinearRegressionBaseline,
    MeanFeatureBaseline,
    get_baseline,
)


def _toy():
    rng = np.random.default_rng(0)
    X = rng.random((40, 4))
    y = X[:, 2] * 0.8 + 0.1  # feature 2 is most informative
    return X, y


def test_registry_and_factory():
    assert set(BASELINES) == {"mean-of-features", "best-single-feature", "linear-regression"}
    for name in BASELINES:
        assert get_baseline(name).name == name


def test_unknown_baseline_raises():
    with pytest.raises(KeyError):
        get_baseline("does-not-exist")


def test_outputs_in_range_and_shape():
    X, y = _toy()
    for name in BASELINES:
        pred = get_baseline(name).fit(X, y).predict(X)
        assert pred.shape == (X.shape[0],)
        assert np.all(pred >= 0.0) and np.all(pred <= 1.0)


def test_best_single_feature_selects_informative_column():
    X, y = _toy()
    model = BestSingleFeatureBaseline().fit(X, y)
    assert model.feature_index_ == 2


def test_mean_of_features_is_row_mean():
    X = np.array([[0.2, 0.4, 0.6, 0.8]])
    assert np.isclose(MeanFeatureBaseline().predict(X)[0], 0.5)


def test_linear_regression_recovers_linear_signal():
    X, y = _toy()
    model = LinearRegressionBaseline().fit(X, y)
    # On (near) linear data the fit should be strongly correlated.
    pred = model.predict(X)
    assert np.corrcoef(pred, y)[0, 1] > 0.95


def test_predict_before_fit_raises():
    with pytest.raises(RuntimeError):
        BestSingleFeatureBaseline().predict(np.zeros((2, 4)))
    with pytest.raises(RuntimeError):
        LinearRegressionBaseline().predict(np.zeros((2, 4)))

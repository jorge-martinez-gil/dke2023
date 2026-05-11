"""Unit tests for neurofuzzy.metrics."""

import numpy as np

from neurofuzzy.metrics import cosine_similarity, evaluate_all, pearson, spearman


def test_pearson_perfect_correlation() -> None:
    y_true = np.array([1.0, 2.0, 3.0])
    y_pred = np.array([2.0, 4.0, 6.0])
    assert np.isclose(pearson(y_true, y_pred), 1.0)


def test_pearson_anti_correlation() -> None:
    y_true = np.array([1.0, 2.0, 3.0])
    y_pred = np.array([3.0, 2.0, 1.0])
    assert np.isclose(pearson(y_true, y_pred), -1.0)


def test_spearman_known_rankings() -> None:
    y_true = np.array([10.0, 20.0, 30.0, 40.0])
    y_pred = np.array([1.0, 2.0, 3.0, 4.0])
    assert np.isclose(spearman(y_true, y_pred), 1.0)


def test_cosine_similarity_identical_vectors() -> None:
    y_true = np.array([1.0, 2.0, 3.0])
    assert np.isclose(cosine_similarity(y_true, y_true), 1.0)


def test_evaluate_all_keys() -> None:
    y_true = np.array([0.1, 0.2, 0.3])
    y_pred = np.array([0.1, 0.2, 0.3])
    result = evaluate_all(y_true, y_pred)
    assert set(result.keys()) == {
        'pearson',
        'spearman',
        'cosine_similarity',
        'mean_absolute_error',
        'root_mean_squared_error',
    }

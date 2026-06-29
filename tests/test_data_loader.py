"""Unit tests for neurofuzzy.data_loader."""

from pathlib import Path

import numpy as np

from neurofuzzy.data_loader import load_dataset, train_test_split

ROOT = Path(__file__).resolve().parents[1]


def test_load_mc_shape() -> None:
    X, y = load_dataset(ROOT / "datasets" / "mc.txt", n_features=4)
    assert X.shape == (30, 4)
    assert y.shape == (30,)


def test_load_geresid_shape() -> None:
    X, y = load_dataset(ROOT / "datasets" / "geresid.txt", n_features=4)
    assert X.shape == (50, 4)
    assert y.shape == (50,)


def test_train_test_split_ratio() -> None:
    X = np.arange(40, dtype=float).reshape(10, 4)
    y = np.arange(10, dtype=float)
    X_train, X_test, y_train, y_test = train_test_split(X, y, train_ratio=0.6, random_state=0)
    assert X_train.shape[0] == 6
    assert X_test.shape[0] == 4
    assert y_train.shape[0] == 6
    assert y_test.shape[0] == 4


def test_malformed_row_raises(tmp_path):
    import pytest

    bad = tmp_path / "bad.txt"
    bad.write_text("1.0, 0.5, 0.4\n")  # too few columns for n_features=4
    with pytest.raises(ValueError):
        load_dataset(bad, n_features=4)


def test_train_test_split_invalid_ratio_raises():
    import pytest

    X = np.zeros((10, 4))
    y = np.zeros(10)
    with pytest.raises(ValueError):
        train_test_split(X, y, train_ratio=1.5)

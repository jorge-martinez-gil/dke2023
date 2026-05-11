"""High-level model API for neurofuzzy semantic similarity."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .data_loader import train_test_split
from .fuzzy_controller import FuzzyController
from .metrics import evaluate_all
from .optimizer import optimize


class NeurofuzzyModel:
    """Neurofuzzy Semantic Similarity Model.

    Combines neural-derived features with an evolutionary-optimized
    fuzzy inference system for word/sentence similarity measurement.

    Parameters
    ----------
    n_features : int, default=4
        Number of input features (neural similarity scores).
    train_ratio : float, default=0.6
        Fraction of data used for training during fit().
    maxiter : int, default=200
        Maximum DE optimization iterations.
    popsize : int, default=15
        DE population size multiplier.
    random_state : int or None, default=42
        Seed used for splitting and optimization.
    verbose : bool, default=True
        Whether to print fit status messages.
    """

    def __init__(
        self,
        n_features: int = 4,
        train_ratio: float = 0.6,
        maxiter: int = 200,
        popsize: int = 15,
        random_state: int | None = 42,
        verbose: bool = True,
    ) -> None:
        self.n_features = n_features
        self.train_ratio = train_ratio
        self.maxiter = maxiter
        self.popsize = popsize
        self.random_state = random_state
        self.verbose = verbose
        self.params_: np.ndarray | None = None
        self.controller_: FuzzyController | None = None

    def fit(self, X: np.ndarray, y: np.ndarray) -> 'NeurofuzzyModel':
        """Fit the model by optimizing controller parameters.

        Parameters
        ----------
        X : numpy.ndarray
            Input feature matrix.
        y : numpy.ndarray
            Ground-truth similarity targets.

        Returns
        -------
        NeurofuzzyModel
            Fitted model instance.
        """
        X_arr = np.asarray(X, dtype=float)
        y_arr = np.asarray(y, dtype=float)
        if X_arr.ndim != 2 or X_arr.shape[1] != self.n_features:
            raise ValueError(
                f'Expected X shape (n_samples, {self.n_features}), got {X_arr.shape}'
            )

        X_train = X_arr
        y_train = y_arr
        if self.train_ratio < 1.0:
            X_train, _, y_train, _ = train_test_split(
                X_arr,
                y_arr,
                train_ratio=self.train_ratio,
                random_state=self.random_state,
            )

        result = optimize(
            X_train,
            y_train,
            n_params=40,
            maxiter=self.maxiter,
            popsize=self.popsize,
            seed=42 if self.random_state is None else self.random_state,
        )
        self.params_ = np.asarray(result.x, dtype=float)
        self.controller_ = FuzzyController(self.params_)
        if self.verbose:
            print(f'Optimization finished: best objective = {result.fun:.6f}')
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Predict similarity scores.

        Parameters
        ----------
        X : numpy.ndarray
            Input feature matrix.

        Returns
        -------
        numpy.ndarray
            Predicted similarity scores.
        """
        if self.controller_ is None:
            raise RuntimeError('Model is not fitted. Call fit() or load() first.')
        return self.controller_.predict(np.asarray(X, dtype=float))

    def evaluate(self, X: np.ndarray, y: np.ndarray) -> dict[str, float]:
        """Evaluate model predictions against ground truth.

        Parameters
        ----------
        X : numpy.ndarray
            Input feature matrix.
        y : numpy.ndarray
            Ground-truth target vector.

        Returns
        -------
        dict
            Full metrics dictionary.
        """
        y_pred = self.predict(X)
        return evaluate_all(np.asarray(y, dtype=float), y_pred)

    def save(self, path: str | Path) -> None:
        """Save optimized parameter vector as ``.npy``.

        Parameters
        ----------
        path : str or pathlib.Path
            Destination file path.
        """
        if self.params_ is None:
            raise RuntimeError('Model has no learned parameters to save.')
        np.save(path, self.params_)

    @classmethod
    def load(cls, path: str | Path) -> 'NeurofuzzyModel':
        """Load model parameters from ``.npy`` file.

        Parameters
        ----------
        path : str or pathlib.Path
            Parameter file path.

        Returns
        -------
        NeurofuzzyModel
            Instantiated model with loaded parameters.
        """
        model = cls()
        model.params_ = np.load(path)
        model.controller_ = FuzzyController(model.params_)
        return model

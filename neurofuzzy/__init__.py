"""neurofuzzy -- Neurofuzzy Semantic Similarity Measurement.

A clean, reproducible Python implementation of a fuzzy-inference approach to
semantic similarity measurement, with a unified benchmark harness, reference
baselines, and publication-quality reporting.

Reference
---------
Martinez-Gil, J., Mokadem, R., Kueng, J., & Hameurlain, A. (2023).
Neurofuzzy semantic similarity measurement.
*Data & Knowledge Engineering*, 145, 102155.
https://doi.org/10.1016/j.datak.2023.102155
"""

from __future__ import annotations

__version__ = "1.1.0"

from .baselines import (
    BASELINES,
    BestSingleFeatureBaseline,
    LinearRegressionBaseline,
    MeanFeatureBaseline,
    get_baseline,
)
from .fuzzy_controller import FuzzyController, trapmf, trimf
from .metrics import (
    cosine_similarity,
    evaluate_all,
    mean_absolute_error,
    pearson,
    root_mean_squared_error,
    spearman,
)
from .model import NeurofuzzyModel

__all__ = [
    "__version__",
    "NeurofuzzyModel",
    "FuzzyController",
    "trimf",
    "trapmf",
    "pearson",
    "spearman",
    "cosine_similarity",
    "mean_absolute_error",
    "root_mean_squared_error",
    "evaluate_all",
    "get_baseline",
    "BASELINES",
    "MeanFeatureBaseline",
    "BestSingleFeatureBaseline",
    "LinearRegressionBaseline",
]

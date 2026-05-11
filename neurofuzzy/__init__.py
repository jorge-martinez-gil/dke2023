"""neurofuzzy — Neurofuzzy Semantic Similarity Measurement.

Reference implementation of:
  Martinez-Gil et al. (2023). Neurofuzzy semantic similarity measurement.
  Data & Knowledge Engineering, 145, 102155.
  https://doi.org/10.1016/j.datak.2023.102155
"""

from .metrics import cosine_similarity, pearson, spearman
from .model import NeurofuzzyModel

__version__ = "1.0.0"
__all__ = ["NeurofuzzyModel", "pearson", "spearman", "cosine_similarity"]

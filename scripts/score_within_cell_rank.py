from __future__ import annotations

import numpy as np
from scipy.stats import rankdata


def score(matrix, program_columns=None):
    values = matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)
    if values.shape[1] == 1:
        return np.ones(values.shape[0])
    ranks = np.apply_along_axis(rankdata, 1, values, method="average")
    percentile = (ranks - 1.0) / (values.shape[1] - 1.0)
    if program_columns is not None:
        percentile = percentile[:, list(program_columns)]
    return np.mean(percentile, axis=1)

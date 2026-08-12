from __future__ import annotations

import numpy as np
from scipy import sparse


def score(matrix):
    if sparse.issparse(matrix):
        return np.asarray((matrix > 0).mean(axis=1)).ravel()
    return np.mean(np.asarray(matrix) > 0, axis=1)


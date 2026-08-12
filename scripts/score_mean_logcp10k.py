from __future__ import annotations

import numpy as np
from scipy import sparse


def counts_to_logcp10k(matrix):
    source = matrix.astype(np.float64)
    totals = np.asarray(source.sum(axis=1)).ravel()
    scale = np.divide(10000.0, totals, out=np.zeros_like(totals), where=totals > 0)
    if sparse.issparse(source):
        result = sparse.diags(scale) @ source
        result.data = np.log1p(result.data)
        return result
    return np.log1p(source * scale[:, None])


def score(matrix):
    return np.asarray(matrix.mean(axis=1)).ravel()


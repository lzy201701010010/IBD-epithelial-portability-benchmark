from __future__ import annotations

import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.genmod.cov_struct import Exchangeable


def fit_partition_adjusted_gee(frame, score_column: str, group_column: str, disease_column: str, partition_column: str):
    """Fit the frozen Gaussian identity-link GEE to an author-supplied donor-partition frame."""
    data = frame[[score_column, group_column, disease_column, partition_column]].copy()
    if data.empty or not np.isfinite(data[score_column].to_numpy(dtype=float)).all():
        raise ValueError("Scores must be present and finite")
    formula = f"{score_column} ~ C({disease_column}, Treatment(reference='NON_IBD_CONTROL')) + C({partition_column})"
    model = smf.gee(
        formula,
        groups=group_column,
        data=data,
        family=sm.families.Gaussian(),
        cov_struct=Exchangeable(),
    )
    if np.linalg.matrix_rank(model.exog) < model.exog.shape[1]:
        raise ValueError("Design matrix is not full rank")
    return model.fit(maxiter=100, cov_type="robust")

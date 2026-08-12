from __future__ import annotations

import numpy as np


def sample_estimability(relevant_cells: int, epithelial_cells: int, mapping_fraction: float, mapped_genes: int) -> tuple[bool, str]:
    if epithelial_cells < 50:
        return False, "LT50_EPITHELIAL_CELLS"
    if relevant_cells < 20:
        return False, "LT20_RELEVANT_LINEAGE_CELLS"
    if mapping_fraction < 0.8:
        return False, "LT80_PERCENT_GENE_COVERAGE"
    if mapped_genes < 3:
        return False, "LT3_MAPPED_GENES"
    return True, ""


def inference_estimability(group_a, group_b) -> tuple[bool, str]:
    if len(group_a) < 8 or len(group_b) < 8:
        return False, "LT8_DONORS_PER_GROUP"
    values = np.asarray(list(group_a) + list(group_b), dtype=float)
    if not np.isfinite(values).all():
        return False, "NON_FINITE_DONOR_SCORE"
    if np.var(group_a, ddof=1) == 0 or np.var(group_b, ddof=1) == 0:
        return False, "ZERO_VARIANCE"
    return True, ""


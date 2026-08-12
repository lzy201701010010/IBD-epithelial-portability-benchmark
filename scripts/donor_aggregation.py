from __future__ import annotations

import numpy as np


def build_sample_to_donor(pairs) -> dict[str, str]:
    result: dict[str, str] = {}
    for sample, donor in pairs:
        if not sample or not donor:
            raise ValueError("MISSING_SAMPLE_OR_DONOR_ID")
        if sample in result and result[sample] != donor:
            raise ValueError("AMBIGUOUS_SAMPLE_TO_DONOR")
        result[sample] = donor
    return result


def aggregate_sample_scores(sample_to_donor: dict[str, str], sample_scores: dict[str, float]) -> dict[str, float]:
    grouped: dict[str, list[float]] = {}
    for sample, value in sample_scores.items():
        if sample not in sample_to_donor or not sample_to_donor[sample]:
            raise ValueError("MISSING_DONOR_ID")
        grouped.setdefault(sample_to_donor[sample], []).append(float(value))
    return {donor: float(np.mean(values)) for donor, values in grouped.items()}

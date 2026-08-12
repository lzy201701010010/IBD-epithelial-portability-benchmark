from __future__ import annotations

import math

import numpy as np
from scipy import optimize, stats


def hedges_g(group_a, group_b) -> float:
    a, b = np.asarray(group_a, dtype=float), np.asarray(group_b, dtype=float)
    pooled = math.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    d = (a.mean() - b.mean()) / pooled
    correction = 1.0 - 3.0 / (4.0 * (len(a) + len(b)) - 9.0)
    return float(correction * d)


def mean_difference(group_a, group_b) -> float:
    return float(np.mean(group_a) - np.mean(group_b))


def hedges_g_ci(group_a, group_b, alpha: float = 0.05) -> tuple[float, float]:
    a, b = np.asarray(group_a, dtype=float), np.asarray(group_b, dtype=float)
    n1, n0 = len(a), len(b)
    pooled = math.sqrt(((n1 - 1) * a.var(ddof=1) + (n0 - 1) * b.var(ddof=1)) / (n1 + n0 - 2))
    t_observed = (a.mean() - b.mean()) / (pooled * math.sqrt(1 / n1 + 1 / n0))
    df = n1 + n0 - 2
    correction = 1.0 - 3.0 / (4.0 * (n1 + n0) - 9.0)

    def solve(target: float) -> float:
        return optimize.brentq(lambda ncp: stats.nct.cdf(t_observed, df, ncp) - target, -100.0, 100.0)

    lower_ncp = solve(1 - alpha / 2)
    upper_ncp = solve(alpha / 2)
    scale = math.sqrt(1 / n1 + 1 / n0)
    return float(correction * lower_ncp * scale), float(correction * upper_ncp * scale)


def welch_ci(group_a, group_b, alpha: float = 0.05) -> tuple[float, float]:
    a, b = np.asarray(group_a, dtype=float), np.asarray(group_b, dtype=float)
    diff = a.mean() - b.mean()
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    se = math.sqrt(va + vb)
    df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    critical = stats.t.ppf(1 - alpha / 2, df)
    return float(diff - critical * se), float(diff + critical * se)

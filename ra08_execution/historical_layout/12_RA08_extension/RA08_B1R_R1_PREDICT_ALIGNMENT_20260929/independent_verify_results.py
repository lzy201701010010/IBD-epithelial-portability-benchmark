import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE / "_deps"))
import h5py
import numpy as np
from scipy import stats

METHODS = [f"SP{i:02d}_{method}" for i in range(1, 12) if i != 8 for method in ("SM01", "SM02", "SM03")]
donor_data = {}
checks = {}
for cohort in ("PREDICT_FGID", "PREDICT_CD", "ELMENTAITE"):
    with np.load(HERE / f"RA08_B1R_R1_{cohort}_EPITHELIAL_CELL_SCORES.npz", allow_pickle=False) as z:
        donor = z["donor"].astype(str)
        group = z["group"].astype(str)
        sample = z["sample"].astype(str)
        assert set(METHODS) <= set(z.files)
        assert len(donor) == len(group) == len(sample)
        sample_ids = defaultdict(set)
        for d, s in zip(donor, sample):
            sample_ids[d].add(s)
        assert all(len(v) == 1 for v in sample_ids.values())
        with (HERE / f"RA08_B1R_R1_{cohort}_DONOR_SCORES.tsv").open(encoding="utf-8", newline="") as f:
            frozen = {r["donor"]: r for r in csv.DictReader(f, delimiter="\t")}
        assert set(frozen) == set(donor)
        comparisons = 0
        for d in frozen:
            mask = donor == d
            assert len(set(group[mask])) == 1
            for member in METHODS:
                actual = float(frozen[d][member])
                expected = float(np.mean(z[member][mask]))
                assert abs(actual - expected) <= 1e-12 * max(1.0, abs(expected))
                comparisons += 1
                donor_data[(d, frozen[d]["group"], cohort, member)] = actual
        checks[cohort] = {"epithelial_cells": len(donor), "donors": len(frozen), "donor_method_mean_checks": comparisons}

with (HERE / "RA08_B1R_R1_EXTERNAL_EVIDENCE_MATRIX.tsv").open(encoding="utf-8", newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
assert len(rows) == 33
numeric = [r for r in rows if r["predict_status"] == "ESTIMABLE"]
assert len(numeric) == 30
p = []
for row in numeric:
    member = f"{row['program']}_{row['method']}"
    cd = np.array([v for (d, g, c, m), v in donor_data.items() if c == "PREDICT_CD" and g == "CD" and m == member])
    fg = np.array([v for (d, g, c, m), v in donor_data.items() if c == "PREDICT_FGID" and g == "FGID" and m == member])
    ec = np.array([v for (d, g, c, m), v in donor_data.items() if c == "ELMENTAITE" and g == "CD" and m == member])
    eh = np.array([v for (d, g, c, m), v in donor_data.items() if c == "ELMENTAITE" and g == "HEALTHY_NON_IBD" and m == member])
    assert len(cd) == 14 and len(fg) == 13 and len(ec) == 7 and len(eh) == 8
    assert abs(float(row["predict_CD_minus_FGID"]) - (cd.mean() - fg.mean())) < 1e-12
    assert abs(float(row["elmentaite_CD_minus_control"]) - (ec.mean() - eh.mean())) < 1e-12
    value = float(stats.ttest_ind(cd, fg, equal_var=False, alternative="two-sided").pvalue)
    assert abs(float(row["P"]) - value) <= 1e-12 * max(1.0, value)
    assert row["elmentaite_status"] == "SUPPORTIVE_DESCRIPTIVE_ONLY_NOT_FORMAL"
    p.append(value)
order = np.argsort(p, kind="stable")
q = np.ones(len(p))
running = 1.0
for pos in range(len(p) - 1, -1, -1):
    idx = int(order[pos])
    running = min(running, p[idx] * len(p) / (pos + 1))
    q[idx] = min(1.0, running)
for row, value in zip(numeric, q):
    assert abs(float(row["Q"]) - value) < 1e-12
    assert int(row["BH_numeric_denominator"]) == 30

h5ad = ROOT / "12_RA08_extension/RA08_B1A_CONTINUATION_20260928/raw/final_pediatric_object_cellxgene_c.h5ad"
explicit_zero = 0
with h5py.File(h5ad, "r") as h:
    data = h["raw/X/data"]
    for start in range(0, len(data), 5_000_000):
        explicit_zero += int(np.count_nonzero(data[start:start + 5_000_000] == 0))
assert explicit_zero == 0
result = {"status": "PASS", "cohort_score_checks": checks, "formal_P_recomputed": len(numeric), "single_BH_Q_recomputed": len(numeric), "Elmentaite_raw_X_explicit_zero_coordinates": explicit_zero, "SP08_terminal_rows": sum(r["program"] == "SP08" for r in rows)}
(HERE / "RA08_B1R_R1_INDEPENDENT_RESULT_VERIFICATION.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))

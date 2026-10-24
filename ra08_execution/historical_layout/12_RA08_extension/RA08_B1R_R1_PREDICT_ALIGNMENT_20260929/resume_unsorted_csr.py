"""Technical continuation after the source H5AD's unsorted, unique CSR indices.

PREDICT score packages are reused byte for byte. Only the H5AD index-order
assumption is removed; all frozen score equations come from score_and_analyze.
"""

import csv
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE / "_deps"))
sys.path.insert(0, str(HERE))
import h5py
import numpy as np
import score_and_analyze as frozen


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def load_predict_donors():
    scores = {}
    for group in ("FGID", "CD"):
        path = HERE / f"RA08_B1R_R1_PREDICT_{group}_DONOR_SCORES.tsv"
        with path.open(encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                key = (row["donor"], row["group"])
                if key in scores:
                    raise ValueError("duplicate frozen PREDICT donor")
                scores[key] = {(program, method): float(row[f"{program}_{method}"]) for program in frozen.PROGRAMS for method in ("SM01", "SM02", "SM03")}
    return scores


def elmentaite_unsorted_csr():
    b1a = ROOT / "12_RA08_extension/RA08_B1A_CONTINUATION_20260928"
    with (b1a / "ELMENTAITE_CELL_SAMPLE_PARTICIPANT_BRIDGE.tsv").open(encoding="utf-8-sig", newline="") as f:
        bridge = {r["source_cell_id"]: r for r in csv.DictReader(f, delimiter="\t")}
    cell_ids, samples, donors, groups, libraries = [], [], [], [], []
    target_rows, positive_rows = [], []
    with h5py.File(b1a / "raw/final_pediatric_object_cellxgene_c.h5ad", "r") as h:
        obs_ids = [v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in h["obs/index"][:]]
        symbols = [v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in h["raw/var/index"][:]]
        target_columns = [symbols.index(gene) for gene in frozen.TARGET_GENES]
        if len(set(symbols)) != len(symbols) or len(obs_ids) != len(bridge):
            raise ValueError("Elmentaite frozen identity drift")
        x = h["raw/X"]
        offsets = x["indptr"][:]
        for i, cell in enumerate(obs_ids):
            record = bridge[cell]
            if record["epithelial_disposition"] != "INCLUDE_EPITHELIAL":
                continue
            start, stop = int(offsets[i]), int(offsets[i + 1])
            indices = x["indices"][start:stop]
            data = x["data"][start:stop]
            if np.any(data < 0) or not np.allclose(data, np.round(data), rtol=0, atol=0):
                raise ValueError("Elmentaite raw count semantics drift")
            if len(np.unique(indices)) != len(indices):
                raise ValueError("Elmentaite duplicate coordinates")
            lookup = dict(zip(indices.tolist(), data.tolist()))
            target_rows.append([lookup.get(column, 0) for column in target_columns])
            positive_rows.append(np.asarray(data, dtype=np.float64))
            cell_ids.append(cell)
            samples.append(record["source_sample_name"])
            donors.append(record["participant_id"])
            if record["ra08_group"] == "CROHN_DISEASE":
                groups.append("CD")
            elif record["ra08_group"] == "HEALTHY_NON_IBD_CONTROL":
                groups.append("HEALTHY_NON_IBD")
            else:
                raise ValueError("unknown frozen Elmentaite group")
            libraries.append(float(np.sum(data, dtype=np.float64)))
    target = np.asarray(target_rows, dtype=np.float64)
    rank_target = np.empty_like(target, dtype=np.float64)
    for i, positive in enumerate(positive_rows):
        rank_target[i] = frozen.ranks_for_cell(positive, target[i], len(symbols))
    scores = frozen.score_arrays(target, np.asarray(libraries), rank_target)
    return cell_ids, samples, donors, groups, scores


def main():
    gate = json.loads((HERE / "RA08_B1R_R1_TECHNICAL_RECOVERY_GATE.json").read_text(encoding="utf-8"))
    if gate["decision"] != "PASS" or sha(Path(__file__)) != gate["recovery_source_sha256"] or sha(HERE / "score_and_analyze.py") != gate["original_scoring_source_sha256"]:
        raise SystemExit("technical recovery gate failed")
    for name, expected in gate["frozen_predict_score_sha256"].items():
        if sha(HERE / name) != expected:
            raise SystemExit(f"frozen PREDICT score drift: {name}")
    predict_donors = load_predict_donors()
    if len(predict_donors) != 27:
        raise SystemExit("frozen PREDICT donor package incomplete")
    cells, samples, donors, groups, scores = elmentaite_unsorted_csr()
    elmentaite_donors, ecounts = frozen.write_scores("ELMENTAITE", cells, samples, donors, groups, scores)
    print(f"Elmentaite supportive scores: {ecounts}", flush=True)
    rows, numeric = frozen.assemble(predict_donors, elmentaite_donors)
    counts = {
        "FGID": {"cells": 31294, "samples": 13, "donors": 13, "reused_frozen_score_package": True},
        "CD": {"cells": 12281, "samples": 14, "donors": 14, "reused_frozen_score_package": True},
        "ELMENTAITE": ecounts,
    }
    result = {"terminal_token": "RA08_B1R_ANALYSIS_COMPLETE", "registered_members": len(rows), "numerically_formal_members": numeric, "terminal_SP08_NA_members": 3, "PREDICT": {"CD_donors": sum(g == "CD" for _, g in predict_donors), "FGID_donors": sum(g == "FGID" for _, g in predict_donors)}, "ELMENTAITE": {"CD_donors": sum(g == "CD" for _, g in elmentaite_donors), "control_donors": sum(g == "HEALTHY_NON_IBD" for _, g in elmentaite_donors), "formal_status": "NOT_ESTIMABLE"}, "score_counts": counts, "BH_family": "one numeric-P family after all 33 terminal states", "technical_recovery": "unsorted unique Elmentaite CSR indices; no PREDICT rescore", "NEXT_PHASE_AUTOSTART": "FORBIDDEN"}
    (HERE / "RA08_B1R_R1_ANALYSIS_STATUS.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()

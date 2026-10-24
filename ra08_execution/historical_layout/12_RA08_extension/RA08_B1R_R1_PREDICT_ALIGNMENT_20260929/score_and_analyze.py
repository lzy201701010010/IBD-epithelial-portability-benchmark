"""Execute the frozen RA08 B1R SP/SM family after the R1 input gate.

The PREDICT loader inverts the verified natural-log1p integer encoding exactly.
The score equations and donor inference follow the frozen RA07A3R1 methods and
RA06C5 numerical implementation; this file changes only source loading.
"""

import csv
import hashlib
import json
import math
import sys
from array import array
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
RAW_PREDICT = ROOT / "12_RA08_extension/RA08_B1R_PROSPECTIVE_EXTERNAL_EVIDENCE_20260928/raw/predict"
B1A = ROOT / "12_RA08_extension/RA08_B1A_CONTINUATION_20260928"
sys.path.insert(0, str(OUT / "_deps"))
import h5py  # noqa: E402
import numpy as np  # noqa: E402
from scipy import optimize, stats  # noqa: E402


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def programs():
    path = ROOT / "03_data_registry/RA01_STATE_PROGRAM_CANDIDATE_REGISTER.csv"
    with path.open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    return {r["program_id"]: r["genes_or_definition"].split(";") for r in rows if r["program_id"] != "SP08"}


PROGRAMS = programs()
TARGET_GENES = sorted({gene for genes in PROGRAMS.values() for gene in genes})
TARGET_INDEX = {gene: i for i, gene in enumerate(TARGET_GENES)}


def score_arrays(target, libraries, rank_target, target_genes=TARGET_GENES):
    output = {}
    gene_to_col = {gene: i for i, gene in enumerate(target_genes)}
    scale = np.divide(10000.0, libraries, out=np.zeros_like(libraries, dtype=np.float64), where=libraries > 0)
    for program, genes in PROGRAMS.items():
        columns = [gene_to_col[g] for g in genes]
        counts = target[:, columns].astype(np.float64)
        output[(program, "SM01")] = np.log1p(counts * scale[:, None]).mean(axis=1)
        output[(program, "SM02")] = rank_target[:, columns].mean(axis=1)
        output[(program, "SM03")] = (counts > 0).mean(axis=1)
    return output


def ranks_for_cell(positive_values, target_values, n_features):
    positive = np.sort(np.asarray(positive_values))
    zero_count = n_features - len(positive)
    if zero_count < 0:
        raise ValueError("more stored values than features")
    targets = np.asarray(target_values)
    result = np.empty(len(targets), dtype=np.float64)
    zero = targets == 0
    result[zero] = ((1.0 + zero_count) / 2.0 - 1.0) / (n_features - 1.0)
    if np.any(~zero):
        active = targets[~zero]
        left = np.searchsorted(positive, active, side="left")
        right = np.searchsorted(positive, active, side="right")
        rank = zero_count + left + (right - left + 1.0) / 2.0
        result[~zero] = (rank - 1.0) / (n_features - 1.0)
    return result


def write_scores(prefix, cell_ids, samples, donors, groups, scores):
    output = OUT / f"RA08_B1R_R1_{prefix}_EPITHELIAL_CELL_SCORES.npz"
    if output.exists():
        raise FileExistsError(output)
    payload = {"cell_id": np.asarray(cell_ids), "sample": np.asarray(samples), "donor": np.asarray(donors), "group": np.asarray(groups)}
    for (program, method), values in scores.items():
        if len(values) != len(cell_ids) or not np.isfinite(values).all():
            raise ValueError(f"nonfinite or length mismatch {prefix} {program} {method}")
        payload[f"{program}_{method}"] = np.asarray(values, dtype=np.float64)
    np.savez_compressed(output, **payload)
    sample_indices = defaultdict(list)
    donor_samples = defaultdict(set)
    for i, (sample, donor, group) in enumerate(zip(samples, donors, groups)):
        sample_indices[(sample, donor, group)].append(i)
        donor_samples[(donor, group)].add(sample)
    sample_scores = {}
    for key, indices in sample_indices.items():
        sample_scores[key] = {member: float(np.mean(values[indices])) for member, values in scores.items()}
    donor_scores = {}
    for donor, group in donor_samples:
        keys = sorted([k for k in sample_scores if k[1] == donor and k[2] == group])
        donor_scores[(donor, group)] = {member: float(np.mean([sample_scores[k][member] for k in keys])) for member in scores}
    sample_path = OUT / f"RA08_B1R_R1_{prefix}_SAMPLE_SCORES.tsv"
    donor_path = OUT / f"RA08_B1R_R1_{prefix}_DONOR_SCORES.tsv"
    fields = ["sample", "donor", "group", "epithelial_cells"] + [f"{p}_{m}" for p, m in scores]
    with sample_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for key in sorted(sample_scores):
            writer.writerow([*key, len(sample_indices[key])] + [sample_scores[key][member] for member in scores])
    with donor_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["donor", "group", "sample_count"] + [f"{p}_{m}" for p, m in scores])
        for key in sorted(donor_scores):
            writer.writerow([*key, len(donor_samples[key])] + [donor_scores[key][member] for member in scores])
    return donor_scores, {"cells": len(cell_ids), "samples": len(sample_scores), "donors": len(donor_scores), "cell_score_file": output.name, "sample_score_file": sample_path.name, "donor_score_file": donor_path.name}


def predict_arm(group, matrix_filename, barcode_filename, gene_filename, metadata_filename):
    with (RAW_PREDICT / barcode_filename).open(encoding="utf-8-sig") as f:
        barcodes = [line.strip() for line in f]
    with (RAW_PREDICT / gene_filename).open(encoding="utf-8-sig", newline="") as f:
        symbols = [row[1] for row in csv.reader(f, delimiter="\t")]
    gene_to_row = {gene: symbols.index(gene) + 1 for gene in TARGET_GENES}
    row_to_target = {row: TARGET_INDEX[gene] for gene, row in gene_to_row.items()}
    with (RAW_PREDICT / metadata_filename).open(encoding="utf-8-sig", newline="") as f:
        metadata = {row["NAME"]: row for row in csv.DictReader(f) if row["NAME"] != "TYPE"}
    study = "SCP1422" if group == "FGID" else "SCP1423"
    with (ROOT / "12_RA08_extension/RA08_B0_8R_FINAL_BRIDGE_ADJUDICATION_20260924/RA08_B0_8R_PREDICT_EPITHELIAL_LABEL_FREEZE_DRAFT.csv").open(encoding="utf-8-sig", newline="") as f:
        labels = {(r["study_id"], r["exact_author_label"]): r for r in csv.DictReader(f)}
    selected_columns = []
    cell_ids, samples, donors, groups, libraries = [], [], [], [], []
    for i, barcode in enumerate(barcodes):
        record = metadata[barcode]
        policy = labels[(study, record["cell_type__custom"])]
        if policy["author_ontology_label"] != record["cell_type__ontology_label"]:
            raise ValueError("annotation ontology mismatch")
        if policy["draft_status"] == "INCLUDE_EPITHELIAL":
            selected_columns.append(i)
            cell_ids.append(barcode)
            samples.append(record["biosample_id"])
            donors.append(record["donor_id"])
            groups.append(group)
            libraries.append(int(float(record["number_of_reads"])))
    selected_map = {col: i for i, col in enumerate(selected_columns)}
    target = np.zeros((len(cell_ids), len(TARGET_GENES)), dtype=np.uint16)
    positive_by_cell = [array("H") for _ in cell_ids]
    value_cache = {}
    with (RAW_PREDICT / matrix_filename).open("rb", buffering=1024 * 1024) as f:
        f.readline()
        f.readline()
        for line in f:
            row_text, col_text, value_text = line.split()
            cell = selected_map.get(int(col_text) - 1)
            if cell is None:
                continue
            count = value_cache.get(value_text)
            if count is None:
                inverse = math.expm1(float(value_text))
                count = round(inverse)
                if count < 1 or count > 65535 or abs(inverse - count) > 1e-8 * count:
                    raise ValueError("PREDICT source value invalid for frozen count loader")
                value_cache[value_text] = count
            positive_by_cell[cell].append(count)
            target_col = row_to_target.get(int(row_text))
            if target_col is not None:
                target[cell, target_col] = count
    rank_target = np.empty_like(target, dtype=np.float64)
    for cell, positive in enumerate(positive_by_cell):
        rank_target[cell] = ranks_for_cell(positive, target[cell], len(symbols))
    scores = score_arrays(target, np.asarray(libraries, dtype=np.float64), rank_target)
    del positive_by_cell
    return cell_ids, samples, donors, groups, scores


def elmentaite():
    bridge_path = B1A / "ELMENTAITE_CELL_SAMPLE_PARTICIPANT_BRIDGE.tsv"
    with bridge_path.open(encoding="utf-8-sig", newline="") as f:
        bridge = {r["source_cell_id"]: r for r in csv.DictReader(f, delimiter="\t")}
    h5ad = B1A / "raw/final_pediatric_object_cellxgene_c.h5ad"
    cell_ids, samples, donors, groups, libraries = [], [], [], [], []
    target_rows = []
    positive_rows = []
    with h5py.File(h5ad, "r") as h:
        obs_ids = [v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in h["obs/index"][:]]
        symbols = [v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in h["raw/var/index"][:]]
        target_columns = [symbols.index(gene) for gene in TARGET_GENES]
        if len(set(symbols)) != len(symbols) or len(obs_ids) != len(bridge):
            raise ValueError("Elmentaite frozen feature or cell identity drift")
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
            if np.any(np.diff(indices) <= 0):
                raise ValueError("Elmentaite CSR row indices not sorted and unique")
            lookup = dict(zip(indices.tolist(), data.tolist()))
            target_rows.append([lookup.get(column, 0) for column in target_columns])
            positive_rows.append(np.asarray(data, dtype=np.float64))
            cell_ids.append(cell)
            samples.append(record["source_sample_name"])
            donors.append(record["participant_id"])
            groups.append("CD" if record["ra08_group"] == "CROHN_DISEASE" else "HEALTHY_NON_IBD")
            libraries.append(float(np.sum(data, dtype=np.float64)))
    target = np.asarray(target_rows, dtype=np.float64)
    rank_target = np.empty_like(target, dtype=np.float64)
    for i, positive in enumerate(positive_rows):
        rank_target[i] = ranks_for_cell(positive, target[i], len(symbols))
    scores = score_arrays(target, np.asarray(libraries), rank_target)
    return cell_ids, samples, donors, groups, scores


def hedges_g(a, b):
    pooled = math.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    correction = 1.0 - 3.0 / (4.0 * (len(a) + len(b)) - 9.0)
    return float(correction * (a.mean() - b.mean()) / pooled)


def hedges_g_ci(a, b, alpha=0.05):
    n1, n0 = len(a), len(b)
    pooled = math.sqrt(((n1 - 1) * a.var(ddof=1) + (n0 - 1) * b.var(ddof=1)) / (n1 + n0 - 2))
    observed = (a.mean() - b.mean()) / (pooled * math.sqrt(1 / n1 + 1 / n0))
    df = n1 + n0 - 2
    correction = 1.0 - 3.0 / (4.0 * (n1 + n0) - 9.0)

    def solve_tail(kind, target):
        function = stats.nct.cdf if kind == "cdf" else stats.nct.sf
        for width in (4.0, 6.0, 8.0, 12.0, 20.0):
            lower, upper = observed - width, observed + width
            f_lower = function(observed, df, lower) - target
            f_upper = function(observed, df, upper) - target
            if np.isfinite(f_lower) and np.isfinite(f_upper) and f_lower * f_upper < 0:
                return float(optimize.brentq(lambda ncp: function(observed, df, ncp) - target, lower, upper))
        raise RuntimeError("NONCENTRAL_T_FINITE_BRACKET_NOT_FOUND")

    lower_ncp = solve_tail("cdf", 1 - alpha / 2)
    upper_ncp = solve_tail("sf", 1 - alpha / 2)
    scale = math.sqrt(1 / n1 + 1 / n0)
    return float(correction * lower_ncp * scale), float(correction * upper_ncp * scale)


def bh_numeric(rows):
    indices = [i for i, r in enumerate(rows) if r.get("P") is not None]
    m = len(indices)
    order = sorted(indices, key=lambda i: (rows[i]["P"], i))
    running = 1.0
    for position in range(m - 1, -1, -1):
        index = order[position]
        running = min(running, rows[index]["P"] * m / (position + 1))
        rows[index]["Q"] = min(1.0, running)
        rows[index]["BH_numeric_denominator"] = m
    return m


def assemble(predict_donors, elmentaite_donors):
    rows = []
    for program in [f"SP{i:02d}" for i in range(1, 12)]:
        for method in ("SM01", "SM02", "SM03"):
            row = {"program": program, "method": method, "predict_role": "FORMAL_SECONDARY", "elmentaite_role": "SUPPORTIVE_DESCRIPTIVE_ONLY"}
            if program == "SP08":
                row.update({"predict_status": "TERMINAL_NA_FROZEN_COMPOSITE", "elmentaite_status": "TERMINAL_NA_FROZEN_COMPOSITE", "reason": "SP08 frozen composite is not a gene-defined score for SM01-SM03"})
                rows.append(row)
                continue
            cd = np.array([s[(program, method)] for (donor, group), s in predict_donors.items() if group == "CD"], dtype=float)
            control = np.array([s[(program, method)] for (donor, group), s in predict_donors.items() if group == "FGID"], dtype=float)
            ecd = np.array([s[(program, method)] for (donor, group), s in elmentaite_donors.items() if group == "CD"], dtype=float)
            ect = np.array([s[(program, method)] for (donor, group), s in elmentaite_donors.items() if group == "HEALTHY_NON_IBD"], dtype=float)
            row.update({"predict_CD_donors": len(cd), "predict_FGID_donors": len(control), "predict_CD_mean": float(cd.mean()), "predict_FGID_mean": float(control.mean()), "predict_CD_minus_FGID": float(cd.mean() - control.mean()), "elmentaite_CD_donors": len(ecd), "elmentaite_healthy_nonIBD_donors": len(ect), "elmentaite_CD_mean": float(ecd.mean()), "elmentaite_control_mean": float(ect.mean()), "elmentaite_CD_minus_control": float(ecd.mean() - ect.mean()), "elmentaite_status": "SUPPORTIVE_DESCRIPTIVE_ONLY_NOT_FORMAL"})
            if len(cd) < 8 or len(control) < 8 or not np.isfinite(cd).all() or not np.isfinite(control).all() or cd.var(ddof=1) <= 0 or control.var(ddof=1) <= 0:
                row.update({"predict_status": "NOT_ESTIMABLE", "reason": "FROZEN_DONOR_OR_VARIANCE_GATE"})
            else:
                test = stats.ttest_ind(cd, control, equal_var=False, alternative="two-sided")
                va, vb = cd.var(ddof=1) / len(cd), control.var(ddof=1) / len(control)
                welch_df = (va + vb) ** 2 / (va ** 2 / (len(cd) - 1) + vb ** 2 / (len(control) - 1))
                try:
                    ci_low, ci_high = hedges_g_ci(cd, control)
                    row.update({"predict_status": "ESTIMABLE", "welch_statistic": float(test.statistic), "welch_df": float(welch_df), "P": float(test.pvalue), "hedges_g": hedges_g(cd, control), "hedges_g_ci_low": ci_low, "hedges_g_ci_high": ci_high})
                except RuntimeError:
                    row.update({"predict_status": "NOT_ESTIMABLE", "reason": "EXACT_NONCENTRAL_T_CI_FAILED"})
            rows.append(row)
    numeric = bh_numeric(rows)
    if len(rows) != 33 or sum(r["program"] == "SP08" for r in rows) != 3:
        raise ValueError("registered family not complete")
    fields = ["program", "method", "predict_role", "predict_status", "reason", "predict_CD_donors", "predict_FGID_donors", "predict_CD_mean", "predict_FGID_mean", "predict_CD_minus_FGID", "hedges_g", "hedges_g_ci_low", "hedges_g_ci_high", "welch_statistic", "welch_df", "P", "Q", "BH_numeric_denominator", "elmentaite_role", "elmentaite_status", "elmentaite_CD_donors", "elmentaite_healthy_nonIBD_donors", "elmentaite_CD_mean", "elmentaite_control_mean", "elmentaite_CD_minus_control"]
    with (OUT / "RA08_B1R_R1_EXTERNAL_EVIDENCE_MATRIX.tsv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return rows, numeric


def main():
    gate_path = OUT / "RA08_B1R_R1_PREDICT_FINAL_INPUT_GATE.json"
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    if gate.get("decision") != "PASS" or not all(gate["gates"][k] == "PASS" for k in gate["mandatory_gate_names"]):
        raise SystemExit("final input gate did not pass")
    if gate.get("scoring_source_sha256") != sha(Path(__file__)):
        raise SystemExit("scoring source hash mismatch")
    if (OUT / "RA08_B1R_R1_EXTERNAL_EVIDENCE_MATRIX.tsv").exists():
        raise SystemExit("scoring package already exists")
    predict_donors = {}
    counts = {}
    for group, files in {
        "FGID": ("gene_sorted-MMRNA_PREDICT_3p_Paper_FGIDFull.mtx", "PrefixedCustomBarcodes_FGID_3p_paper.tsv", "CustomGenes_FGID_3p_paper.tsv", "AlexandriaUpload_AlexandriaMetadataConvention_PREDICT_3p_Paper_FGID_final_cellOntologyLabels2.csv"),
        "CD": ("gene_sorted-MMRNA_PREDICT_3p_Paper_CD.mtx", "BarcodesCD_3p_paper.tsv", "GenesCD_3p_paper.tsv", "AlexandriaUpload_AlexandriaMetadataConvention_PREDICT_3p_Paper_CD_final_cellOntologyLabels2.csv"),
    }.items():
        cells, samples, donors, groups, scores = predict_arm(group, *files)
        d, counts[group] = write_scores(f"PREDICT_{group}", cells, samples, donors, groups, scores)
        predict_donors.update(d)
        print(f"PREDICT {group} scored: {counts[group]}", flush=True)
    cells, samples, donors, groups, scores = elmentaite()
    elmentaite_donors, counts["ELMENTAITE"] = write_scores("ELMENTAITE", cells, samples, donors, groups, scores)
    print(f"Elmentaite scored: {counts['ELMENTAITE']}", flush=True)
    rows, numeric = assemble(predict_donors, elmentaite_donors)
    result = {"terminal_token": "RA08_B1R_ANALYSIS_COMPLETE", "registered_members": len(rows), "numerically_formal_members": numeric, "terminal_SP08_NA_members": 3, "PREDICT": {"CD_donors": sum(g == "CD" for _, g in predict_donors), "FGID_donors": sum(g == "FGID" for _, g in predict_donors)}, "ELMENTAITE": {"CD_donors": sum(g == "CD" for _, g in elmentaite_donors), "control_donors": sum(g == "HEALTHY_NON_IBD" for _, g in elmentaite_donors), "formal_status": "NOT_ESTIMABLE"}, "score_counts": counts, "BH_family": "one numeric-P family after all 33 terminal states", "NEXT_PHASE_AUTOSTART": "FORBIDDEN"}
    (OUT / "RA08_B1R_R1_ANALYSIS_STATUS.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()

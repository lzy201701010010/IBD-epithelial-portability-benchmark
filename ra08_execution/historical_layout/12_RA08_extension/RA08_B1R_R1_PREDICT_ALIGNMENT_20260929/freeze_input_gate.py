import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
RA08 = ROOT / "12_RA08_extension"
RAW = RA08 / "RA08_B1R_PROSPECTIVE_EXTERNAL_EVIDENCE_20260928/raw/predict"
AMENDMENT = RA08 / "RA08_B1R_PROSPECTIVE_EXTERNAL_EVIDENCE_20260928/RA08_B1R_PROSPECTIVE_AMENDMENT.json"
LOCK = RA08 / "RA08_B1R_PROSPECTIVE_EXTERNAL_EVIDENCE_20260928/RA08_B1R_PREANALYSIS_LOCK.json"
DRAFT = RA08 / "RA08_B0_8R_FINAL_BRIDGE_ADJUDICATION_20260924/RA08_B0_8R_PREDICT_EPITHELIAL_LABEL_FREEZE_DRAFT.csv"


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def load(name):
    return json.loads((OUT / name).read_text(encoding="utf-8-sig"))


source = load("RA08_B1R_R1_SOURCE_HASH_AUDIT.json")
alignment = load("RA08_B1R_R1_ALIGNMENT_AUDIT.json")
cohort = load("RA08_B1R_R1_PREDICT_COHORT_AUDIT.json")
features = load("RA08_B1R_R1_FEATURE_AUDIT.json")
dryrun = load("RA08_B1R_R1_METHOD_DRYRUN.json")
scans = [load(f"RA08_B1R_R1_{g}_MATRIX_SCAN.json") for g in ("FGID", "CD")]
semantics = [load(f"RA08_B1R_R1_{g}_SEMANTICS_AUDIT.json") for g in ("FGID", "CD")]
amendment = json.loads(AMENDMENT.read_text(encoding="utf-8"))
lock = json.loads(LOCK.read_text(encoding="utf-8"))

authority = []
for item in amendment["source_authority_hashes"]:
    path = ROOT / item["path"]
    actual = sha(path)
    authority.append({"path": item["path"], "sha256": actual, "pass": actual == item["sha256"]})
for name, expected in [("RA08_B1R_PROSPECTIVE_AMENDMENT.md", lock["amendment_md_sha256"]), ("RA08_B1R_PROSPECTIVE_AMENDMENT.json", lock["amendment_json_sha256"])]:
    path = AMENDMENT.parent / name
    authority.append({"path": str(path.relative_to(ROOT)).replace("\\", "/"), "sha256": sha(path), "pass": sha(path) == expected})
acq_manifest = RA08 / "RA08_B1R_ACQ1_PREDICT_INPUT_RECOVERY_20260928/RA08_B1R_PREDICT_ACQ1_FROZEN_SOURCE_MANIFEST.tsv"
acq_supplement = json.loads((acq_manifest.parent / "RA08_B1R_PREDICT_ACQ1_INPUT_GATE_SUPPLEMENT.json").read_text(encoding="utf-8"))
authority.append({"path": str(acq_manifest.relative_to(ROOT)).replace("\\", "/"), "sha256": sha(acq_manifest), "pass": sha(acq_manifest) == acq_supplement["frozen_source_manifest_sha256"]})

fields = {
    "matrix_cell_id": "barcode TSV row, exactly equals metadata NAME and CellID",
    "sample": "biosample_id",
    "participant_donor": "donor_id",
    "source_phenotype": "disease",
    "RA08_group": {"SCP1422": "FGID / non-IBD / non-inflammatory comparator", "SCP1423": "Crohn disease"},
    "treatment_status": "disease__treated exactly FALSE",
    "baseline_timepoint": "biosample_id contains exact T0D component",
    "author_annotation": ["cell_type__custom", "cell_type__ontology_label"],
    "tissue": "organ exactly UBERON_0000331; source paper and portal describe ileal mucosa",
    "source": {"SCP1422": "FGID", "SCP1423": "CD"},
    "structural_noncell_record": "CSV line 2 NAME=TYPE, CellID=group, donor_id=group, biosample_id=group in each study; excluded from expression join as a schema declaration",
    "identifier_normalization": "NONE; exact byte-decoded equality",
}
(OUT / "RA08_B1R_R1_SOURCE_TO_RA08_FIELD_MAPPING.json").write_text(json.dumps(fields, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

counts = Counter()
for study, filename in [
    ("SCP1422", "AlexandriaUpload_AlexandriaMetadataConvention_PREDICT_3p_Paper_FGID_final_cellOntologyLabels2.csv"),
    ("SCP1423", "AlexandriaUpload_AlexandriaMetadataConvention_PREDICT_3p_Paper_CD_final_cellOntologyLabels2.csv"),
]:
    with (RAW / filename).open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            if row["NAME"] != "TYPE":
                counts[(study, row["cell_type__custom"], row["cell_type__ontology_label"])] += 1
with DRAFT.open(encoding="utf-8-sig", newline="") as f:
    draft_rows = list(csv.DictReader(f))
with (OUT / "RA08_B1R_R1_PREDICT_AUTHOR_LABEL_DISPOSITION_FROZEN.tsv").open("w", encoding="utf-8", newline="") as f:
    writer = csv.writer(f, delimiter="\t", lineterminator="\n")
    writer.writerow(["study_id", "exact_author_label", "author_ontology_label", "disposition", "observed_cells", "source_draft_sha256", "freeze_basis"])
    for r in draft_rows:
        writer.writerow([r["study_id"], r["exact_author_label"], r["author_ontology_label"], r["draft_status"], counts[(r["study_id"], r["exact_author_label"], r["author_ontology_label"])], sha(DRAFT), "unchanged author-linked ontology disposition; full-object exact-match pre-score freeze"])

gates = {
    "SOURCE_HASHES": "PASS" if source["all_pass"] and len(source["files"]) == 8 and all(x["pass"] for x in authority) else "FAIL",
    "MATRIX_STRUCTURE": "PASS" if all(x["structure_pass"] for x in scans) and [x["feature_rows"] for x in scans] == [22431, 22431] and [x["barcode_columns"] for x in scans] == [89849, 107432] else "FAIL",
    "BARCODE_ALIGNMENT": "PASS" if all(x["alignment_pass"] for x in alignment) else "FAIL",
    "METADATA_EXTRA_RECORD_RESOLVED": "PASS" if all(x["metadata_type_declaration_count"] == 1 and x["type_declaration"][0]["csv_line"] == 2 and x["type_declaration"][0]["NAME"] == "TYPE" and x["type_declaration"][0]["CellID"] == "group" for x in alignment) else "FAIL",
    "FEATURE_ALIGNMENT": "PASS" if features["feature_gate_pass"] and features["byte_identical"] else "FAIL",
    "MATRIX_SEMANTICS": "PASS" if all(x["raw_count_interpretation_pass"] for x in semantics) else "FAIL",
    "BASELINE_SELECTION": "PASS" if sum(x["cells"] for x in cohort) == 197281 and all(not x["non_T0D_samples"] and x["disease_treated_counts"] == {"FALSE": x["cells"]} for x in cohort) else "FAIL",
    "DONOR_BRIDGE": "PASS" if [x["donor_count"] for x in cohort] == [13, 14] and all(not x["samples_with_multiple_donors"] and not x["donors_with_multiple_samples"] and x["missing_sample_or_donor"] == 0 and x["identity_sample_prefix_failures"] == 0 for x in cohort) else "FAIL",
    "AUTHOR_ANNOTATION": "PASS" if all(not x["unmatched_label_counts"] and not x["ontology_mismatch_counts"] for x in cohort) and len(draft_rows) == 423 else "FAIL",
    "EPITHELIAL_FREEZE": "PASS" if sum(x["category_counts"].get("INCLUDE_EPITHELIAL", 0) for x in cohort) == 43575 and all(x["category_counts"].get("UNRESOLVED", 0) == 0 for x in cohort) else "FAIL",
    "PROGRAM_FEATURE_COVERAGE": "PASS" if all(p.get("coverage_pass", True) for p in features["programs"]) and sum(p.get("coverage_pass", False) for p in features["programs"]) == 10 else "FAIL",
    "METHOD_TECHNICAL_ESTIMABILITY": "PASS" if dryrun["status"] == "PASS" and all(x["raw_count_interpretation_pass"] for x in semantics) else "FAIL",
    "DONOR_ESTIMABILITY": "PASS" if [x["eligible_donors_by_sample_cell_gate"] for x in cohort] == [13, 14] else "FAIL",
}

result = {
    "gate_id": "RA08_B1R_R1_PREDICT_FINAL_INPUT_GATE",
    "frozen_utc": datetime.now(timezone.utc).isoformat(),
    "mandatory_gate_names": list(gates),
    "gates": gates,
    "decision": "PASS" if all(v == "PASS" for v in gates.values()) else "FAIL",
    "FGID_matrix_cells": 89849,
    "CD_matrix_cells": 107432,
    "baseline_cells": 197281,
    "epithelial_cells": 43575,
    "eligible_FGID_donors": cohort[0]["eligible_donors_by_sample_cell_gate"],
    "eligible_CD_donors": cohort[1]["eligible_donors_by_sample_cell_gate"],
    "matrix_semantics": "natural log1p of integer raw per-gene counts; inverse expm1 rounds to exact integers whose per-cell sum equals source number_of_reads and nonzero count equals nFeature_RNA in all 197281 cells",
    "loader_adaptation": "exact inverse of source log1p encoding before unchanged frozen SM01-SM03 formulas",
    "annotation_disposition_sha256": sha(OUT / "RA08_B1R_R1_PREDICT_AUTHOR_LABEL_DISPOSITION_FROZEN.tsv"),
    "source_to_RA08_mapping_sha256": sha(OUT / "RA08_B1R_R1_SOURCE_TO_RA08_FIELD_MAPPING.json"),
    "scoring_source_sha256": sha(OUT / "score_and_analyze.py"),
    "authority_checks": authority,
    "source_files_sha256": {x["official_filename"]: x["sha256"] for x in source["files"]},
    "preanalysis_lock_utc": lock["timestamp_utc"],
    "postscore_variance_and_finite_rules": "apply per registered member before Welch; failures receive terminal nonnumeric row, no pseudo-P or BH rank",
    "SP08": "three frozen terminal N/A rows, no numeric score",
    "NEXT_PHASE_AUTOSTART": "FORBIDDEN",
}
path = OUT / "RA08_B1R_R1_PREDICT_FINAL_INPUT_GATE.json"
if path.exists():
    raise FileExistsError(path)
path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(json.dumps({"decision": result["decision"], "gates": gates, "authority_pass": sum(x["pass"] for x in authority), "authority_total": len(authority), "scoring_source_sha256": result["scoring_source_sha256"]}, indent=2))

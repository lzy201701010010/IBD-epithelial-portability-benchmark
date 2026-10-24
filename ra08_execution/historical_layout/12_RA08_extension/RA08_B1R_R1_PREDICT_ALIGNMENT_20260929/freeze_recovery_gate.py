import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

OUT = Path(__file__).resolve().parent


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


input_gate = json.loads((OUT / "RA08_B1R_R1_PREDICT_FINAL_INPUT_GATE.json").read_text(encoding="utf-8"))
csr = json.loads((OUT / "RA08_B1R_R1_ELMENTAITE_CSR_ORDER_AUDIT.json").read_text(encoding="utf-8"))
files = [OUT / f"RA08_B1R_R1_PREDICT_{group}_{suffix}" for group in ("FGID", "CD") for suffix in ("EPITHELIAL_CELL_SCORES.npz", "SAMPLE_SCORES.tsv", "DONOR_SCORES.tsv")]
if input_gate["decision"] != "PASS" or not csr["unsorted_only_adaptation_permitted"] or any(not path.is_file() for path in files):
    raise SystemExit("technical recovery preconditions failed")
result = {
    "gate_id": "RA08_B1R_R1_TECHNICAL_RECOVERY_GATE",
    "frozen_utc": datetime.now(timezone.utc).isoformat(),
    "decision": "PASS",
    "failure": "Original scorer rejected unsorted unique Elmentaite raw.X CSR indices after PREDICT score files were written and before any formal inference",
    "csr_order_audit_sha256": sha(OUT / "RA08_B1R_R1_ELMENTAITE_CSR_ORDER_AUDIT.json"),
    "PREDICT_final_input_gate_sha256": sha(OUT / "RA08_B1R_R1_PREDICT_FINAL_INPUT_GATE.json"),
    "original_scoring_source_sha256": sha(OUT / "score_and_analyze.py"),
    "recovery_source_sha256": sha(OUT / "resume_unsorted_csr.py"),
    "frozen_predict_score_sha256": {path.name: sha(path) for path in files},
    "scientific_rule_change": False,
    "predict_rescore_permitted": False,
    "elmentaite_loader_change": "allow unsorted unique CSR feature indices; exact index lookup and frozen value-rank function unchanged",
    "post_result_rescue": False,
    "NEXT_PHASE_AUTOSTART": "FORBIDDEN",
}
path = OUT / "RA08_B1R_R1_TECHNICAL_RECOVERY_GATE.json"
if path.exists():
    raise FileExistsError(path)
path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps({"decision": result["decision"], "PREDICT_scores_frozen": len(files), "recovery_source_sha256": result["recovery_source_sha256"]}, indent=2))

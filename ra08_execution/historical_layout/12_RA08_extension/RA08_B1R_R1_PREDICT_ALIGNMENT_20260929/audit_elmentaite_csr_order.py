import json
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
sys.path.insert(0, str(OUT / "_deps"))
import h5py
import numpy as np

path = ROOT / "12_RA08_extension/RA08_B1A_CONTINUATION_20260928/raw/final_pediatric_object_cellxgene_c.h5ad"
with h5py.File(path, "r") as f:
    x = f["raw/X"]
    indptr = x["indptr"][:]
    indices = x["indices"][:]
    n_features = int(x.attrs["shape"][1])
unsorted_rows = 0
duplicate_rows = 0
out_of_bounds_rows = 0
first_unsorted = None
first_duplicate = None
for i in range(len(indptr) - 1):
    row = indices[indptr[i]:indptr[i + 1]]
    if len(row) and (int(row.min()) < 0 or int(row.max()) >= n_features):
        out_of_bounds_rows += 1
    if len(row) > 1 and np.any(np.diff(row) < 0):
        unsorted_rows += 1
        if first_unsorted is None:
            first_unsorted = i
    if len(row) > 1 and len(np.unique(row)) != len(row):
        duplicate_rows += 1
        if first_duplicate is None:
            first_duplicate = i
result = {"source": str(path.relative_to(ROOT)).replace("\\", "/"), "rows": len(indptr) - 1, "stored_coordinates": len(indices), "unsorted_rows": unsorted_rows, "duplicate_coordinate_rows": duplicate_rows, "out_of_bounds_rows": out_of_bounds_rows, "first_unsorted_row_zero_based": first_unsorted, "first_duplicate_row_zero_based": first_duplicate, "unsorted_only_adaptation_permitted": unsorted_rows > 0 and duplicate_rows == 0 and out_of_bounds_rows == 0}
(OUT / "RA08_B1R_R1_ELMENTAITE_CSR_ORDER_AUDIT.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "manifests" / "SHA256SUMS.csv"


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest().upper()


def main() -> None:
    failures = []
    with MANIFEST.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        path = ROOT / Path(row["path"])
        if not path.is_file() or path.stat().st_size != int(row["size"]) or digest(path) != row["sha256"]:
            failures.append(row["path"])
    expected_tables = {"Table_1.csv"} | {f"Supplementary_Table_S{i}.csv" for i in range(1, 8)}
    actual_tables = {path.name for path in (ROOT / "tables").glob("*.csv")}
    if actual_tables != expected_tables:
        failures.append("publication-table inventory")
    if failures:
        raise SystemExit("RELEASE_VERIFICATION=FAIL: " + ", ".join(failures))
    print(f"RELEASE_VERIFICATION=PASS files={len(rows)} tables={len(actual_tables)}")


if __name__ == "__main__":
    main()

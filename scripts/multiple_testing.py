from __future__ import annotations

from statsmodels.stats.multitest import multipletests


def adjust_families(rows: list[dict]) -> list[dict]:
    output = [dict(row) for row in rows]
    families = sorted({row["family_id"] for row in output})
    for family in families:
        positions = [index for index, row in enumerate(output) if row["family_id"] == family]
        adjusted = multipletests([float(output[index]["raw_p"]) for index in positions], method="fdr_bh")[1]
        for index, value in zip(positions, adjusted):
            output[index]["adjusted_p"] = float(value)
    return output


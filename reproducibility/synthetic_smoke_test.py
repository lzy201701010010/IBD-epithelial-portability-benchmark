from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from classify_portability import classify_method, classify_program
from score_detection_fraction import score as detection_score
from score_mean_logcp10k import counts_to_logcp10k, score as mean_score
from score_within_cell_rank import score as rank_score


def main() -> None:
    counts = np.asarray([[0, 1, 3], [2, 2, 2]], dtype=float)
    assert np.allclose(detection_score(counts), [2 / 3, 1])
    transformed = counts_to_logcp10k(counts)
    assert np.isfinite(mean_score(transformed)).all()
    assert np.allclose(rank_score(counts), [0.5, 0.5])
    assert classify_method(1, 1, True, True) == "STRONG_PORTABLE"
    assert classify_method(1, -1, True, True) == "DISCORDANT"
    assert classify_program(["STRONG_PORTABLE", "STRONG_PORTABLE", "CONSISTENT_INSUFFICIENT_EVIDENCE"]) == "PORTABLE"
    assert classify_program(["STRONG_PORTABLE", "DISCORDANT", "STRONG_PORTABLE"]) == "DISCORDANT"
    print("SYNTHETIC_SMOKE_TEST=PASS")


if __name__ == "__main__":
    main()

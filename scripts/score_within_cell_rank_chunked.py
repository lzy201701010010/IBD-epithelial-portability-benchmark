from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.stats import rankdata


def _normalized_columns(program_columns: Iterable[int] | None) -> list[int] | None:
    return None if program_columns is None else list(program_columns)


def _validate_batch_size(batch_size: int) -> int:
    if isinstance(batch_size, bool) or not isinstance(batch_size, (int, np.integer)):
        raise TypeError("batch_size must be a positive integer")
    if int(batch_size) <= 0:
        raise ValueError("batch_size must be a positive integer")
    return int(batch_size)


def _score_dense_block(values: np.ndarray, program_columns: list[int] | None) -> np.ndarray:
    if values.shape[1] == 1:
        return np.ones(values.shape[0])
    ranks = np.apply_along_axis(rankdata, 1, values, method="average")
    percentile = (ranks - 1.0) / (values.shape[1] - 1.0)
    if program_columns is not None:
        percentile = percentile[:, program_columns]
    return np.mean(percentile, axis=1)


def _score_block(matrix, start: int, stop: int, program_columns: list[int] | None) -> np.ndarray:
    block = matrix[start:stop]
    values = block.toarray() if hasattr(block, "toarray") else np.asarray(block)
    return _score_dense_block(values, program_columns)


def score(matrix, program_columns=None, batch_size: int = 256) -> np.ndarray:
    """Numerically equivalent, row-chunked implementation of frozen SM02.

    Each row batch is densified across the complete input feature universe. Ranking,
    average-tie handling, percentile normalization, program-column selection and the
    final mean exactly follow the frozen ``score_within_cell_rank.py`` implementation.
    ``batch_size`` is a technical parameter only.
    """

    batch_size = _validate_batch_size(batch_size)
    source = matrix if hasattr(matrix, "shape") else np.asarray(matrix)
    if len(source.shape) != 2:
        raise ValueError("matrix must be two-dimensional")
    n_rows, n_features = source.shape
    columns = _normalized_columns(program_columns)
    if n_features == 1:
        return np.ones(n_rows)
    output = np.empty(n_rows, dtype=np.float64)
    for start in range(0, n_rows, batch_size):
        stop = min(start + batch_size, n_rows)
        output[start:stop] = _score_block(source, start, stop, columns)
    return output


def _atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = Path(str(path) + ".partial")
    partial.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    os.replace(partial, path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def score_to_npy_checkpointed(
    matrix,
    output_path: str | Path,
    checkpoint_path: str | Path,
    program_columns=None,
    batch_size: int = 256,
    stop_after_batches: int | None = None,
) -> Path:
    """Write SM02 scores with resumable batches and atomic final promotion.

    ``stop_after_batches`` is a synthetic-test hook. It raises ``InterruptedError``
    only after the current batch and checkpoint have both been flushed.
    """

    batch_size = _validate_batch_size(batch_size)
    source = matrix if hasattr(matrix, "shape") else np.asarray(matrix)
    if len(source.shape) != 2:
        raise ValueError("matrix must be two-dimensional")
    n_rows, n_features = map(int, source.shape)
    columns = _normalized_columns(program_columns)
    output_path = Path(output_path)
    checkpoint_path = Path(checkpoint_path)
    partial_output = Path(str(output_path) + ".partial")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.exists():
        raise FileExistsError(f"final output already exists: {output_path}")

    if checkpoint_path.exists():
        state = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        expected = {
            "n_rows": n_rows,
            "n_features": n_features,
            "program_columns": columns,
            "dtype": "float64",
        }
        for key, value in expected.items():
            if state.get(key) != value:
                raise ValueError(f"checkpoint mismatch for {key}")
        if state.get("completed") is True:
            raise ValueError("checkpoint says completed but final output is absent")
        if not partial_output.exists():
            raise FileNotFoundError("checkpoint exists but partial output is absent")
        next_row = int(state["next_row"])
        output = np.lib.format.open_memmap(partial_output, mode="r+")
        if output.shape != (n_rows,) or output.dtype != np.float64:
            raise ValueError("partial output shape or dtype mismatch")
    else:
        if partial_output.exists():
            raise FileExistsError("untracked partial output exists")
        next_row = 0
        output = np.lib.format.open_memmap(
            partial_output, mode="w+", dtype=np.float64, shape=(n_rows,)
        )
        _atomic_json(
            checkpoint_path,
            {
                "completed": False,
                "next_row": 0,
                "n_rows": n_rows,
                "n_features": n_features,
                "program_columns": columns,
                "dtype": "float64",
                "last_batch_size": batch_size,
            },
        )

    batches_completed = 0
    for start in range(next_row, n_rows, batch_size):
        stop = min(start + batch_size, n_rows)
        output[start:stop] = _score_block(source, start, stop, columns)
        output.flush()
        _atomic_json(
            checkpoint_path,
            {
                "completed": False,
                "next_row": stop,
                "n_rows": n_rows,
                "n_features": n_features,
                "program_columns": columns,
                "dtype": "float64",
                "last_batch_size": batch_size,
            },
        )
        batches_completed += 1
        if stop_after_batches is not None and batches_completed >= stop_after_batches:
            del output
            raise InterruptedError("synthetic checkpoint interruption requested")

    del output
    os.replace(partial_output, output_path)
    _atomic_json(
        checkpoint_path,
        {
            "completed": True,
            "next_row": n_rows,
            "n_rows": n_rows,
            "n_features": n_features,
            "program_columns": columns,
            "dtype": "float64",
            "last_batch_size": batch_size,
            "output_sha256": _sha256(output_path),
        },
    )
    return output_path

from __future__ import annotations

from collections import Counter


def map_program(program_genes: list[str], feature_symbols: list[str], minimum_fraction: float = 0.8, minimum_genes: int = 3) -> dict:
    normalized = [value.upper() for value in feature_symbols if value]
    counts = Counter(normalized)
    conflicts = sorted(gene for gene in program_genes if counts[gene.upper()] > 1)
    mapped = [gene for gene in program_genes if counts[gene.upper()] == 1]
    fraction = len(mapped) / len(program_genes) if program_genes else 0.0
    estimable = not conflicts and len(mapped) >= minimum_genes and fraction >= minimum_fraction
    return {
        "mapped_genes": mapped,
        "mapping_fraction": fraction,
        "duplicate_conflicts": conflicts,
        "estimable": estimable,
    }


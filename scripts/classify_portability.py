from __future__ import annotations


SUPPORTED = {"STRONG_PORTABLE"}
CONCORDANT = {"STRONG_PORTABLE", "DIRECTIONALLY_CONSISTENT_LIMITED", "CONSISTENT_INSUFFICIENT_EVIDENCE"}


def classify_method(direction_a: int, direction_b: int, supported_a: bool, supported_b: bool) -> str:
    """Apply the frozen direction-and-within-cohort-support rule to one method."""
    if direction_a == 0 or direction_b == 0:
        return "ZERO_OR_UNDEFINED"
    if direction_a != direction_b:
        return "DISCORDANT"
    if supported_a and supported_b:
        return "STRONG_PORTABLE"
    if supported_a or supported_b:
        return "DIRECTIONALLY_CONSISTENT_LIMITED"
    return "CONSISTENT_INSUFFICIENT_EVIDENCE"


def classify_program(method_classes: list[str], role: str = "PRIMARY_HELDOUT", estimable: bool = True) -> str:
    """Apply role/non-estimability precedence followed by frozen program-level rules."""
    if role in {"REFERENCE_ONLY", "NEGATIVE_COMPARATOR", "FAILURE_AUDIT", "GOVERNANCE_ONLY"}:
        return role
    if not estimable:
        return "NOT_ESTIMABLE"
    if "DISCORDANT" in method_classes:
        return "DISCORDANT"
    strong = sum(value in SUPPORTED for value in method_classes)
    concordant = sum(value in CONCORDANT for value in method_classes)
    if strong >= 2:
        return "PORTABLE"
    if concordant >= 2 and strong >= 1:
        return "PARTIALLY_PORTABLE_CONTEXT_DEPENDENT"
    if concordant >= 2:
        return "DIRECTIONALLY_CONSISTENT_INSUFFICIENT"
    return "INCONCLUSIVE"

# Historical workflow and expected artifacts

The following is documentation, not authorization to run. This task ran none of these scripts.

1. freeze_input_gate.py consumes historical source, alignment, feature, cohort, semantics and dry-run audit records plus the prospective amendment/lock. It hashes registered sources and writes the final input gate, field map and annotation disposition. All 13 mandatory gates must pass; the frozen scoring-source hash is required. Do not regenerate a passing record to avoid missing provenance.
2. score_and_analyze.py checks that gate and its source hash, reads the exact PREDICT sources, writes FGID/CD cell NPZ plus sample/donor TSVs, then attempts Elmentaite. The historical original rejected unsorted CSR indices before final inference. This is a preserved historical failure, not a runnable final single-command pipeline.
3. audit_elmentaite_csr_order.py inspects the authorized Elmentaite raw.X CSR ordering, duplicates and bounds. freeze_recovery_gate.py freezes the order audit, original gate, both code hashes and all six PREDICT checkpoint hashes. It does not authorize new data access.
4. resume_unsorted_csr.py checks the recovery gate and all six checkpoints, reads both PREDICT donor TSVs and the authorized Elmentaite container/bridge, and permits unsorted unique indices with unchanged equations. It requires exactly 27 PREDICT donors. It writes the Elmentaite scores, 33-member external evidence matrix and analysis status. The frozen family has 30 numerical members and three SP08 terminal entries.
5. independent_verify_results.py reads cell NPZ and donor TSV for all three arms, the evidence matrix and Elmentaite raw.X. It historically checks cell-to-donor means, differences, P and q values; it is not an additional biological validation. Its imports perform work. Do not import it for inspection.

Expected results: 14 CD/13 FGID PREDICT donors; 7 CD/8 non-IBD Elmentaite donors, formal NOT_ESTIMABLE; original portability and eligibility unchanged. No numerical outputs are regenerated or shipped in this code candidate. Existing public aggregate Table S8 remains a separate publication-support artifact.

All code hashes, count/feature identities, full rank universe, annotation policy, donor aggregation, no-fallback rules and checkpoint gates stay intact. Copying this folder alone does not supply the missing historical inputs, audits or authorization. No fresh synthetic or smoke execution was undertaken because it would import the historical scientific implementation.

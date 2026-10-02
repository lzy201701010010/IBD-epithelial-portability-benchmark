# RA08 reader-facing reproducibility update

Version 1.1.0 adds publication-support materials for the post-classification RA08 external-evidence layer. It does not change the historical Phase 3G classification, program definitions, held-out eligibility, or the scientific results frozen for the manuscript. RA08 is a separate later evidence layer; this release is not a prospective registration.

## Files

- `Supplementary_Table_S3_data.csv`: complete historical 90-member multiplicity-family mapping plus the zero-member closure record. The v1.0.0 `tables/Supplementary_Table_S3.csv` remains the family summary.
- `RA08_B4_Supplementary_Table_S8_complete_exact.tsv`: exact 33-row post-classification evidence table, including 30 numerical PREDICT comparisons and three terminal SP08 rows. Elmentaite is descriptive and non-estimable for formal inference.
- `RA08_B3_EXTERNAL_EVIDENCE_OVERVIEW_SOURCE.tsv`: the 66 display-level signs and statuses used for Figure S3; no cell- or donor-level source rows.
- `RA08_B4_Supplementary_Figure_S3_source.R`: rendering source for that display TSV. Run with the TSV in the same directory and R packages `ggplot2` and `svglite`; pass the `ra08/` directory as its argument.
- `RA08_B1R_R1_SOURCE_TO_RA08_FIELD_MAPPING.json`: reader-facing source-field and rule mapping, without source annotation rows.
- `SHA256SUMS.tsv`: SHA-256 values for the five materials above and this README.

These materials are author-created derived tables, display data, code and documentation. No source expression datasets, source barcodes, source sample/donor mapping rows, protected SCP1884 row-level material, or third-party raw data are redistributed. Obtain third-party data from their original custodians under their own access terms. SCP1884 remains subject to its repository access conditions. The RA08 execution code, input files and exact controlled-data reconstruction pathway are not included, so this release does not provide turnkey reproduction of every later result.

The repository's MIT License covers author-owned code and documentation; it does not relicense third-party data. The derived tables and display data are supplied as publication-support material. Source identifiers in these files refer to original records and do not grant rights over those records.

# Donor-level benchmarking of cross-cohort portability for epithelial programs in inflammatory bowel disease

This repository is the version 1.2.0 publication reproducibility release for a donor-level benchmark of predefined epithelial-program portability across inflammatory bowel disease transcriptomic cohorts.

## Objective and scope

The study evaluates whether frozen epithelial programs retain donor-level direction and within-cohort multiplicity-adjusted support across prespecified cohorts and scoring methods. This repository supports inspection and verification of publication-level outputs. It is not an entire project dump and does not claim a turnkey reconstruction of third-party source datasets.

## Source datasets

The benchmark used repository records GSE305528, S-BSST2944, E-MTAB-16998, E-MTAB-16999, GSE266546, GSE214695 and SCP1884. Controlled raw IBDverse materials under EGAD00001015692 were not used. Obtain source data from their original custodians under the applicable terms.

## Included

- Frozen program, scoring, estimability, multiplicity and classification specifications.
- Author-owned deterministic scoring and statistical helper modules.
- Publication-facing Table 1 and Supplementary Tables S1-S7.
- Safety-filtered figure display/source data; no figure files are modified or redistributed here.
- Environment records, synthetic checks and SHA-256 verification.

## Version history and RA08 update

Version 1.0.0 is the original benchmark release and remains available unchanged under its historical tag and Zenodo record. Version 1.1.0 adds `ra08/` with the complete historical S3 member mapping, exact post-classification Table S8, Figure S3 display data and R rendering source, a reader-facing RA08 field/rule map, and checksums. See [`ra08/README.md`](ra08/README.md) for file roles and access limits.

The update does not change the Phase 3G classification or held-out eligibility. RA08 is a separate post-classification evidence layer, not a new prospective registration. No third-party source dataset is redistributed; SCP1884 remains repository governed. Version 1.1.0 did not include the RA08 execution code. Version 1.2.0 adds the six frozen historical scripts and safe program/schema configuration in `ra08_execution/` under an explicit MIT licence. It includes run-order and dependency documentation, with all original code hashes preserved; no scientific analysis was rerun. The release does not supply source datasets, bridges or checkpoints and does not demonstrate end-to-end reproduction.

## Intentionally excluded

No raw or duplicated source datasets, restricted row-level metadata, participant or sample identifiers, access-control artifacts, private correspondence, internal project-control packages, caches, credentials or local paths are included.

## Verification

Create the documented environment, then run:

```text
python reproducibility/synthetic_smoke_test.py
python reproducibility/verify_release.py
```

The first command tests core deterministic rules on synthetic inputs. The second verifies the original release manifest and publication-table schema. For the version 1.1.0 update, compare `ra08/SHA256SUMS.tsv` with the six files in `ra08/`. Reproducing analyses from source data requires obtaining the relevant third-party datasets and respecting their access terms.

## Citation and license

Use `CITATION.cff` for creator and version metadata. Author-owned code and original repository documentation are released under the MIT License. Third-party datasets, source publications and repository records remain under their own terms and are not relicensed by this repository.

## Version 1.2.0 scope and current manuscript

See `ra08_execution/README.md`, `LICENSE_SCOPE.md` and `FILE_MANIFEST.csv`. Existing scientific code, specifications, aggregate tables and display data are preserved byte for byte. The historic table CSVs remain archival publication-support versions; the corrected manuscript Supplementary Table S1 takes precedence for its cohort-role wording. No manuscript or new scientific result is published in this code release. `manifests/SHA256SUMS.csv` describes this version; `manifests/historical_v1.1.0_SHA256SUMS.csv` describes the prior version only.

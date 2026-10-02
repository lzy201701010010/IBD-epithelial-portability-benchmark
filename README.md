# Donor-level benchmarking of cross-cohort portability for epithelial programs in inflammatory bowel disease

This repository is the version 1.1.0 publication reproducibility release for a donor-level benchmark of predefined epithelial-program portability across inflammatory bowel disease transcriptomic cohorts.

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

The update does not change the Phase 3G classification or held-out eligibility. RA08 is a separate post-classification evidence layer, not a new prospective registration. No third-party source dataset is redistributed; SCP1884 remains repository governed. The later RA08 execution code is still local, so this version does not enable complete turnkey reproduction from controlled source data.

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

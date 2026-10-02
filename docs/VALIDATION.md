# Validation

Run checks from the repository root. Synthetic fixtures are not biological results; real SCC acceptance additionally checks actual source access, exact reference/allele identity and native outputs. This repository does not contain participant-level data or test results copied from controlled cohorts.

## Runtime preflight

In SCC RStudio select R 4.5.2. The shared packages and native tool modules are used by both the launcher and the Run App button. A quick Terminal check is:

```bash
module load R/4.5.2 samtools/1.23 regtools/1.0.0 bcftools/1.23
Rscript scripts/check_access.R
```

Loading modules in the RStudio Terminal does not change the already-running R Console. `runtime_setup.R` performs the corresponding native-path initialization inside the R process. It does not install software or import the complete login environment.

## Synthetic checks

The native suites generate their own small BAM/VCF fixtures and exercise the same tools used by the application:

```bash
Rscript tests/run_tests.R
Rscript tests/test_contigs.R
Rscript tests/test_variants.R
Rscript tests/test_rna_bases.R
Rscript tests/test_all_samples.R --native
Rscript tests/test_variant_jobs.R
Rscript tests/test_annotations.R
Rscript tests/test_splice_evidence.R --native
Rscript tests/test_data_sources.R --native
```

The UI suites check immutable source switching, asynchronous stale-result suppression, failure retention, matched counts and plotting semantics:

```bash
Rscript tests/test_autoload.R app.R
Rscript tests/test_variant_ui.R
Rscript tests/test_annotation_ui.R
Rscript tests/test_splice_evidence_ui.R
Rscript tests/test_data_sources_ui.R
Rscript tests/test_single_bam_ui.R
```

Run any runtime/launcher-specific tests listed in `tests/` alongside these suites. Generated fixtures stay in temporary directories. Do not run uncontrolled full-cohort analyses on a login node.

## SCC acceptance

A release check should use a fresh clone/checkout in an independent writable path, clear personal `REGSHINY_*` deployment overrides, then prepare FHS from its original CSV and VCF/BAM sources. Validate a real coordinate lookup and bounded RNA preview, and also connect custom BAM/VCF data through explicit and read-group sample matching. Confirm that genotype groups, read-depth units, known-zero/missing states and data-source identity remain correct.

For scheduler runs, success requires a matching application receipt, artifact/source checks and terminal `qacct` with both `failed=0` and `exit_status=0`. Queue disappearance is not sufficient. Browser verification should use a real SCC RStudio session, check default FHS and a custom source, and save screenshots outside tracked source.

Private run-specific receipts, manifests, data and screenshots belong in the user's validation/state location. Do not commit them or label a current checkout verified based only on an unrelated historical `validation_status.json` or stale checksum file.

Repository access, SCC source-file permissions, software correctness and biological interpretation are separate checks. Other members use their own authenticated GitHub/SCC accounts; an owner's successful run does not claim an independent member login was tested.

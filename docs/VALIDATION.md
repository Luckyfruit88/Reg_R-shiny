# Validation

Run checks from the repository root. Synthetic fixtures are not biological results; real SCC acceptance additionally checks actual source access, exact reference/allele identity and native outputs. This repository does not contain participant-level data or test results copied from controlled cohorts.

## Runtime preflight

In SCC RStudio select R 4.5.2. The shared packages and native tool modules are used by both the launcher and the Run App button. A quick Terminal check is:

```bash
module load R/4.5.2 samtools/1.23 regtools/1.0.0 bcftools/1.23
Rscript scripts/check_access.R
```

Loading modules in the RStudio Terminal does not change the already-running R Console. `src/runtime_setup.R` performs the corresponding native-path initialization inside the R process. It does not install software or import the complete login environment.

## Synthetic checks

One command runs every synthetic R suite and the independent Python arithmetic oracle, each in its own process, and fails if any suite fails or if the run changes the git worktree (new, removed or rewritten status entries, compared by content):

```bash
Rscript tests/run_tests.R                  # full: native suites need samtools/1.23 regtools/1.0.0 bcftools/1.23
Rscript tests/run_tests.R --unit           # native-free subset
Rscript tests/run_tests.R --unit --require-clean   # as run in CI: worktree must be clean before and after
```

Use `--require-clean` from a clean checkout for release acceptance. Without it, local changes are allowed and only reported as unaltered. A git query that fails is reported as unverified and fails the run. Full per-suite logs are kept in the ignored `test-artifacts/test-runs/<run>/` directory; the terminal shows a summary and the tail of each failing log.

The native suites generate their own small BAM/VCF fixtures and exercise the same tools used by the application. The UI suites check immutable source switching, asynchronous stale-result suppression, failure retention, matched counts and plotting semantics. A single suite can still be run directly, for example `Rscript tests/test_variants.R --unit`; the suite list and per-suite flags are in `tests/run_tests.R`.

Generated fixtures stay in temporary directories. Do not run uncontrolled full-cohort analyses on a login node.

## Frontend and browser checks

`tests/test_frontend_contract.R` renders the production WGS UI, checks unique namespaced controls and executes the production display-control and display-range expressions. `tests/test_variant_ui.R` uses synthetic job mocks for cancellation, stale batch confirmation, replay, source deactivation, independent per-line outcomes and viewed-result identity. Neither opens cohort data or submits SCC jobs, and neither replaces the native `tests/test_single_bam_ui.R` suite.

The audit-table browser regression starts a synthetic fixture (`tests/audit_layout_app.R`) that uses the production `variant_ui()` and the production audit-table renderer with 75 synthetic rows and long synthetic paths. Chromium (Playwright) measures each DT output, its footer and the next heading at 10/25/50/100 rows, with pagination, horizontal scrolling, server-side search, tab switching, laptop/narrow viewports and 2x CSS zoom, first without `www/frontend.css` (reproduction control) and then with it. The CSS must keep overriding DT's `.html-fill-container > .html-fill-item.datatables { flex-basis: 400px; }` for the namespaced `samples` and `genotypes` outputs only; it must not impose fixed row heights or change global card/plot sizing. Long paths stay on one line and use horizontal scrolling.

```bash
python -m pip install playwright==1.58.0
python -m playwright install chromium
mkdir -p test-artifacts/audit-layout
Rscript tests/audit_layout_app.R > test-artifacts/audit-layout/shiny.log 2>&1 &
app_pid=$!
trap 'kill "$app_pid" 2>/dev/null || true' EXIT
python tests/test_audit_layout_browser.py || tail -n 100 test-artifacts/audit-layout/shiny.log
kill "$app_pid"; trap - EXIT
```

This mirrors the CI step: the fixture's log is kept and the fixture is stopped afterwards, so a later run cannot connect to a stale instance on port 8765.

Screenshots, geometry JSON and the fixture log are written to the ignored `test-artifacts/` directory (and uploaded as a short-lived CI artifact). Synthetic browser checks do not certify the SCC session's installed package versions.

Manual browser acceptance after a frontend change, in a real SCC RStudio session: drag the Single BAM display slider and verify the plot axes change without changing metrics; open and cancel the batch review; review batch A, modify inputs and verify nothing is submitted; confirm once and check independent job receipts; open saved job A while the form contains B; inspect the layout at laptop width and 200% zoom; verify modal focus and keyboard navigation.

## SCC acceptance

A release check should use a fresh clone/checkout in an independent writable path, clear personal `REGSHINY_*` deployment overrides, then prepare FHS from its original CSV and VCF/BAM sources. Validate a real coordinate lookup and bounded RNA preview, and also connect custom BAM/VCF data through explicit and read-group sample matching. Confirm that genotype groups, read-depth units, known-zero/missing states and data-source identity remain correct.

For scheduler runs, success requires a matching application receipt, artifact/source checks and terminal `qacct` with both `failed=0` and `exit_status=0`. Queue disappearance is not sufficient. Browser verification should use a real SCC RStudio session, check default FHS and a custom source, and save screenshots outside tracked source.

For parallel-job changes, submit at least three distinct exact variants from one prepared dataset. Verify separate scheduler job IDs and directories, `slots=16` in each terminal accounting record, and `allocated_slots=16` / `workers=15` in each result. Compare each result and selected sample count with its own frozen request. Check scheduler start/end times to distinguish observed overlap from merely submitting multiple queued jobs. Exercise invalid or ambiguous batch lines without cancelling valid lines, and verify that changing inputs or opening another job cannot relabel a result. A failed or interrupted job must not block unrelated submissions; resuming the same job still requires terminal accounting and an exclusive per-job lock.

Private run-specific receipts, manifests, data and screenshots belong in the user's validation/state location. Do not commit them or label a current checkout verified based only on an unrelated historical `validation_status.json` or stale checksum file.

Repository access, SCC source-file permissions, software correctness and biological interpretation are separate checks. Other members use their own authenticated GitHub/SCC accounts; an owner's successful run does not claim an independent member login was tested.

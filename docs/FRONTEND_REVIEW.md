# Frontend safety and navigation changes

This branch is based on main commit `a23833ffd5ae64155b4e964398ef6644294ccab9`.
It does not change RNA calculations, genotype grouping, the all-matched default,
per-job 16-core allocation, checkpoint rules, or source identity validation.

## Implemented

- Fix the Single BAM dynamic display slider's module namespace. Reject missing,
  non-finite, reversed, zero-width and out-of-range display windows. Display
  changes leave the fixed analysis interval and its denominator unchanged.
- Move batch submission and saved-job controls from the WGS analysis sidebar
  into the existing Jobs tab. Keep one instance of every input. The sidebar
  provides a direct navigation button; batch jobs remain full-cohort regardless
  of the single-variant preview setting.
- Add a side-effect-free resource review and explicit server-enforced batch
  confirmation. It shows the dataset, assembly, frozen filters and an upper
  bound of unique query strings times 16 cores. This is not an exact validated
  variant count: VCF resolution, allele disambiguation, record-level deduplication
  and source checks still run after confirmation. Prior batches are not deduplicated.
- Reject confirmation without a review, after inputs change, after cancellation,
  after source deactivation, and after an approval has already been consumed.
  Existing per-line submission outcomes and uncertain-receipt handling remain.
- Add a viewed-snapshot context panel, separate from the current query form.
  Keep full provenance and state codes available. Add polite status live regions,
  responsive action spacing and collapsible sample-accounting definitions.

## Verification

Run from the repository root with the project's R dependencies installed:

```sh
Rscript tests/test_frontend_contract.R
Rscript tests/test_variant_ui.R
```

The contract suite parses the changed R sources, renders the production WGS UI,
checks unique namespaced controls, and executes the production display-control
and display-range expressions in a native-free Shiny test module. The state suite
uses synthetic job mocks and exercises cancellation, stale confirmation, replay,
source deactivation, independent outcomes, and viewed-result identity. Neither
suite opens cohort data or submits SCC jobs.

Also run the existing native `tests/test_single_bam_ui.R` suite in the supported
SCC environment. Mock and structural tests do not replace native or browser tests.

Browser acceptance remains required: drag the Single BAM slider; verify the
plot axes change without changing metrics; open and cancel the batch review;
review A then modify inputs and verify no submission; confirm once and check
independent job receipts; open saved job A while the form contains B; inspect the
layout at laptop width and 200% zoom; verify modal focus and keyboard navigation.

## Deliberately deferred

Exact VCF preflight before resource confirmation, cross-batch duplicate detection,
shared genotype colors across all plots, linked junction highlighting, truncated
cross-window arcs, data-source form redesign and a top-level task-center page are
not implemented in this first frontend change. The existing Jobs tab is reorganized,
not replaced by a new global navigation architecture.

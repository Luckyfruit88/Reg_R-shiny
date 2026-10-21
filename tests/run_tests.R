# Run from the repository root: Rscript tests/run_tests.R [--unit]
# Single entry point for the synthetic R and Python suites. Each suite runs in
# its own process, so mocks and environment changes cannot leak between suites.
# --unit runs only checks that need no samtools, bcftools or RegTools (the CI
# subset). Without it, native sections run too: use an SCC compute node with
# samtools/1.23 regtools/1.0.0 bcftools/1.23 loaded.
# Browser layout checks need Playwright; see docs/VALIDATION.md.
unit_only <- "--unit" %in% commandArgs(trailingOnly = TRUE)
if (!file.exists("app.R") || !dir.exists("tests")) stop("Run from the repository root.")

# unit: arguments in --unit mode (NULL = skipped there); full: arguments otherwise.
suite <- function(script, unit = character(), full = unit) list(script = script, unit = unit, full = full)
suites <- list(
  suite("tests/test_backend.R", unit = "--unit", full = character()),
  suite("tests/test_variants.R", unit = "--unit", full = character()),
  suite("tests/test_all_samples.R", full = "--native"),
  suite("tests/test_variant_jobs.R"),
  suite("tests/test_parallel_variant_jobs.R"),
  suite("tests/test_annotations.R"),
  suite("tests/test_splice_evidence.R", full = "--native"),
  suite("tests/test_data_sources.R", full = "--native"),
  suite("tests/test_autoload.R", unit = "app.R"),
  suite("tests/test_frontend_contract.R"),
  suite("tests/test_variant_ui.R"),
  suite("tests/test_annotation_ui.R"),
  suite("tests/test_splice_evidence_ui.R"),
  suite("tests/test_data_sources_ui.R"),
  suite("tests/test_single_bam_ui.R", unit = NULL, full = character()),
  # Independent arithmetic oracle for the synthetic Single BAM fixture.
  suite("tests/oracle.py")
)

git_status <- function() {
  out <- tryCatch(suppressWarnings(system2("git", c("status", "--porcelain=v1", "--untracked-files=all"),
    stdout = TRUE, stderr = FALSE)), error = function(e) NULL)
  if (is.null(out) || !is.null(attr(out, "status"))) NULL else out
}

run_suite <- function(s) {
  args <- if (unit_only) s$unit else s$full
  if (is.null(args)) return(list(status = "SKIP", seconds = 0, log = NULL))
  command <- if (endsWith(s$script, ".py")) unname(Sys.which("python3")) else file.path(R.home("bin"), "Rscript")
  if (!nzchar(command)) return(list(status = "FAIL", seconds = 0, log = "python3 is not on PATH."))
  if (endsWith(s$script, ".py")) args <- c("-I", s$script, args) else args <- c(s$script, args)
  log <- tempfile("regshiny-test-", fileext = ".log")
  started <- Sys.time()
  code <- system2(command, args, stdout = log, stderr = log)
  list(status = if (identical(code, 0L)) "PASS" else "FAIL",
    seconds = as.numeric(difftime(Sys.time(), started, units = "secs")), log = readLines(log, warn = FALSE))
}

cat("Mode:", if (unit_only) "unit (no native tools)" else "full (native tools required)", "\n")
before <- git_status()
results <- lapply(suites, function(s) {
  r <- run_suite(s)
  cat(sprintf("%-4s %6.1fs  %s\n", r$status, r$seconds, s$script))
  if (r$status == "FAIL") cat(paste0("    | ", tail(r$log, 40L)), sep = "\n")
  r
})
statuses <- vapply(results, `[[`, character(1), "status")

# Tests must not leave files that git would offer to commit.
after <- git_status()
new_files <- if (is.null(before) || is.null(after)) character() else setdiff(after, before)
if (length(new_files)) cat("FAIL worktree: tests left new git status entries\n",
  paste0("    | ", new_files, collapse = "\n"), "\n", sep = "")

note <- if (is.null(before)) "; git worktree check unavailable" else if (length(new_files)) "" else "; worktree unchanged"
cat(sprintf("\n%d passed, %d failed, %d skipped%s\n", sum(statuses == "PASS"), sum(statuses == "FAIL"),
  sum(statuses == "SKIP"), note))
if (any(statuses == "FAIL") || length(new_files)) quit(status = 1L)

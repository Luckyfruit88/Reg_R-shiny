# Run from the repository root: Rscript tests/run_tests.R [--unit] [--require-clean]
# Single entry point for the synthetic R and Python suites. Each suite runs in
# its own process, so mocks and environment changes cannot leak between suites.
# --unit runs only checks that need no samtools, bcftools or RegTools (the CI
# subset). Without it, native sections run too: use an SCC compute node with
# samtools/1.23 regtools/1.0.0 bcftools/1.23 loaded.
# --require-clean (CI and release acceptance) fails unless the worktree is
# clean before and after the run. Otherwise a dirty worktree is allowed and the
# run fails only if the tests change it: new, removed or rewritten status entries.
# Full per-suite logs go to the ignored test-artifacts/test-runs/<run>/.
# Browser layout checks need Playwright; see docs/VALIDATION.md.

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
  suite("tests/test_run_tests.R"),
  # Independent arithmetic oracle for the synthetic Single BAM fixture.
  suite("tests/oracle.py"),
  # Independent arithmetic oracle for splice-site bin counts (shared JSON fixture).
  suite("tests/oracle_bin_ratio.py")
)

# Worktree snapshot: git status entries plus the content hash of every listed
# path, because an unchanged status code (" M", "??") says nothing about content.
# NULL means git could not report the state.
worktree_snapshot <- function(dir = ".") {
  status <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(dir), "status", "--porcelain=v1",
    "--untracked-files=all"), stdout = TRUE, stderr = FALSE)), error = function(e) NULL)
  if (is.null(status) || !is.null(attr(status, "status"))) return(NULL)
  paths <- sub("^.* -> ", "", substring(status, 4L))
  paths <- ifelse(startsWith(paths, "\""), gsub("^\"|\"$", "", paths), paths)
  hashes <- unname(tools::md5sum(file.path(dir, paths)))
  data.frame(status = substr(status, 1L, 2L), path = paths, md5 = hashes, stringsAsFactors = FALSE)
}

# Returns list(ok, message, changes). A failed git query is never a pass.
worktree_verdict <- function(before, after, require_clean = FALSE) {
  if (is.null(before) && require_clean)
    return(list(ok = FALSE, message = "worktree UNVERIFIED: git status unavailable before the run", changes = character()))
  if (is.null(before))
    return(list(ok = TRUE, message = "worktree not checked (not a git checkout)", changes = character()))
  if (require_clean && nrow(before))
    return(list(ok = FALSE, message = "worktree was not clean before the run (--require-clean)",
      changes = paste(before$status, before$path)))
  if (is.null(after))
    return(list(ok = FALSE, message = "worktree UNVERIFIED: git status failed after the run", changes = character()))
  key <- function(x) paste(x$status, x$path, x$md5, sep = "\t")
  changed <- union(setdiff(key(after), key(before)), setdiff(key(before), key(after)))
  paths <- unique(sub("\t.*$", "", sub("^[^\t]*\t", "", changed)))
  if (length(paths)) return(list(ok = FALSE, message = "tests changed the worktree", changes = paths))
  list(ok = TRUE, changes = character(), message = if (nrow(before))
    "worktree had local changes; the tests did not alter them" else "worktree clean before and after")
}

run_suite <- function(s, unit_only, log_dir) {
  args <- if (unit_only) s$unit else s$full
  if (is.null(args)) return(list(status = "SKIP", seconds = 0, log = NULL))
  python <- endsWith(s$script, ".py")
  command <- if (python) unname(Sys.which("python3")) else file.path(R.home("bin"), "Rscript")
  log <- file.path(log_dir, paste0(tools::file_path_sans_ext(basename(s$script)), ".log"))
  if (!nzchar(command)) {
    writeLines("python3 is not on PATH.", log)
    return(list(status = "FAIL", seconds = 0, log = log))
  }
  args <- c(if (python) "-I", s$script, args)
  started <- Sys.time()
  code <- system2(command, args, stdout = log, stderr = log)
  list(status = if (identical(code, 0L)) "PASS" else "FAIL",
    seconds = as.numeric(difftime(Sys.time(), started, units = "secs")), log = log)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  unit_only <- "--unit" %in% args
  require_clean <- "--require-clean" %in% args
  if (!file.exists("app.R") || !dir.exists("tests")) stop("Run from the repository root.")
  log_dir <- file.path("test-artifacts", "test-runs", format(Sys.time(), "%Y%m%d-%H%M%S"))
  before <- worktree_snapshot()
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Mode:", if (unit_only) "unit (no native tools)" else "full (native tools required)",
    "\nLogs:", log_dir, "\n")
  results <- lapply(suites, function(s) {
    r <- run_suite(s, unit_only, log_dir)
    cat(sprintf("%-4s %6.1fs  %s\n", r$status, r$seconds, s$script))
    if (r$status == "FAIL") {
      cat(paste0("    | ", tail(readLines(r$log, warn = FALSE), 40L)), sep = "\n")
      cat("    full log:", r$log, "\n")
    }
    r
  })
  statuses <- vapply(results, `[[`, character(1), "status")
  tree <- worktree_verdict(before, worktree_snapshot(), require_clean)
  cat(if (tree$ok) "" else "FAIL ", tree$message, "\n", sep = "")
  if (length(tree$changes)) cat(paste0("    | ", tree$changes), sep = "\n")
  cat(sprintf("\n%d passed, %d failed, %d skipped\n", sum(statuses == "PASS"), sum(statuses == "FAIL"),
    sum(statuses == "SKIP")))
  if (any(statuses == "FAIL") || !tree$ok) quit(status = 1L)
}

if (sys.nframe() == 0L) main()

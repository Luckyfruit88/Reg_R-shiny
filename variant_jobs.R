# Persistent, per-user SCC jobs for all-matched-BAM Reg_Shiny analyses.
# Source backend.R and variant_backend.R first. UI probes never block on qacct.
.variant_job_probes <- new.env(parent = emptyenv())
.variant_job_uid <- NULL

job_uid <- function() {
  if (!is.null(.variant_job_uid)) return(.variant_job_uid)
  p <- processx::run("id", "-u", timeout = 5)
  x <- trimws(p$stdout)
  if (!grepl("^[0-9]+$", x)) stop("Cannot determine the current Unix user.")
  .variant_job_uid <<- x
  x
}
job_atomic_json <- function(value, path) {
  tmp <- tempfile(".publish-", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  jsonlite::write_json(value, tmp, auto_unbox = TRUE, pretty = TRUE, na = "null")
  Sys.chmod(tmp, "0600")
  if (!file.rename(tmp, path)) stop("Cannot publish job state.")
}
job_read_json <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(jsonlite::read_json(path, simplifyVector = TRUE), error = function(e) NULL)
}
job_now <- function() format(Sys.time(), tz = "UTC", usetz = TRUE)
job_private_dir <- function(path, create = FALSE) {
  if (isTRUE(!is.na(Sys.readlink(path)) && nzchar(Sys.readlink(path)))) stop("Job directories must not be symbolic links.")
  if (create && !dir.exists(path)) dir.create(path, recursive = FALSE, mode = "0700")
  if (isTRUE(!is.na(Sys.readlink(path)) && nzchar(Sys.readlink(path)))) stop("Job directories must not be symbolic links.")
  if (!dir.exists(path)) stop("Job directory is unavailable.")
  info <- file.info(path)
  if (is.null(info$uid) || as.character(info$uid) != job_uid()) stop("Job directory belongs to another user.")
  if (bitwAnd(as.integer(info$mode), 63L) != 0L) stop("Job directory must be private (0700).")
  normalizePath(path, mustWork = TRUE)
}
variant_user_job_root <- function(job_root = Sys.getenv("REGSHINY_JOB_ROOT")) {
  if (length(job_root) != 1L || !nzchar(job_root) || !startsWith(job_root, "/") || !dir.exists(job_root))
    stop("REGSHINY_JOB_ROOT must name the configured persistent project job directory.")
  base <- normalizePath(job_root, mustWork = TRUE)
  mode <- as.integer(file.info(base)$mode)
  if (bitwAnd(mode, 16L) != 0L && bitwAnd(mode, 512L) == 0L)
    stop("A shared group-writable job root must have the sticky bit (03770).")
  job_private_dir(file.path(base, paste0("uid-", job_uid())), create = TRUE)
}
validate_variant_job_dir <- function(job_dir, job_root = Sys.getenv("REGSHINY_JOB_ROOT")) {
  root <- variant_user_job_root(job_root)
  if (length(job_dir) != 1L || !startsWith(job_dir, paste0(root, "/")) ||
      dirname(job_dir) != root || !grepl("^run-[A-Za-z0-9_-]+$", basename(job_dir)))
    stop("Select a job from your own Reg_Shiny job history.")
  job_private_dir(job_dir)
}
job_request <- function(job_dir) {
  x <- readRDS(file.path(job_dir, "request.rds"))
  if (!is.list(x) || !identical(x$schema, "regshiny-full-job-v1")) stop("Invalid saved job request.")
  x
}
job_result_valid <- function(job_dir) {
  receipt <- job_read_json(file.path(job_dir, "application_receipt.json"))
  path <- file.path(job_dir, "result.rds")
  state <- job_read_json(file.path(job_dir, "state.json"))
  !is.null(receipt) && !is.null(state) && identical(as.character(receipt$job_id), as.character(state$job_id)) &&
    identical(receipt$status, "APPLICATION_COMPLETE") && file.exists(path) &&
    identical(unname(tools::md5sum(path)), receipt$result_md5)
}
job_parse_accounting <- function(lines, expected_job_id) {
  value <- function(field) {
    hit <- lines[grepl(paste0("^", field, "[[:space:]]+"), lines)]
    if (length(hit) != 1L) return(NA_character_)
    strsplit(trimws(hit), "[[:space:]]+")[[1L]][[2L]]
  }
  if (!identical(value("jobnumber"), as.character(expected_job_id))) return(NULL)
  failed <- suppressWarnings(as.integer(value("failed")))
  exit_status <- suppressWarnings(as.integer(value("exit_status")))
  if (is.na(failed) || is.na(exit_status)) return(NULL)
  list(job_id = as.character(expected_job_id), failed = failed, exit_status = exit_status)
}
# Polling commands run asynchronously and are throttled. Closing the UI may stop
# a harmless probe, never the submitted qsub analysis.
job_scheduler_probe <- function(job_dir, state) {
  if (is.null(state$job_id) || !grepl("^[0-9]+$", as.character(state$job_id))) return(state)
  key <- paste(job_dir, state$job_id, sep = "|")
  now <- as.numeric(Sys.time())
  probe <- if (exists(key, .variant_job_probes, inherits = FALSE)) get(key, .variant_job_probes) else NULL
  if (!is.null(probe) && !is.null(probe$process)) {
    if (probe$process$is_alive()) {
      if (now - probe$started > 45) {
        probe$process$kill_tree(); probe$process <- NULL; probe$last_checked <- now
        assign(key, probe, .variant_job_probes)
      }
      return(state)
    }
    exit <- probe$process$get_exit_status()
    output <- if (file.exists(probe$stdout)) readLines(probe$stdout, warn = FALSE) else character()
    if (probe$kind == "qstat") {
      if (identical(exit, 0L)) {
        if (!identical(state$status, "FAILED")) state$status <- if (identical(state$status, "RUNNING") || file.exists(file.path(job_dir, "progress.json"))) "RUNNING" else "QUEUED"
        probe$kind <- "qstat"; probe$last_checked <- now
      } else {
        probe$kind <- "qacct"; probe$last_checked <- 0
        if (!identical(state$status, "FAILED")) state$status <- "ACCOUNTING_PENDING"
      }
    } else {
      accounting <- job_parse_accounting(output, state$job_id)
      if (!is.null(accounting)) {
        state$accounting <- accounting
        if (accounting$failed == 0L && accounting$exit_status == 0L && job_result_valid(job_dir)) {
          receipt <- job_read_json(file.path(job_dir, "application_receipt.json"))
          state$status <- if (receipt$failed_samples > 0) "COMPLETE_WITH_FAILURES" else "COMPLETE"
          state$result_ready <- TRUE
        } else {
          state$status <- if (file.exists(file.path(job_dir, "error.json"))) "FAILED" else "INTERRUPTED"
          state$result_ready <- FALSE
          state$error <- "Scheduler termination or missing/invalid application receipt. Resume the saved request after resolving the cause."
        }
        state$verified_at <- job_now()
      } else {
        # A qstat error may be transient. Recheck liveness instead of remaining
        # on qacct for the entire duration of a still-running job.
        probe$kind <- "qstat"
        state$scheduler_note <- "Terminal accounting is not available yet; job liveness will be checked again."
      }
      probe$last_checked <- now
    }
    probe$process <- NULL
    job_atomic_json(state, file.path(job_dir, "state.json"))
    assign(key, probe, .variant_job_probes)
  }
  if (isTRUE(state$result_ready) || !is.null(state$accounting)) return(state)
  if (is.null(probe)) probe <- list(kind = "qstat", last_checked = 0, process = NULL)
  interval <- if (probe$kind == "qacct") 60 else 30
  if (now - probe$last_checked < interval) return(state)
  executable <- unname(Sys.which(probe$kind))
  if (!nzchar(executable)) return(state)
  probe$stdout <- file.path(job_dir, paste0(probe$kind, "-", state$job_id, ".txt"))
  probe$stderr <- file.path(job_dir, paste0(probe$kind, "-", state$job_id, ".err"))
  probe$started <- now
  probe$process <- processx::process$new(executable, c("-j", as.character(state$job_id)),
    stdout = probe$stdout, stderr = probe$stderr, cleanup = TRUE, cleanup_tree = TRUE)
  assign(key, probe, .variant_job_probes)
  state
}
read_variant_job <- function(job_dir, job_root = Sys.getenv("REGSHINY_JOB_ROOT")) {
  job_dir <- validate_variant_job_dir(job_dir, job_root)
  request <- job_request(job_dir)
  state <- job_read_json(file.path(job_dir, "state.json"))
  if (is.null(state)) state <- list(status = "FAILED", result_ready = FALSE, error = "Submission state is missing.")
  state <- job_scheduler_probe(job_dir, state)
  progress <- job_read_json(file.path(job_dir, "progress.json"))
  if (is.null(progress)) progress <- list(total = nrow(request$selected_plan), completed = 0, success = 0, failed = 0, resumed = 0)
  progress$pending <- max(0, progress$total - progress$completed)
  error <- job_read_json(file.path(job_dir, "error.json"))
  list(job_dir = job_dir, job_id = state$job_id, status = state$status,
    progress = progress, result_ready = isTRUE(state$result_ready),
    error = if (!is.null(error)) error$message else state$error,
    scheduler_terminal = !is.null(state$accounting), accounting = state$accounting,
    request = list(variant = request$variant, cfg = request$cfg, ui_controls = request$ui_controls,
                   created_at = request$created_at, selected_n = nrow(request$selected_plan)))
}
list_variant_jobs <- function(job_root = Sys.getenv("REGSHINY_JOB_ROOT")) {
  root <- variant_user_job_root(job_root)
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  rows <- lapply(dirs, function(d) tryCatch({
    d <- validate_variant_job_dir(d, job_root); x <- job_request(d)
    state <- job_read_json(file.path(d, "state.json"))
    data.frame(job_dir = d, label = paste(x$created_at, x$variant$record_id,
      if (is.null(state)) "UNKNOWN" else state$status, sep = " | "), stringsAsFactors = FALSE)
  }, error = function(e) NULL))
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame(job_dir = character(), label = character()))
  out <- do.call(rbind, rows); out[order(out$job_dir, decreasing = TRUE), , drop = FALSE]
}
job_assert_no_active <- function(job_root, except = NULL) {
  jobs <- list_variant_jobs(job_root)
  for (d in setdiff(jobs$job_dir, except)) {
    state <- read_variant_job(d, job_root)
    if (!state$status %in% c("COMPLETE", "COMPLETE_WITH_FAILURES", "FAILED", "INTERRUPTED", "SUBMISSION_REJECTED") ||
        (state$status %in% c("FAILED", "INTERRUPTED") && !isTRUE(state$scheduler_terminal)))
      stop("One full-cohort job is already active or awaiting scheduler verification. Open it from job history before submitting another.")
  }
}
job_submit_saved <- function(job_dir, job_root) {
  request <- job_request(job_dir)
  qsub <- unname(Sys.which("qsub"))
  if (!nzchar(qsub)) stop("qsub is unavailable. Start Reg_Shiny inside SCC with the documented modules.")
  number <- length(list.files(job_dir, pattern = "^submission-[0-9]+\\.json$")) + 1L
  state <- list(status = "SUBMITTING", result_ready = FALSE, attempt = number, submitted_at = job_now())
  job_atomic_json(state, file.path(job_dir, "state.json"))
  p <- tryCatch(processx::run(qsub, c("-terse", "-P", request$project, "-pe", "omp", request$cores,
      "-l", paste0("h_rt=", request$walltime), "-l", "mem_per_core=4G", "-N", "RegShiny_all",
      "-o", file.path(job_dir, paste0("attempt-", number, ".out")),
      "-e", file.path(job_dir, paste0("attempt-", number, ".err")), file.path(job_dir, "run.sh")),
      timeout = 30, error_on_status = FALSE), error = function(e) e)
  if (!inherits(p, "error")) {
    writeLines(p$stdout, file.path(job_dir, paste0("submission-", number, "-qsub.out")))
    writeLines(p$stderr, file.path(job_dir, paste0("submission-", number, "-qsub.err")))
    Sys.chmod(file.path(job_dir, paste0("submission-", number, c("-qsub.out", "-qsub.err"))), "0600")
  }
  if (inherits(p, "error") || p$status != 0L) {
    state$status <- if (inherits(p, "error")) "SUBMISSION_UNKNOWN" else "SUBMISSION_REJECTED"
    state$error <- if (inherits(p, "error")) paste("Submission outcome is unknown; inspect qstat before any retry.", conditionMessage(p)) else p$stderr
    job_atomic_json(state, file.path(job_dir, "state.json")); stop(state$error)
  }
  id <- trimws(p$stdout)
  if (!grepl("^[0-9]+$", id)) {
    state$status <- "SUBMISSION_UNKNOWN"; state$error <- "qsub returned an unrecognized job receipt; inspect this request before resubmitting."
    job_atomic_json(state, file.path(job_dir, "state.json")); stop(state$error)
  }
  state$job_id <- id; state$status <- "QUEUED"
  job_atomic_json(state, file.path(job_dir, paste0("submission-", number, ".json")))
  job_atomic_json(state, file.path(job_dir, "state.json"))
  list(job_dir = job_dir, job_id = id, status = "QUEUED")
}
job_acquire_submit_lock <- function(job_root) {
  root <- variant_user_job_root(job_root)
  lock <- file.path(root, ".submit-lock")
  if (!dir.create(lock, mode = "0700", showWarnings = FALSE))
    stop("A submission is already in progress for this account. If its process was interrupted, an administrator must verify it ended before clearing .submit-lock.")
  job_atomic_json(list(pid = Sys.getpid(), host = Sys.info()[["nodename"]], at = job_now()), file.path(lock, "owner.json"))
  lock
}
submit_variant_job <- function(resources, variant, cfg, backend_file, variant_file, job_source,
                              job_root = Sys.getenv("REGSHINY_JOB_ROOT"), ui_controls = NULL) {
  root <- variant_user_job_root(job_root)
  lock <- job_acquire_submit_lock(job_root); on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  job_assert_no_active(job_root)
  variant_check_resources(resources)
  cfg <- validate_config(cfg)
  if (cfg$demo || nrow(variant) != 1L || cfg$chrom != variant$chrom ||
      variant$pos1 < cfg$start1 || variant$pos1 > cfg$end1) stop("Select one real variant inside the requested interval.")
  current <- lookup_variants(resources, paste(variant$chrom, variant$pos1, variant$ref, variant$alt, sep = ":"))
  if (nrow(current) != 1L || current$record_id != variant$record_id) stop("Variant identity changed; repeat lookup.")
  info <- variant_vcf_info(current$vcf, current$chrom)
  calls <- variant_get_genotypes(resources, current)
  plan <- calls[calls$linked & calls$call_status == "CALLED", c("vcf_sample", "raw_gt", "genotype", "bam"), drop = FALSE]
  plan <- plan[order(plan$genotype, plan$vcf_sample, method = "radix"), , drop = FALSE]; rownames(plan) <- NULL
  if (!nrow(plan)) stop("No called DNA genotype has a matched RNA BAM.")
  cores <- suppressWarnings(as.integer(Sys.getenv("REGSHINY_FULL_CORES", "8")))
  cores <- int_scalar(cores, "Full-job cores", 1, 16)
  walltime <- Sys.getenv("REGSHINY_FULL_WALLTIME", "12:00:00")
  if (!grepl("^([0-9]|[1-4][0-9]):[0-5][0-9]:[0-5][0-9]$", walltime)) stop("Invalid full-job walltime.")
  project <- Sys.getenv("REGSHINY_SGE_PROJECT", "mtdna-alcohol")
  if (!grepl("^[A-Za-z0-9_-]+$", project)) stop("Invalid scheduler project.")
  job_dir <- tempfile(paste0("run-", format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), "-"), tmpdir = root)
  job_private_dir(job_dir, create = TRUE)
  dir.create(file.path(job_dir, "source"), mode = "0700")
  sources <- setNames(c(backend_file, variant_file, job_source), c("backend.R", "variant_backend.R", "variant_jobs.R"))
  for (name in names(sources)) {
    path <- file.path(job_dir, "source", name)
    if (!file.copy(sources[[name]], path)) stop("Cannot snapshot application source.")
    Sys.chmod(path, "0400")
  }
  hashes <- unname(tools::md5sum(file.path(job_dir, "source", names(sources)))); names(hashes) <- names(sources)
  cfg$source_fingerprint <- paste(hashes, collapse = ":")
  request <- list(schema = "regshiny-full-job-v1", created_at = job_now(), owner_uid = job_uid(),
    resources = resources, variant = current, cfg = cfg, ui_controls = ui_controls,
    selected_plan = plan, vcf_state = info$state, bam_states = lapply(plan$bam, variant_bam_state), source_hashes = hashes,
    cores = cores, workers = max(1L, cores-1L), walltime = walltime, project = project)
  saveRDS(request, file.path(job_dir, "request.rds")); Sys.chmod(file.path(job_dir, "request.rds"), "0400")
  entry <- c('args <- commandArgs(trailingOnly=TRUE)', 'job_dir <- args[[1L]]',
    'for (f in c("backend.R","variant_backend.R","variant_jobs.R")) source(file.path(job_dir,"source",f))',
    'run_variant_job(job_dir)')
  writeLines(entry, file.path(job_dir, "worker.R")); Sys.chmod(file.path(job_dir, "worker.R"), "0400")
  script <- c("#!/bin/bash -l", "set -euo pipefail", "umask 077",
    "module load R/4.5.2 samtools/1.23 regtools/1.0.0 bcftools/1.23",
    "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1",
    paste0("export REGSHINY_JOB_ROOT=", shQuote(normalizePath(job_root))),
    paste("exec Rscript --vanilla", shQuote(file.path(job_dir, "worker.R")), shQuote(job_dir)))
  writeLines(script, file.path(job_dir, "run.sh")); Sys.chmod(file.path(job_dir, "run.sh"), "0500")
  job_submit_saved(job_dir, job_root)
}
resume_variant_job <- function(job_dir, job_root = Sys.getenv("REGSHINY_JOB_ROOT")) {
  job_dir <- validate_variant_job_dir(job_dir, job_root)
  lock <- job_acquire_submit_lock(job_root); on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  state <- read_variant_job(job_dir, job_root)
  if (!state$status %in% c("FAILED", "INTERRUPTED") || !isTRUE(state$scheduler_terminal))
    stop("Only a failed/interrupted job with terminal scheduler accounting can be resumed.")
  job_assert_no_active(job_root, except = job_dir)
  for (name in c("error.json", "application_receipt.json", "result.rds", "progress.json")) {
    path <- file.path(job_dir, name)
    if (file.exists(path)) file.rename(path, paste0(path, ".prior-", format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")))
  }
  job_submit_saved(job_dir, job_root)
}
load_variant_job_result <- function(job_dir, job_root = Sys.getenv("REGSHINY_JOB_ROOT")) {
  state <- read_variant_job(job_dir, job_root)
  if (!isTRUE(state$result_ready) || !job_result_valid(state$job_dir)) stop("The full result is not yet verified and ready.")
  readRDS(file.path(state$job_dir, "result.rds"))
}
run_variant_job <- function(job_dir) {
  job_dir <- validate_variant_job_dir(job_dir)
  tryCatch({
    request <- job_request(job_dir)
    actual <- unname(tools::md5sum(file.path(job_dir, "source", names(request$source_hashes))))
    if (!identical(actual, unname(request$source_hashes))) stop("Saved source snapshot changed.")
    variant_check_resources(request$resources)
    if (!identical(request$vcf_state, variant_file_state(request$vcf_state$path))) stop("VCF or index changed after submission.")
    calls <- variant_get_genotypes(request$resources, request$variant)
    selected <- calls[calls$linked & calls$call_status == "CALLED", c("vcf_sample", "raw_gt", "genotype", "bam"), drop = FALSE]
    selected <- selected[order(selected$genotype, selected$vcf_sample, method = "radix"), , drop = FALSE]; rownames(selected) <- NULL
    if (!identical(selected, request$selected_plan)) stop("Matched genotype/BAM selection changed after submission.")
    if (!identical(request$bam_states, lapply(selected$bam, variant_bam_state)))
      stop("A BAM or index changed while the request was queued. Submit a new request for the new inputs.")
    slots <- suppressWarnings(as.integer(Sys.getenv("NSLOTS", "1")))
    if (is.na(slots) || slots < request$cores) stop("Allocated slots are smaller than the frozen request.")
    current_job_id <- Sys.getenv("JOB_ID")
    if (!grepl("^[0-9]+$", current_job_id)) stop("The worker must run inside its SCC scheduler allocation.")
    state <- job_read_json(file.path(job_dir, "state.json")); state$status <- "RUNNING"; state$job_id <- current_job_id
    job_atomic_json(state, file.path(job_dir, "state.json"))
    result <- analyze_variant(request$resources, request$variant, request$cfg,
      sampling_mode = "all", progress_file = file.path(job_dir, "progress.json"),
      checkpoint_dir = file.path(job_dir, "checkpoints"), resume = TRUE, workers = request$workers)
    if (sum(result$group_summary$selected_n) != nrow(request$selected_plan) ||
        sum(result$group_summary$analyzed_n + result$group_summary$failed_n) != nrow(request$selected_plan))
      stop("Full-cohort accounting is incomplete.")
    tmp <- file.path(job_dir, ".result.rds.partial"); saveRDS(result, tmp); Sys.chmod(tmp, "0600")
    if (!file.rename(tmp, file.path(job_dir, "result.rds"))) stop("Cannot publish complete result.")
    job_atomic_json(list(status = "APPLICATION_COMPLETE", job_id = current_job_id, selected_n = nrow(request$selected_plan),
      successful_samples = sum(result$group_summary$analyzed_n), failed_samples = sum(result$group_summary$failed_n),
      result_md5 = unname(tools::md5sum(file.path(job_dir, "result.rds"))), completed_at = job_now()),
      file.path(job_dir, "application_receipt.json"))
    cat("FULL_APPLICATION_COMPLETE selected=", nrow(request$selected_plan),
        " success=", sum(result$group_summary$analyzed_n), " failed=", sum(result$group_summary$failed_n), "\n", sep = "")
  }, error = function(e) {
    job_atomic_json(list(message = conditionMessage(e), at = job_now()), file.path(job_dir, "error.json"))
    state <- job_read_json(file.path(job_dir, "state.json")); state$status <- "FAILED"; state$result_ready <- FALSE
    job_atomic_json(state, file.path(job_dir, "state.json"))
    stop(e)
  })
  invisible(TRUE)
}

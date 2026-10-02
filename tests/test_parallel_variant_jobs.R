# Independent regression: an old scheduler probe must not overwrite a resumed
# attempt published by another RStudio process. Synthetic state only; no SGE.
source("variant_jobs.R")
work <- tempfile("regshiny-parallel-probe-")
dir.create(work, mode = "0700")

exercise_stale_probe <- function(kind) {
  job_dir <- file.path(work, paste0("run-", kind))
  dir.create(job_dir, mode = "0700")
  current <- list(status = "RUNNING", result_ready = FALSE, attempt = 2L,
                  job_id = "222", attempt_token = "new-attempt")
  stale <- list(status = "RUNNING", result_ready = FALSE, attempt = 1L,
                job_id = "111", attempt_token = "old-attempt")
  job_atomic_json(current, file.path(job_dir, "state.json"))
  output <- file.path(work, paste0(kind, "-old.txt"))
  writeLines(if (kind == "qacct") c("jobnumber 111", "failed 0", "exit_status 1") else
               "job_number: 111", output)
  process <- list(is_alive = function() FALSE, get_exit_status = function() 0L,
                  kill = function(...) invisible(NULL), kill_tree = function(...) invisible(NULL))
  assign(paste(job_dir, stale$job_id, sep = "|"),
         list(kind = kind, process = process, stdout = output,
              started = as.numeric(Sys.time()) - 1, last_checked = 0), .variant_job_probes)
  answer <- job_scheduler_probe(job_dir, stale)
  disk <- job_read_json(file.path(job_dir, "state.json"))
  if (!identical(as.character(disk$job_id), "222") ||
      !identical(as.character(disk$attempt_token), "new-attempt") ||
      !identical(disk$status, "RUNNING")) {
    stop("A stale ", kind, " probe overwrote the newer scheduler attempt.", call. = FALSE)
  }
  if (!identical(as.character(answer$job_id), "222"))
    stop("A stale ", kind, " probe returned old-attempt state to the UI.", call. = FALSE)
}

tryCatch({
  exercise_stale_probe("qstat")
  exercise_stale_probe("qacct")
  cat("Cross-attempt scheduler probe publication isolation: PASS\n")
}, finally = unlink(work, recursive = TRUE))

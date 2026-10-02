# Synthetic scheduler tests. No SCC submission or study data are used.
source("backend.R"); source("variant_backend.R"); source("variant_jobs.R")
expect_error <- function(expr, pattern = NULL) {
  error <- tryCatch({ force(expr); NULL }, error = function(e) e)
  stopifnot(inherits(error, "error"))
  if (!is.null(pattern)) stopifnot(grepl(pattern, conditionMessage(error), fixed = TRUE))
  invisible(error)
}
local({
  root <- tempfile("regshiny_jobs_"); dir.create(root, mode = "0700")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  bin <- file.path(root, "bin"); dir.create(bin)
  mode <- file.path(root, "scheduler_mode"); writeLines("live", mode)
  submit_mode <- file.path(root, "submission_mode"); writeLines("success", submit_mode)
  counter <- file.path(root, "counter"); writeLines("10000", counter)
  script <- function(name, body, interpreter = "/bin/sh") {
    path <- file.path(bin, name); writeLines(c(paste0("#!", interpreter), body), path); Sys.chmod(path, "0700")
  }
  script("qsub", c("import sys,fcntl", "from pathlib import Path",
    paste0("counter=Path(", encodeString(counter, quote = "'"), ")"),
    paste0("mode=Path(", encodeString(submit_mode, quote = "'"), ").read_text().strip()"),
    "with counter.open('r+') as f:", " fcntl.flock(f,fcntl.LOCK_EX)",
    " n=int(f.read().strip())+1", " f.seek(0);f.truncate();f.write(str(n)+'\\n');f.flush()",
    "Path(str(counter)+'-'+str(n)+'.args').write_text('\\n'.join(sys.argv[1:])+'\\n')",
    "if mode=='reject': sys.exit('synthetic scheduler rejection')",
    "print(n if mode=='success' else 'ambiguous submission receipt')"), unname(Sys.which("python3")))
  script("qstat", paste0("if [ \"$(cat ", shQuote(mode), ")\" = live ]; then exit 0; else exit 1; fi"))
  script("qacct", c(paste0("case \"$(cat ", shQuote(mode), ")\" in"),
    "success) printf 'jobnumber %s\\nfailed 0\\nexit_status 0\\n' \"$2\";;",
    "failure) printf 'jobnumber %s\\nfailed 0\\nexit_status 1\\n' \"$2\";;", "*) exit 1;; esac"))
  variables <- c("PATH", "REGSHINY_JOB_ROOT", "REGSHINY_FULL_CORES", "JOB_ID", "REGSHINY_ATTEMPT_TOKEN", "NSLOTS")
  old <- Sys.getenv(variables, unset = NA_character_)
  on.exit({ Sys.unsetenv(variables[is.na(old)]); for (name in names(old)[!is.na(old)]) do.call(Sys.setenv, setNames(list(old[[name]]), name)) }, add = TRUE)
  Sys.setenv(PATH = paste(bin, Sys.getenv("PATH"), sep = .Platform$path.sep), REGSHINY_JOB_ROOT = root,
             REGSHINY_FULL_CORES = "8", NSLOTS = "16")
  uidroot <- variant_user_job_root(root)
  outside <- file.path(root, "outside"); dir.create(outside, mode = "0700")
  link <- file.path(uidroot, "run-link"); file.symlink(outside, link)
  expect_error(validate_variant_job_dir(link, root), "symbolic")
  expect_error(validate_variant_job_dir(outside, root), "own Reg_Shiny")
  vcf <- file.path(root, "synthetic.vcf.gz"); writeLines("synthetic", vcf); writeLines("index", paste0(vcf, ".tbi"))
  bams <- file.path(root, paste0("sample", seq_len(36), ".bam"))
  for (path in bams) { writeLines("synthetic", path); writeLines("index", paste0(path, ".bai")) }
  registry <- file.path(root, "vcf.tsv"); manifest <- file.path(root, "samples.tsv")
  write.table(data.frame(chrom = "chrDemo", vcf = vcf, build = "SYNTHETIC"), registry, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(data.frame(vcf_sample = paste0("S", seq_len(36)), bam = bams, build = "SYNTHETIC"), manifest, sep = "\t", row.names = FALSE, quote = FALSE)
  resources <- load_variant_resources(registry, manifest)
  variants <- data.frame(record_id = paste0("SYNTHETIC:chrDemo:", 190:193, ":A:G"), chrom = "chrDemo", pos1 = 190:193,
    ref = "A", alt = "G", build = "SYNTHETIC", vcf = vcf)
  calls <- data.frame(vcf_sample = paste0("S", seq_len(36)), raw_gt = rep(c("0/0", "0/1", "1/1"), each = 12),
    genotype = rep(c("0/0", "0/1", "1/1"), each = 12), bam = bams, linked = TRUE, call_status = "CALLED")
  replaced <- c("lookup_variants", "variant_vcf_info", "variant_get_genotypes", "analyze_variant")
  originals <- mget(replaced, .GlobalEnv)
  on.exit(list2env(originals, .GlobalEnv), add = TRUE)
  assign("lookup_variants", function(resources, query) {
    p <- as.integer(strsplit(query, ":", fixed = TRUE)[[1L]][[2L]])
    variants[variants$pos1 == p, , drop = FALSE]
  }, .GlobalEnv)
  assign("variant_vcf_info", function(path, chrom) list(state = variant_file_state(c(path, paste0(path, ".tbi")))), .GlobalEnv)
  assign("variant_get_genotypes", function(resources, variant) {
    if (variant$pos1 == 193L) stop("Synthetic matching failure for this variant")
    if (variant$pos1 == 191L) calls[1:18, , drop = FALSE] else calls
  }, .GlobalEnv)
  cfg <- default_demo_config(); cfg$demo <- FALSE; cfg$data_source_id <- "synthetic-profile"
  submit <- function(i) submit_variant_job(resources, variants[i, , drop = FALSE], cfg,
    normalizePath("backend.R"), normalizePath("variant_backend.R"), normalizePath("variant_jobs.R"), root,
    ui_controls = list(sampling_mode = "all"))
  first <- submit(1L); second <- submit(2L)
  request <- job_request(first$job_dir); request2 <- job_request(second$job_dir)
  stopifnot(first$job_id != second$job_id, first$job_dir != second$job_dir,
    request$cores == 16L, request$workers == 15L, request2$cores == 16L,
    nrow(request$selected_plan) == 36L, nrow(request2$selected_plan) == 18L,
    length(request$bam_states) == 36L, length(request$source_hashes) == 3L,
    !identical(first$attempt_token, second$attempt_token))
  args <- readLines(paste0(counter, "-", first$job_id, ".args"))
  pe <- match("-pe", args); stopifnot(identical(args[(pe + 1L):(pe + 2L)], c("omp", "16")), "-v" %in% args)
  expect_error(job_submit_saved(first$job_dir, root), "already submitted")
  # A busy lock for one job never excludes a different variant's submission.
  lock <- job_acquire_submit_lock(first$job_dir)
  expect_error(job_acquire_submit_lock(first$job_dir), "already in progress")
  third <- submit(3L); unlink(lock, recursive = TRUE)
  stopifnot(third$job_dir != first$job_dir)
  failed_match <- expect_error(submit(4L), "matching failure")
  stopifnot(failed_match$status == "PREPARATION_FAILED", dir.exists(failed_match$job_dir),
    file.exists(file.path(failed_match$job_dir, "preparation_error.json")),
    identical(job_request(failed_match$job_dir)$prepared, FALSE))
  expect_error(job_submit_saved(failed_match$job_dir, root), "already submitted")
  dashboard <- list_variant_job_statuses(root)
  stopifnot(nrow(dashboard) == 4L, all(dashboard$cores == 16L), all(dashboard$workers == 15L),
    dashboard$total[dashboard$job_dir == second$job_dir] == 18L,
    is.na(dashboard$total[dashboard$job_dir == failed_match$job_dir]))
  # Immutable metadata caches are compact: no genotype/BAM plan is retained.
  cached <- get(first$job_dir, .variant_job_metadata)
  stopifnot(!any(c("selected_plan", "bam_states", "resources") %in% names(cached$value)))
  logs <- read_variant_job_logs(failed_match$job_dir, root, max_lines = 2L)
  stopifnot(logs$logs[["preparation_error.json"]]$exists,
    length(logs$logs[["preparation_error.json"]]$lines) <= 2L,
    all(vapply(logs$logs, function(x) startsWith(x$path, paste0(failed_match$job_dir, "/")), logical(1))))
  # Per-attempt worker output/progress and source/selection validation.
  execution <- new.env(); execution$last <- NULL
  assign("analyze_variant", function(resources, variant, cfg, sampling_mode, progress_file, checkpoint_dir, resume, workers) {
    stopifnot(sampling_mode == "all", workers == 15L, resume)
    n <- if (variant$pos1 == 191L) 18L else 36L
    execution$last <- list(progress = progress_file, checkpoints = checkpoint_dir)
    job_atomic_json(list(state = "COMPLETE", total = n, completed = n, success = n, failed = 0L, resumed = 0L), progress_file)
    list(synthetic = TRUE, group_summary = data.frame(selected_n = n, analyzed_n = n, failed_n = 0L),
      provenance = list(execution = list(workers = workers, allocated_slots = as.integer(Sys.getenv("NSLOTS")))))
  }, .GlobalEnv)
  Sys.setenv(JOB_ID = first$job_id, REGSHINY_ATTEMPT_TOKEN = first$attempt_token)
  uncertain <- job_read_json(file.path(first$job_dir, "state.json"))
  uncertain$status <- "SUBMISSION_UNKNOWN"; uncertain$error <- "Submission outcome unknown"; uncertain$job_id <- NULL
  job_atomic_json(uncertain, file.path(first$job_dir, "state.json"))
  run_variant_job(first$job_dir)
  stopifnot(grepl(first$attempt_token, basename(execution$last$progress), fixed = TRUE),
    !file.exists(file.path(first$job_dir, "progress.json")), job_result_valid(first$job_dir),
    job_read_json(file.path(first$job_dir, "application_receipt.json"))$attempt_token == first$attempt_token,
    !read_variant_job(first$job_dir, root)$result_ready,
    is.null(job_read_json(file.path(first$job_dir, "state.json"))$error))
  drive <- function(directory, expected) {
    for (i in seq_len(150L)) {
      for (key in ls(.variant_job_probes)) { p <- get(key, .variant_job_probes); p$last_checked <- 0; assign(key, p, .variant_job_probes) }
      state <- read_variant_job(directory, root)
      if (state$status == expected) return(state)
      Sys.sleep(.025)
    }
    stop("Mock scheduler state not reached: ", expected)
  }
  writeLines("success", mode); done <- drive(first$job_dir, "COMPLETE")
  stopifnot(done$result_ready, done$scheduler_terminal, load_variant_job_result(first$job_dir, root)$synthetic,
    file.exists(file.path(first$job_dir, paste0("qacct-", first$job_id, ".txt"))))
  native_receipt <- job_read_json(file.path(first$job_dir, "application_receipt.json"))
  altered <- native_receipt; altered$attempt_token <- "obsolete"
  job_atomic_json(altered, file.path(first$job_dir, "application_receipt.json")); stopifnot(!job_result_valid(first$job_dir))
  job_atomic_json(native_receipt, file.path(first$job_dir, "application_receipt.json"))
  stopifnot(is.null(job_parse_accounting(c("jobnumber 999", "failed 0", "exit_status 0"), first$job_id)))
  # Resume can overlap other active variants, but needs terminal accounting for
  # this exact job; two clients resuming it must issue exactly one qsub call.
  state <- job_read_json(file.path(first$job_dir, "state.json"))
  state$status <- "FAILED"; state$result_ready <- FALSE; state$accounting <- NULL
  job_atomic_json(state, file.path(first$job_dir, "state.json"))
  writeLines("live", mode); expect_error(resume_variant_job(first$job_dir, root), "terminal scheduler")
  state$accounting <- list(job_id = first$job_id, failed = 0L, exit_status = 1L)
  job_atomic_json(state, file.path(first$job_dir, "state.json"))
  before <- as.integer(readLines(counter))
  attempts <- parallel::mclapply(1:2, function(i) tryCatch(resume_variant_job(first$job_dir, root), error = function(e) e),
    mc.cores = 2L, mc.preschedule = FALSE)
  ok <- which(!vapply(attempts, inherits, logical(1), "error"))
  stopifnot(length(ok) == 1L, as.integer(readLines(counter)) == before + 1L)
  resumed <- attempts[[ok]]
  stopifnot(resumed$job_id != first$job_id, resumed$attempt == 2L, resumed$cores == 16L,
    !identical(resumed$attempt_token, first$attempt_token), file.exists(file.path(first$job_dir, "submission-2.json")))
  # A stale worker cannot alter the resumed attempt or its result/progress.
  before_state <- readBin(file.path(first$job_dir, "state.json"), "raw", n = 1000000L)
  expect_error(run_variant_job(first$job_dir), "obsolete or different")
  stopifnot(identical(before_state, readBin(file.path(first$job_dir, "state.json"), "raw", n = 1000000L)),
    !file.exists(file.path(first$job_dir, "result.rds")))
  Sys.setenv(JOB_ID = resumed$job_id, REGSHINY_ATTEMPT_TOKEN = resumed$attempt_token)
  run_variant_job(first$job_dir)
  stopifnot(grepl(resumed$attempt_token, basename(execution$last$progress), fixed = TRUE), job_result_valid(first$job_dir))
  # Ambiguous submission never retries the same request and never blocks others.
  writeLines("unknown", submit_mode)
  unknown <- expect_error(submit(3L), "unrecognized job receipt")
  stopifnot(unknown$status == "SUBMISSION_UNKNOWN")
  before <- as.integer(readLines(counter))
  expect_error(job_submit_saved(unknown$job_dir, root), "unknown submission outcome")
  expect_error(resume_variant_job(unknown$job_dir, root), "terminal scheduler")
  stopifnot(as.integer(readLines(counter)) == before)
  writeLines("success", submit_mode); independent <- submit(2L)
  stopifnot(independent$status == "QUEUED")
  writeLines("reject", submit_mode)
  rejected <- expect_error(submit(1L), "synthetic scheduler rejection")
  stopifnot(rejected$status == "SUBMISSION_REJECTED")
  # An old immutable request retains its original resource contract on resume.
  writeLines("success", submit_mode)
  legacy_dir <- file.path(uidroot, "run-legacy"); dir.create(legacy_dir, mode = "0700")
  legacy <- request; legacy$cores <- 8L; legacy$workers <- 7L; legacy$worker_contract <- NULL
  job_write_request(legacy, legacy_dir); file.copy(file.path(first$job_dir, "run.sh"), file.path(legacy_dir, "run.sh"))
  job_atomic_json(list(status = "FAILED", result_ready = FALSE, job_id = "900", attempt = 1L,
    accounting = list(job_id = "900", failed = 0L, exit_status = 1L)), file.path(legacy_dir, "state.json"))
  old_request <- unname(tools::md5sum(file.path(legacy_dir, "request.rds")))
  old_resumed <- resume_variant_job(legacy_dir, root)
  stopifnot(old_resumed$cores == 8L, old_resumed$workers == 7L,
    identical(old_request, unname(tools::md5sum(file.path(legacy_dir, "request.rds")))))
  for (key in ls(.variant_job_probes)) {
    p <- get(key, .variant_job_probes)
    if (!is.null(p$process) && p$process$is_alive()) p$process$kill()
  }
})
cat("Parallel independent 16-core jobs, per-attempt isolation, accounting, logs and resume: PASS\n")

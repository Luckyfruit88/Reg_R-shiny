# Deterministic Shiny state tests. All records/results below are synthetic;
# no VCF, BAM, external executables or controlled data are used.
needed <- c("shiny", "bslib", "DT", "future", "promises", "later", "jsonlite")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("UI_TEST_DEPENDENCIES_MISSING: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(bslib))
source("annotation_module.R")
source("splice_evidence_module.R")
source("variant_module.R")
backend_file <- normalizePath("backend.R", mustWork = TRUE)
splice_file <- normalizePath("splice_evidence_backend.R", mustWork = TRUE)
backend <- new.env(parent = globalenv())
sys.source(backend_file, envir = backend)
fixture_file <- tempfile("regshiny-ui-fixture-", fileext = ".R")
writeLines(c(
  'lookup_variants <- function(resources, query) {',
  '  if (query == "chr1:999") stop("Synthetic lookup failure")',
  '  n <- if (query == "chr1:100") 2L else if (query == "chr1:200") 1L else 0L',
  '  data.frame(record_id = seq_len(n), chrom = rep("chr1", n), pos1 = rep(if (query == "chr1:200") 200L else 100L, n),',
  '    id = rep("synthetic", n), ref = rep("A", n), alt = c("G", "T")[seq_len(n)],',
  '    qual = rep("50", n), filter = rep("PASS", n), build = rep("SYNTHETIC", n), vcf = rep("synthetic.vcf.gz", n), contig_length = rep(220L, n))',
  '}',
  'analyze_variant <- function(resources, variant, cfg, max_per_group = 3L, sampling_mode = "all") {',
  '  list(variant = variant,',
  '    group_summary = data.frame(genotype = c("0/1", "NO_CALL"), call_status = c("CALLED", "NO_CALL"), dna_n = c(2L, 1L), available_n = c(2L, 0L),',
  '      selected_n = c(2L, 0L), analyzed_n = c(1L, 0L), failed_n = c(1L, 0L)),',
  '    genotypes = data.frame(vcf_sample = c("s1", "s2", "s3"), raw_gt = c("0|1", "1/0", "./."), genotype = c("0/1", "0/1", "NO_CALL"),',
  '      call_status = c("CALLED", "CALLED", "NO_CALL"), linked = c(TRUE, TRUE, FALSE)),',
  '    samples = data.frame(vcf_sample = c("s1", "s2"), genotype = "0/1", retained_reads = c(0, NA), status = c("SUCCESS", "FAILED"), error = c(NA, "Synthetic missing BAM")),',
  '    depth = data.frame(genotype = "0/1", chrom = cfg$chrom, pos1 = seq.int(cfg$start1, cfg$end1), mean_depth = 0, analyzed_n = 1L),',
  '    junctions = data.frame(genotype = character(), chrom = character(), intron_start1 = integer(), intron_end1 = integer(), strand = character(), count_sum = numeric(), count_mean = numeric(), analyzed_n = integer()),',
  '    rna_bases = data.frame(genotype = "0/1", base = c("A", "C", "G", "T", "N"), count = 0, analyzed_n = 1L),',
  '    provenance = list(source = "synthetic deterministic UI test", config = cfg, selection = list(mode = sampling_mode), rna_evidence_status = "SNV_BASE_COUNTS", completed_at = "synthetic"), warnings = "Synthetic selected failure retained")',
  '}'
), fixture_file)
future::plan(future::sequential)
html_text <- function(x) {
  value <- if (is.list(x) && !is.null(x$html)) as.character(x$html) else as.character(x)
  # HTML formatting inserts whitespace between text nodes; browsers collapse it.
  gsub("[[:space:]]+", " ", value)
}
await_task <- function(task, session) {
  for (i in seq_len(100L)) {
    later::run_now(0.05)
    session$flushReact()
    if (shiny::isolate(task$status()) != "running") {
      # ExtendedTask completion can occur during the preceding flush. Drain the
      # dependent output invalidations before asserting rendered UI text.
      later::run_now(0.01)
      session$flushReact()
      return(invisible(TRUE))
    }
  }
  stop("Synthetic async UI test timed out")
}
stopifnot(inherits(variant_ui("test"), "shiny.tag"))
shiny::testServer(variant_server,
  args = list(backend = backend, resources = list(build = "SYNTHETIC", vcf_registry = data.frame(chrom = "chr1"), samples = data.frame(vcf_sample = "s1")),
    backend_file = backend_file, variant_file = fixture_file, splice_file = splice_file), {
    session$setInputs(query = "chr1:100", record = "", sampling_mode = "preview", flank = 50, sample_cap = 3,
      mapq = 20, baseq = 0, anchor = 8, min_intron = 70, max_intron = 500000,
      strand_mode = "XS", exclude_dup = FALSE, nh1 = FALSE)
    session$setInputs(lookup = 1)
    await_task(lookup_task, session)
    stopifnot(isolate(lookup_task$status()) == "success", nrow(isolate(records())) == 2L)
    stopifnot(grepl("Choose the intended REF / ALT record explicitly", html_text(output$lookup_status), fixed = TRUE))
    # A position with two records cannot launch before an explicit allele choice.
    session$setInputs(run = 1)
    stopifnot(isolate(analysis_task$status()) == "initial")
    session$setInputs(record = "2", run = 2)
    await_task(analysis_task, session)
    stopifnot(isolate(analysis_task$status()) == "success")
    stopifnot(isolate(result())$variant$alt == "T", !isolate(results_stale()))
    stopifnot(grepl("Some selected samples failed", html_text(output$status), fixed = TRUE))
    stopifnot(isolate(result())$group_summary$failed_n[[1L]] == 1L)
    stopifnot(isolate(result())$genotypes$raw_gt[[3L]] == "./.")
    # Changed coordinate cannot silently reuse the old selected record.
    session$setInputs(query = "chr1:200", run = 3)
    stopifnot(isolate(results_stale()), isolate(result())$variant$pos1 == 100L)
    stopifnot(grepl("Inputs have changed", html_text(output$status), fixed = TRUE))
    stopifnot(grepl("Run Find variant", html_text(output$lookup_status), fixed = TRUE))
    # A failed fresh lookup invalidates selection; old results remain labelled stale.
    session$setInputs(query = "chr1:999", lookup = 2)
    await_task(lookup_task, session)
    stopifnot(isolate(lookup_task$status()) == "error", isolate(results_stale()))
    session$setInputs(record = "2", run = 4)
    stopifnot(isolate(result())$variant$pos1 == 100L)
    stopifnot(grepl("Synthetic lookup failure", html_text(output$lookup_status), fixed = TRUE))
    # No record is a separate state, never an inferred hom-reference call.
    session$setInputs(query = "chr1:300", lookup = 3)
    await_task(lookup_task, session)
    stopifnot(nrow(isolate(records())) == 0L)
    stopifnot(grepl("not evidence of a reference genotype", html_text(output$lookup_status), fixed = TRUE))
    # Repeating a prior coordinate still invalidates the prior snapshot revision.
    session$setInputs(query = "chr1:100", lookup = 4)
    await_task(lookup_task, session)
    session$setInputs(record = "2")
    stopifnot(isolate(results_stale()))
    session$setInputs(run = 5)
    await_task(analysis_task, session)
    stopifnot(!isolate(results_stale()), isolate(result())$variant$alt == "T")
    # A valid locus near the chromosome end is clipped to the declared contig.
    session$setInputs(query = "chr1:200", lookup = 5)
    await_task(lookup_task, session)
    session$setInputs(record = "1", flank = 50, run = 6)
    await_task(analysis_task, session)
    stopifnot(isolate(launched_config())$start1 == 150L, isolate(launched_config())$end1 == 220L)
    stopifnot(isolate(result())$variant$pos1 == 200L)
    # Unconfigured full mode never falls back to sampling, nor relabels old preview.
    session$setInputs(sampling_mode = "all", run = 7)
    stopifnot(isolate(active_mode()) == "preview", isolate(result_mode()) == "preview", isolate(results_stale()))
  })
shiny::testServer(variant_server,
  args = list(backend = backend, resources = simpleError("Synthetic unconfigured WGS resources"),
    backend_file = backend_file, variant_file = fixture_file, splice_file = splice_file), {
    stopifnot(grepl("Single BAM / demo tab remains available", html_text(output$resource_status), fixed = TRUE))
    session$setInputs(query = "chr1:100", lookup = 1, run = 1)
    stopifnot(isolate(lookup_task$status()) == "initial", isolate(analysis_task$status()) == "initial")
  })
# The job manager mock tests UI state only. It makes no scheduler submissions.
fixture <- new.env(parent = globalenv())
sys.source(fixture_file, envir = fixture)
job_fixture_file <- tempfile("regshiny-ui-jobs-", fileext = ".R")
writeLines(c(
  'submit_variant_job <- function(...) job_mock$submit_variant_job(...)',
  'read_variant_job <- function(...) job_mock$read_variant_job(...)',
  'resume_variant_job <- function(...) job_mock$resume_variant_job(...)'
), job_fixture_file)
job_mock <- new.env(parent = globalenv())
job_memory <- new.env(parent = emptyenv())
job_memory$submitted <- 0L; job_memory$resumed <- 0L; job_memory$loaded <- 0L
job_memory$state <- NULL; job_memory$result <- NULL
job_mock$list_variant_jobs <- function(job_root) {
  if (is.null(job_memory$state)) data.frame(job_dir = character(), label = character()) else
    data.frame(job_dir = "/synthetic/jobs/uid-test/job-1", label = "Synthetic saved full job")
}
job_mock$submit_variant_job <- function(resources, variant, cfg, backend_file, variant_file, job_source, job_root, ui_controls) {
  job_memory$submitted <- job_memory$submitted + 1L
  job_memory$state <- list(status = "QUEUED", job_id = "SYNTHETIC-1", result_ready = FALSE,
    progress = list(total = 2L, completed = 0L, success = 0L, failed = 0L, resumed = 0L,
      missing_link_n = 1L, incomplete_gt_n = 1L,
      group_summary = data.frame(genotype = c("0/1", "NO_CALL"), dna_n = c(2L, 1L),
        linked_n = c(2L, 0L), selected_n = c(2L, 0L), completed = 0L, success = 0L, failed = 0L)),
    request = list(variant = variant, cfg = cfg, ui_controls = ui_controls), error = NULL)
  job_memory$result <- fixture$analyze_variant(resources, variant, cfg, sampling_mode = "all")
  list(job_dir = "/synthetic/jobs/uid-test/job-1", job_id = "SYNTHETIC-1", status = "QUEUED")
}
job_mock$read_variant_job <- function(job_dir, job_root) job_memory$state
job_mock$load_variant_job_result <- function(job_dir, job_root) {
  stopifnot(isTRUE(job_memory$state$result_ready))
  job_memory$loaded <- job_memory$loaded + 1L
  job_memory$result
}
job_mock$resume_variant_job <- function(job_dir, job_root) {
  job_memory$resumed <- job_memory$resumed + 1L
  job_memory$state$status <- "QUEUED"
  list(job_dir = job_dir, job_id = "SYNTHETIC-2", status = "QUEUED")
}
full_args <- list(backend = backend,
  resources = list(build = "SYNTHETIC", vcf_registry = data.frame(chrom = "chr1"), samples = data.frame(vcf_sample = "s1")),
  backend_file = backend_file, variant_file = fixture_file, job_api = job_mock,
  job_source = job_fixture_file, job_root = "/synthetic/jobs", splice_file = splice_file)
shiny::testServer(variant_server, args = full_args, {
  # NULL scope (before the radio initializes) must also default to all.
  session$setInputs(query = "chr1:200", record = "", flank = 50, sample_cap = 1,
    mapq = 20, baseq = 0, anchor = 8, min_intron = 70, max_intron = 500000,
    strand_mode = "XS", exclude_dup = FALSE, nh1 = FALSE, lookup = 1)
  await_task(lookup_task, session)
  session$setInputs(record = "1", run = 1)
  await_task(single_submit_task, session)
  session$flushReact()
  stopifnot(isolate(active_mode()) == "all", job_memory$submitted == 1L)
  stopifnot(isolate(analysis_task$status()) == "initial", is.null(isolate(full_result())))
  stopifnot(is.null(job_memory$state$request$ui_controls$sample_cap))
  stopifnot(grepl("Selected: 2", html_text(output$job_progress), fixed = TRUE))
  stopifnot(grepl("Missing RNA link: 1", html_text(output$job_progress), fixed = TRUE))
  stopifnot(!is.null(output$progress_groups))
  # Artifact existence without the manager's validation gate cannot expose it.
  job_memory$state$status <- "ACCOUNTING_PENDING"
  session$elapse(2100); session$flushReact()
  stopifnot(job_memory$loaded == 0L, is.null(isolate(full_result())))
  session$setInputs(query = "chr1:100")
  stopifnot(grepl("Current inputs differ", html_text(output$job_progress), fixed = TRUE))
  job_memory$state$status <- "COMPLETE_WITH_FAILURES"
  job_memory$state$result_ready <- TRUE
  job_memory$state$progress <- list(total = 2L, completed = 2L, success = 1L, failed = 1L, resumed = 0L)
  session$elapse(2100); session$flushReact()
  stopifnot(job_memory$loaded == 1L, isolate(result_mode()) == "all")
  stopifnot(grepl("RNA evidence is partial", html_text(output$status), fixed = TRUE))
  stopifnot(grepl("ALL MATCHED BAMs", html_text(output$selected_variant), fixed = TRUE))
  stopifnot(isolate(result())$variant$pos1 == 200L)
  stopifnot(isolate(display_groups())$missing_rna_link_n[[2L]] == 1L)
})
# A new browser session can open the persisted result without launching work.
shiny::testServer(variant_server, args = full_args, {
  session$flushReact()
  session$setInputs(saved_job = "/synthetic/jobs/uid-test/job-1", open_job = 1)
  session$flushReact()
  stopifnot(job_memory$submitted == 1L, isolate(restored_job()))
  stopifnot(isolate(result())$variant$pos1 == 200L, isolate(launched_config())$end1 == 220L)
  stopifnot(grepl("saved job's submitted variant", html_text(output$status), fixed = TRUE))
  session$setInputs(open_job = 2)
  stopifnot(!is.null(isolate(full_result())), isolate(result())$variant$pos1 == 200L)
  # Interrupted jobs can be resumed from the explicit saved-job action.
  job_memory$state$result_ready <- FALSE; job_memory$state$status <- "INTERRUPTED"
  session$setInputs(resume_job = 1)
  await_task(resume_task, session)
  stopifnot(job_memory$resumed == 1L, job_memory$submitted == 1L, is.null(isolate(full_result())))
  # A fully successful result has a distinct, accurate completion message.
  job_memory$result$group_summary$analyzed_n[[1L]] <- 2L
  job_memory$result$group_summary$failed_n[[1L]] <- 0L
  job_memory$result$samples$status <- "SUCCESS"
  job_memory$result$warnings <- character()
  job_memory$state$status <- "COMPLETE"; job_memory$state$result_ready <- TRUE
  session$elapse(2100); session$flushReact()
  stopifnot(grepl("All matched BAMs completed successfully", html_text(output$status), fixed = TRUE))
  stopifnot(!grepl("RNA evidence is partial", html_text(output$status), fixed = TRUE))
})
# Dataset boundaries are independent of a variant coordinate or job label.
enabled <- reactiveVal(TRUE)
job_memory$reads <- 0L
job_mock$read_variant_job <- function(job_dir, job_root) {job_memory$reads <- job_memory$reads + 1L; job_memory$state}
job_memory$state$request$cfg$data_source_id <- "dataset-other"
source_args <- c(full_args, list(source_identity = list(id = "dataset-A", label = "Dataset A"), is_active = reactive(enabled())))
shiny::testServer(variant_server, args = source_args, {
  session$flushReact()
  before_resume <- job_memory$resumed
  session$setInputs(saved_job = "/synthetic/jobs/uid-test/job-1", open_job = 1, resume_job = 1)
  stopifnot(is.null(isolate(active_job())), is.null(isolate(full_result())), job_memory$resumed == before_resume)
  job_memory$state$request$cfg$data_source_id <- "dataset-A"
  job_memory$result$provenance$config$data_source_id <- "dataset-A"
  session$setInputs(open_job = 2)
  stopifnot(!is.null(isolate(full_result())))
  stopifnot(grepl("Dataset A", html_text(output$selected_variant), fixed = TRUE))
  enabled(FALSE); session$flushReact()
  stopifnot(is.null(isolate(full_result())), is.null(isolate(active_job())))
  reads_before <- job_memory$reads
  session$elapse(6000); session$flushReact()
  session$setInputs(open_job = 3, resume_job = 2, query = "chr1:100", lookup = 1, run = 1)
  stopifnot(job_memory$reads == reads_before, job_memory$resumed == before_resume)
  hidden <- tryCatch({isolate(result()); FALSE}, error = function(e) inherits(e, "shiny.silent.error"))
  stopifnot(hidden)
})
# Independent batch jobs and late asynchronous receipts use a separate fixture
# store. The job source is loaded by the same worker path used in production.
batch_variant_file <- tempfile("regshiny-batch-variants-", fileext = ".R")
writeLines(c(
  'lookup_variants <- function(resources, query) {',
  '  p <- strsplit(query, ":", fixed=TRUE)[[1L]]',
  '  if (!length(p) %in% c(2L,4L) || !grepl("^chr1$",p[[1L]]) || !grepl("^[0-9]+$",p[[2L]])) stop("Invalid variant coordinate")',
  '  position <- as.integer(p[[2L]])',
  '  alleles <- if(position==100L)c("G","T") else if(position %in% c(200L,210L))"G" else character()',
  '  x <- data.frame(chrom=rep("chr1",length(alleles)),pos1=rep(position,length(alleles)),ref=rep("A",length(alleles)),alt=alleles,',
  '    build=rep("SYNTHETIC",length(alleles)),id=rep("synthetic",length(alleles)),qual=rep("50",length(alleles)),filter=rep("PASS",length(alleles)),contig_length=rep(1000L,length(alleles)))',
  '  if(length(p)==4L)x<-x[x$ref==p[[3L]] & x$alt==p[[4L]],,drop=FALSE]',
  '  x$record_id<-paste(x$build,x$chrom,x$pos1,x$ref,x$alt,sep=":");x',
  '}'
), batch_variant_file)
batch_memory <- new.env(parent = emptyenv()); batch_memory$counter <- 0L; batch_memory$jobs <- list(); batch_memory$resumes <- 0L
batch_mock <- new.env(parent = globalenv())
batch_mock$submit_variant_job <- function(resources, variant, cfg, backend_file, variant_file, job_source, job_root, ui_controls) {
  batch_memory$counter <- batch_memory$counter + 1L
  path <- paste0("/synthetic/independent/job-", batch_memory$counter)
  failed <- variant$pos1[[1L]] == 210L
  state <- list(status = if (failed) "PREPARATION_FAILED" else "QUEUED", job_id = if (failed) NA_character_ else paste0("JOB-",batch_memory$counter),
    result_ready = FALSE, error = if(failed) "Synthetic matching failure" else NULL,
    progress = list(total=2L,completed=0L,success=0L,failed=0L,resumed=0L),
    request = list(variant=variant,cfg=cfg,ui_controls=ui_controls,cores=16L,workers=15L))
  batch_memory$jobs[[path]] <- list(state=state,result=fixture$analyze_variant(resources,variant,cfg,sampling_mode="all"))
  if(failed)stop(structure(list(message="Synthetic matching failure",call=NULL,job_dir=path,record_id=variant$record_id,
    status="PREPARATION_FAILED"),class=c("regshiny_job_error","error","condition")))
  list(job_dir=path,job_id=state$job_id,status="QUEUED",cores=16L,workers=15L)
}
batch_mock$list_variant_jobs <- function(job_root) {
  paths <- names(batch_memory$jobs)
  data.frame(job_dir=paths,label=paths,stringsAsFactors=FALSE)
}
batch_mock$list_variant_job_statuses <- function(job_root) {
  paths <- names(batch_memory$jobs)
  if(!length(paths))return(data.frame(job_dir=character(),label=character(),job_id=character(),status=character(),cores=integer()))
  do.call(rbind,lapply(paths,function(path){s<-batch_memory$jobs[[path]]$state
    data.frame(job_dir=path,label=path,job_id=s$job_id,status=s$status,variant=s$request$variant$record_id,
      cores=16L,workers=15L,total=2L,completed=s$progress$completed,success=s$progress$success,failed=s$progress$failed,
      result_ready=s$result_ready,stringsAsFactors=FALSE)}))
}
batch_mock$read_variant_job <- function(job_dir,job_root) batch_memory$jobs[[job_dir]]$state
batch_mock$load_variant_job_result <- function(job_dir,job_root) {
  stopifnot(isTRUE(batch_memory$jobs[[job_dir]]$state$result_ready));batch_memory$jobs[[job_dir]]$result
}
batch_mock$resume_variant_job <- function(job_dir,job_root) {
  batch_memory$resumes <- batch_memory$resumes+1L
  batch_memory$jobs[[job_dir]]$state$status <- "QUEUED"
  list(job_dir=job_dir,job_id="RESUMED",status="QUEUED")
}
job_mock <- batch_mock
batch_args <- full_args; batch_args$job_api <- batch_mock; batch_args$variant_file <- batch_variant_file
shiny::testServer(variant_server, args=batch_args, {
  session$setInputs(flank=50,mapq=20,baseq=0,anchor=8,min_intron=70,max_intron=500000,
    strand_mode="XS",exclude_dup=FALSE,nh1=FALSE,sampling_mode="all",batch_queries=paste(c(
      "chr1:200","chr1:100","chr1:100:A:T","chr1:200:A:G","bad coordinate","chr1:300","chr1:210"),collapse="\n"))
  session$setInputs(batch_submit=1)
  await_task(batch_task,session)
  out <- isolate(batch_outcomes())
  stopifnot(identical(out$outcome,c("SUBMITTED","AMBIGUOUS_ALLELES","SUBMITTED","DUPLICATE_IN_BATCH","LOOKUP_FAILED","NOT_FOUND","PREPARATION_FAILED")))
  stopifnot(batch_memory$counter==3L,length(unique(out$job_dir[out$outcome=="SUBMITTED"]))==2L)
  stopifnot(all(out$requested_cores[out$outcome=="SUBMITTED"]==16L))
  stopifnot(grepl("chr1:100:A:G",out$candidates[[2L]],fixed=TRUE),grepl("chr1:100:A:T",out$candidates[[2L]],fixed=TRUE))
  stopifnot(!is.na(out$job_dir[[7L]]),is.null(isolate(active_job())))
  dirs <- out$job_dir[out$outcome=="SUBMITTED"]
  configs <- lapply(dirs,function(path)batch_memory$jobs[[path]]$state$request$cfg)
  stopifnot(configs[[1L]]$start1==150L,configs[[2L]]$start1==50L)
  stopifnot(all(vapply(dirs,function(path)batch_memory$jobs[[path]]$state$request$cores,integer(1))==16L))
  stopifnot(grepl("2 independent jobs submitted",html_text(output$batch_status),fixed=TRUE))
})
# The receipt for A must not adopt B's controls or force the view back to A.
shiny::testServer(variant_server,args=batch_args,{
  session$setInputs(query="chr1:100:A:T",record="",flank=50,mapq=20,baseq=0,anchor=8,min_intron=70,max_intron=500000,
    strand_mode="XS",exclude_dup=FALSE,nh1=FALSE,sampling_mode="all",lookup=1)
  await_task(lookup_task,session)
  session$setInputs(record="SYNTHETIC:chr1:100:A:T",run=1)
  session$setInputs(query="chr1:200")
  await_task(single_submit_task,session)
  value<-isolate(single_receipt())
  stopifnot(isTRUE(value$ok),is.null(isolate(active_job())))
  request<-batch_memory$jobs[[value$receipt$job_dir]]$state$request
  stopifnot(request$variant$pos1==100L,request$cfg$start1==50L,request$ui_controls$query=="chr1:100:A:T")
})
# Explicitly opening an existing result while another submission is pending
# wins, even when the coordinate/filter controls have not changed.
old_path<-names(batch_memory$jobs)[[1L]]
batch_memory$jobs[[old_path]]$state$status<-"COMPLETE_WITH_FAILURES"
batch_memory$jobs[[old_path]]$state$result_ready<-TRUE
shiny::testServer(variant_server,args=batch_args,{
  session$setInputs(query="chr1:100:A:T",record="",flank=50,mapq=20,baseq=0,anchor=8,min_intron=70,max_intron=500000,
    strand_mode="XS",exclude_dup=FALSE,nh1=FALSE,sampling_mode="all",lookup=1)
  await_task(lookup_task,session)
  session$setInputs(record="SYNTHETIC:chr1:100:A:T",run=1)
  session$setInputs(saved_job=old_path,open_job=1)
  await_task(single_submit_task,session)
  stopifnot(isolate(active_job())==old_path,isolate(result())$variant$pos1==200L)
  stopifnot(isolate(single_receipt())$variant$pos1==100L)
})
unlink(c(fixture_file, job_fixture_file, batch_variant_file))
cat("VARIANT_UI_STATE_TESTS_OK\n")

# Synthetic per-session source connection tests. No study files are accessed.
needed <- c("shiny", "bslib", "future", "promises", "later", "jsonlite")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("DATA_SOURCE_UI_DEPENDENCIES_MISSING: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(shiny))
source("data_sources_module.R")
future::plan(future::sequential)
scratch <- tempfile("source-ui-"); dir.create(scratch)
source_file <- file.path(scratch, "sources_fixture.R")
writeLines(c(
  'default_fhs_spec <- function() list(kind="fhs",label="FHS fixture",build="GRCh38")',
  'custom_source_spec <- function(...) c(list(kind="custom"),list(...))',
  'prepare_data_profile <- function(spec,app_dir,state_root=NULL,progress_file=NULL,build_annotation=TRUE) {',
  '  if (!is.null(progress_file)) jsonlite::write_json(list(stage="VERIFYING",state="RUNNING",message=spec$label),progress_file,auto_unbox=TRUE)',
  '  if (spec$label == "FAIL") stop("Synthetic invalid source; exact mapping did not validate")',
  '  if (spec$label == "Slow old request") Sys.sleep(.25)',
  '  id <- paste0(spec$kind,"-",gsub("[^A-Za-z0-9]","_",spec$label))',
  '  root <- file.path(app_dir,"profiles");dir.create(root,showWarnings=FALSE)',
  '  path <- file.path(root,id);dir.create(path,showWarnings=FALSE)',
  '  bam <- if (!is.null(spec$bam) && nzchar(spec$bam)) spec$bam else paste0("/synthetic/",id,".bam")',
  '  bundle <- list(id=id,label=spec$label,variant_resources=list(build=spec$build,vcf_registry=data.frame(chrom="chr1"),',
  '    samples=data.frame(vcf_sample="synthetic",bam=bam)),bam_root=dirname(bam),bam_choices=setNames(bam,basename(bam)),',
  '    annotation_resources=simpleError("Synthetic optional annotation unavailable"),splice_reference=simpleError("Synthetic optional FASTA unavailable"),',
  '    job_root=file.path(path,"jobs"),profile_dir=path,config=list(spec=spec),issues=character())',
  '  saveRDS(bundle,file.path(path,"bundle.rds"));bundle',
  '}',
  'list_data_profiles <- function(app_dir,state_root=NULL) {',
  '  files<-list.files(file.path(app_dir,"profiles"),pattern="bundle.rds",recursive=TRUE,full.names=TRUE)',
  '  if (!length(files)) return(data.frame(profile_dir=character(),id=character(),label=character()))',
  '  do.call(rbind,lapply(files,function(f){b<-readRDS(f);data.frame(profile_dir=dirname(f),id=b$id,label=b$label)}))',
  '}',
  'load_data_profile <- function(profile_dir,app_dir,state_root=NULL) {',
  '  if(file.exists(file.path(profile_dir,"INVALID")))stop("Synthetic cached source changed")',
  '  readRDS(file.path(profile_dir,"bundle.rds"))',
  '}'
), source_file)
source_api <- new.env(parent = globalenv()); sys.source(source_file, envir = source_api)
stopifnot(inherits(data_sources_ui("source_test"), "shiny.tag"))
stopifnot(regshiny_ui_worker_count("1", "") == 1L, regshiny_ui_worker_count("2", "") == 1L,
  regshiny_ui_worker_count("3", "") == 2L, regshiny_ui_worker_count("32", "") == 2L,
  regshiny_ui_worker_count("unknown", "") == 1L, regshiny_ui_worker_count("3.5", "") == 1L,
  regshiny_ui_worker_count("3", "1") == 1L)
stopifnot(inherits(tryCatch(regshiny_ui_worker_count("1", "2"), error = function(e) e), "error"))
html_text <- function(x) gsub("[[:space:]]+", " ", if (is.list(x) && !is.null(x$html)) as.character(x$html) else as.character(x))
await_source <- function(task, session) {
  for (i in seq_len(160L)) {
    later::run_now(.05); session$flushReact()
    status <- isolate(task$status())
    if (status %in% c("success", "error")) {
      later::run_now(.01); session$flushReact(); return(invisible(status))
    }
  }
  stop("Synthetic preparation timed out")
}
original_environment <- Sys.getenv("REGSHINY_VCF_MANIFEST", unset = NA_character_)
Sys.setenv(REGSHINY_VCF_MANIFEST = "/synthetic/process-global-sentinel.tsv")
last_bundle <- NULL
shiny::testServer(data_sources_server,
  args = list(source_api = source_api, source_file = source_file, app_dir = scratch), {
    session$flushReact(); await_source(source_task, session)
    stopifnot(isolate(active_bundle())$label == "FHS fixture")
    fhs <- isolate(active_bundle())
    session$setInputs(kind="custom",label="Dataset B",build="GRCh38",vcf="/synthetic/cohort.vcf.gz",vcf_registry="",
      mapping="/synthetic/mapping.tsv",mapping_mode="explicit",bam_dir="",bam="/synthetic/B.bam",
      annotation_db="",annotation_gtf="",reference_fasta="",reference_receipt="",prepare=1)
    await_source(source_task, session)
    b <- isolate(active_bundle())
    stopifnot(b$label == "Dataset B", b$id != fhs$id, b$job_root != fhs$job_root)
    stopifnot(unname(b$bam_choices) == "/synthetic/B.bam", b$config$spec$mapping == "/synthetic/mapping.tsv")
    stopifnot(Sys.getenv("REGSHINY_VCF_MANIFEST") == "/synthetic/process-global-sentinel.tsv")
    # Invalid custom preparation must never replace a usable connection.
    session$setInputs(label="FAIL",prepare=2); await_source(source_task, session)
    stopifnot(isolate(source_task$status()) == "error", identical(isolate(active_bundle()), b))
    stopifnot(grepl("previous dataset remains connected", html_text(output$status), fixed=TRUE))
    # Cached profiles are explicitly revalidated before activation.
    file.create(file.path(fhs$profile_dir,"INVALID"))
    session$setInputs(saved=fhs$profile_dir,load=1); await_source(source_task, session)
    stopifnot(isolate(source_task$status()) == "error", identical(isolate(active_bundle()), b))
    unlink(file.path(fhs$profile_dir,"INVALID"))
    session$setInputs(load=2); await_source(source_task, session)
    stopifnot(isolate(active_bundle())$id == fhs$id)
    last_bundle <<- b
  })
# Another browser session starts with its own FHS choice; no global selection was mutated.
shiny::testServer(data_sources_server,
  args = list(source_api = source_api, source_file = source_file, app_dir = scratch), {
    session$flushReact(); await_source(source_task, session)
    stopifnot(isolate(active_bundle())$label == "FHS fixture", isolate(active_bundle())$id != last_bundle$id)
  })
# A single dedicated worker is genuinely asynchronous, and superseded work cannot connect old data.
future::plan(future::multisession, workers = I(1L))
stopifnot(future::value(future::future(Sys.getpid())) != Sys.getpid())
shiny::testServer(data_sources_server,
  args = list(source_api = source_api, source_file = source_file, app_dir = scratch), {
    session$flushReact(); await_source(source_task, session)
    session$setInputs(kind="custom",label="Slow old request",build="GRCh38",vcf="",vcf_registry="",mapping="",
      mapping_mode="explicit",bam_dir="",bam="/synthetic/old.bam",annotation_db="",annotation_gtf="",
      reference_fasta="",reference_receipt="",prepare=1)
    session$setInputs(label="Latest request",bam="/synthetic/new.bam",prepare=2)
    await_source(source_task, session)
    # Drain the queued invocation if its predecessor just completed.
    for (i in seq_len(80L)) {
      if (isolate(active_bundle())$label == "Latest request") break
      later::run_now(.05);session$flushReact()
    }
    stopifnot(isolate(active_bundle())$label == "Latest request", unname(isolate(active_bundle())$bam_choices) == "/synthetic/new.bam")
  })
future::plan(future::sequential)
# Exercise the actual app controller with lightweight viewer stubs. Every
# successful switch must instantiate new namespaces and deactivate old ones.
registered <- new.env(parent = emptyenv())
variant_server <- function(id, backend, resources, backend_file, variant_file, ..., source_identity = NULL, is_active = reactive(TRUE)) {
  force(resources); force(source_identity)
  registered[[id]] <- list(identity = source_identity, resources = resources, active = is_active)
}
single_bam_server <- function(id, bundle, is_active = reactive(TRUE)) {
  force(bundle)
  registered[[id]] <- list(identity = list(id = bundle$id), resources = bundle, active = is_active)
}
variant_ui <- function(id) shiny::div(id = id, "Synthetic WGS viewer")
single_bam_ui <- function(id, bam_choices) shiny::div(id = id, "Synthetic Single BAM viewer")
app_dir <- scratch; startup_spec <- source_api$default_fhs_spec()
backend <- new.env(); job_api <- new.env()
backend_file <- variant_file <- job_source <- annotation_file <- splice_file <- "/synthetic/source.R"
for (expression in parse("app.R")) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.symbol(expression[[2L]]) && identical(as.character(expression[[2L]]), "server")) eval(expression)
}
shiny::testServer(server, {
  session$flushReact(); await_source(sources$task, session)
  first <- isolate(current_view())
  stopifnot(first$generation == 1L, isolate(registered[[first$variant_id]]$active()))
  first_identity <- registered[[first$variant_id]]$identity$id
  session$setInputs(`sources-kind`="custom",`sources-label`="Controller B",`sources-build`="GRCh38",
    `sources-vcf`="/synthetic/B.vcf.gz",`sources-vcf_registry`="",`sources-mapping`="/synthetic/B.tsv",
    `sources-mapping_mode`="explicit",`sources-bam_dir`="",`sources-bam`="/synthetic/controller-B.bam",
    `sources-annotation_db`="",`sources-annotation_gtf`="",`sources-reference_fasta`="",`sources-reference_receipt`="",`sources-prepare`=1)
  await_source(sources$task, session)
  second <- isolate(current_view())
  stopifnot(second$generation == 2L, second$variant_id != first$variant_id, second$single_id != first$single_id)
  stopifnot(!isolate(registered[[first$variant_id]]$active()), !isolate(registered[[first$single_id]]$active()))
  stopifnot(isolate(registered[[second$variant_id]]$active()), registered[[first$variant_id]]$identity$id == first_identity)
  stopifnot(registered[[second$variant_id]]$identity$id != first_identity,
    unname(registered[[second$single_id]]$resources$bam_choices) == "/synthetic/controller-B.bam")
})
if (is.na(original_environment)) Sys.unsetenv("REGSHINY_VCF_MANIFEST") else Sys.setenv(REGSHINY_VCF_MANIFEST = original_environment)
unlink(scratch, recursive = TRUE)
cat("DATA_SOURCES_UI_TESTS_OK\n")

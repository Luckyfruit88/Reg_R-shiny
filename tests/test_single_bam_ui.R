# Native synthetic demo through the session-scoped Single BAM module.
# This checks that namespacing/resource isolation did not alter the old metrics.
source("runtime_setup.R")
regshiny_runtime_setup(strict_native = TRUE)
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(bslib))
source("annotation_module.R")
backend_file <- normalizePath("backend.R", mustWork = TRUE)
backend <- new.env(parent = globalenv()); sys.source(backend_file, envir = backend)
annotation_file <- normalizePath("annotation_backend.R", mustWork = TRUE)
annotation_backend <- new.env(parent = globalenv()); sys.source(annotation_file, envir = annotation_backend)
for (expression in parse("app.R")) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) && is.symbol(expression[[2L]]) &&
      as.character(expression[[2L]]) %in% c("single_bam_ui", "single_bam_server")) eval(expression)
}
future::plan(future::sequential)
unavailable <- simpleError("Synthetic profile has no DNA or GENCODE resource")
bundle <- list(id = "synthetic-demo-profile", label = "Synthetic demo profile", bam_choices = character(),
  variant_resources = unavailable, annotation_resources = unavailable)
enabled <- reactiveVal(TRUE)
stopifnot(grepl("demo_test-source", as.character(single_bam_ui("demo_test")), fixed = TRUE))
shiny::testServer(single_bam_server, args = list(bundle = bundle, is_active = reactive(enabled())), {
  session$setInputs(source = "demo", chrom = "chrDemo", roi_start = 101, roi_end = 450,
    mapq = 20, baseq = 0, anchor = 8, min_intron = 70, max_intron = 500000,
    strand_mode = "XS", exclude_dup = FALSE, nh1 = FALSE,
    target_start = 201, target_end = 300, delta = 5)
  session$flushReact()
  session$setInputs(run = 1)
  for (i in seq_len(160L)) {
    later::run_now(.05); session$flushReact()
    if (isolate(task$status()) %in% c("success", "error")) break
  }
  stopifnot(isolate(task$status()) == "success")
  m <- isolate(metrics())
  stopifnot(m$denominator_reads == 4L, m$junction_exact_count == 2L,
    m$junction_near_count == 1L, m$junction_per_100_reads == 50)
  stopifnot(isolate(result())$config$data_source_id == bundle$id)
  enabled(FALSE); session$flushReact()
  session$setInputs(run = 2)
  stopifnot(isolate(task$status()) == "success")
  hidden <- tryCatch({isolate(result()); FALSE}, error = function(e) inherits(e, "shiny.silent.error"))
  stopifnot(hidden)
})
cat("SINGLE_BAM_SESSION_UI_TESTS_OK\n")

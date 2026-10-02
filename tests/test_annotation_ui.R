# Synthetic annotation UI tests only: no GTF, BAM, VCF, or controlled inputs.
needed <- c("shiny", "bslib", "DT", "future", "promises", "later", "jsonlite")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("ANNOTATION_UI_DEPENDENCIES_MISSING: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(bslib))
source("annotation_module.R")
future::plan(future::sequential)
fixture_dir <- tempfile("regshiny-annotation-ui-"); dir.create(fixture_dir)
fixture_file <- file.path(fixture_dir, "annotation_fixture.R")
query_log <- file.path(fixture_dir, "queries.txt")
writeLines(character(), query_log)
writeLines(c(
  'query_annotation <- function(resources, chrom, start1, end1, build) {',
  '  cat(paste(chrom, start1, end1), "\\n", file = resources$query_log, append = TRUE)',
  '  if (chrom == "chrERROR") stop("Synthetic annotation query failure")',
  '  ids <- sprintf("ENST900000000%02d.1", seq_len(26L))',
  '  tx <- data.frame(transcript_id = ids, gene_id = "ENSG_SYNTHETIC.1", gene_name = "SYNTHETIC",',
  '    transcript_name = ids, transcript_type = "synthetic", chrom = chrom, start1 = 80L, end1 = 450L,',
  '    strand = rep(c("+", "-"), 13L), tags = "synthetic", transcript_support_level = "")',
  '  genes <- unique(tx[, c("gene_id", "gene_name", "chrom", "start1", "end1", "strand")])',
  '  exons <- do.call(rbind, lapply(seq_len(nrow(tx)), function(i) data.frame(transcript_id = ids[[i]],',
  '    gene_id = "ENSG_SYNTHETIC.1", chrom = chrom, start1 = c(80L, 250L), end1 = c(150L, 450L),',
  '    strand = tx$strand[[i]], exon_id = paste0(ids[[i]], c("-E1", "-E2")), exon_number = c(1L, 2L))))',
  '  features <- do.call(rbind, lapply(seq_len(nrow(tx)), function(i) data.frame(transcript_id = ids[[i]],',
  '    gene_id = "ENSG_SYNTHETIC.1", chrom = chrom, strand = tx$strand[[i]],',
  '    feature = c("UTR", "CDS", "CDS", "UTR"), start1 = c(80L, 120L, 250L, 361L), end1 = c(119L, 150L, 360L, 450L))))',
  '  introns <- data.frame(transcript_id = ids, chrom = chrom, strand = tx$strand, start1 = 151L, end1 = 249L)',
  '  if (chrom == "chr2") {genes <- genes[FALSE, ]; tx <- tx[FALSE, ]; exons <- exons[FALSE, ]; features <- features[FALSE, ]; introns <- introns[FALSE, ]}',
  '  list(genes = genes, transcripts = tx, exons = exons, features = features, introns = introns,',
  '    metadata = list(source = "SYNTHETIC_GENCODE_TEST_ONLY", release = "48", build = build, roi = list(chrom = chrom, start1 = start1, end1 = end1)))',
  '}',
  'annotate_junctions <- function(junctions, annotation) {',
  '  out <- junctions; out$annotation_status <- if (nrow(annotation$transcripts)) "EXACT_MATCH" else "NOT_IN_ANNOTATION"',
  '  list(junctions = out, matches = data.frame(transcript_id = annotation$transcripts$transcript_id[annotation$transcripts$strand == "+"]))',
  '}'
), fixture_file)
fixture <- new.env(parent = globalenv()); sys.source(fixture_file, envir = fixture)
resources <- list(build = "GRCh38", release = "48", query_log = query_log)
native <- data.frame(chrom = "chr1", intron_start1 = 151L, intron_end1 = 249L, strand = "+", score = 7L)
base_region <- list(chrom = "chr1", start1 = 100L, end1 = 400L, build = "GRCh38", demo = FALSE, source_id = "synthetic-bam-1")
model <- fixture$query_annotation(resources, "chr1", 100L, 400L, "GRCh38")
last_id <- tail(model$transcripts$transcript_id, 1L)
chosen <- annotation_pick_transcripts(model$transcripts, last_id)
stopifnot(nrow(chosen) == 20L, chosen$transcript_id[[1L]] == last_id, nrow(model$transcripts) == 26L)
stopifnot(nrow(annotation_pick_transcripts(model$transcripts, window = c(500L, 600L))) == 0L)
stopifnot(all(c("+", "-") %in% chosen$strand), all(grepl("\\.1$", chosen$transcript_id)))
ranked_models <- model$transcripts
ranked_models$tags[c(24L, 23L, 22L)] <- c("basic;MANE_Select", "Ensembl_canonical;basic", "basic")
ranked <- annotation_pick_transcripts(ranked_models, last_id)
stopifnot(identical(ranked$transcript_id[1:4], ranked_models$transcript_id[c(26L, 24L, 23L, 22L)]))
stopifnot(nrow(ranked_models) == 26L, !anyDuplicated(ranked$transcript_id))
png_file <- file.path(fixture_dir, "synthetic-model.png")
grDevices::png(png_file, width = 1200, height = 950)
plotted <- annotation_plot_models(model, c(100L, 400L), last_id, marker = 125L)
grDevices::dev.off()
stopifnot(file.info(png_file)$size > 1000L, identical(plotted$transcript_id, chosen$transcript_id))
stopifnot(inherits(annotation_track_ui("track"), "shiny.tag"))
stopifnot(inherits(annotation_detail_ui("detail", TRUE), "shiny.tag.list"))
html_text <- function(x) gsub("[[:space:]]+", " ", if (is.list(x) && !is.null(x$html)) as.character(x$html) else as.character(x))
await_task <- function(task, session) {
  for (i in seq_len(100L)) {
    later::run_now(.05); session$flushReact()
    if (shiny::isolate(task$status()) != "running") {
      later::run_now(.01); session$flushReact()
      return(invisible(TRUE))
    }
  }
  stop("Synthetic annotation task timed out")
}
region_value <- reactiveVal(base_region)
window_value <- reactiveVal(c(100L, 400L))
queries_before <- length(readLines(query_log))
shiny::testServer(annotation_server, args = list(resources = resources, annotation_file = fixture_file,
  region = reactive(region_value()), display_range = reactive(window_value()),
  junctions = reactive(native), marker = reactive(125L)), {
  session$flushReact(); await_task(annotation_task, session)
  stopifnot(isolate(annotation_task$status()) == "success", nrow(isolate(annotation())$transcripts) == 26L)
  stopifnot(nrow(isolate(selected_transcripts())) == 20L)
  stopifnot(identical(isolate(answer())$junctions$score, native$score), native$score == 7L)
  stopifnot(min(isolate(annotation())$exons$start1) == 80L, max(isolate(annotation())$exons$end1) == 450L)
  stopifnot(grepl("Drawing 20 of 26", html_text(output$display_note), fixed = TRUE))
  session$setInputs(priority = last_id)
  stopifnot(isolate(selected_transcripts())$transcript_id[[1L]] == last_id)
  window_value(c(300L, 380L)); session$flushReact()
  stopifnot(length(readLines(query_log)) == queries_before + 1L)
  stopifnot(grepl("Selected variant at 125 lies outside", html_text(output$display_note), fixed = TRUE))
  # A new result interval fetches reference data without changing RNA counts.
  r <- base_region; r$chrom <- "chr2"; region_value(r)
  session$flushReact(); await_task(annotation_task, session)
  stopifnot(nrow(isolate(annotation())$transcripts) == 0L, isolate(answer())$junctions$score == 7L)
  r$chrom <- "chrERROR"; region_value(r)
  session$flushReact(); await_task(annotation_task, session)
  stopifnot(isolate(annotation_task$status()) == "error")
  stopifnot(grepl("RNA metrics are unchanged", html_text(output$status), fixed = TRUE))
})
shiny::testServer(annotation_server, args = list(resources = simpleError("Synthetic missing annotation DB"),
  annotation_file = fixture_file, region = reactive(base_region), display_range = reactive(c(100L, 400L)),
  junctions = reactive(native)), {
  session$flushReact()
  stopifnot(isolate(annotation_task$status()) == "initial")
  stopifnot(grepl("RNA metrics remain available", html_text(output$status), fixed = TRUE))
})
demo <- base_region; demo$demo <- TRUE; demo$chrom <- "chrDemo"
shiny::testServer(annotation_server, args = list(resources = resources, annotation_file = fixture_file,
  region = reactive(demo), display_range = reactive(c(100L, 400L)), junctions = reactive(native)), {
  session$flushReact()
  stopifnot(isolate(annotation_task$status()) == "initial")
  stopifnot(grepl("synthetic chrDemo", html_text(output$status), fixed = TRUE))
})
wrong_build <- base_region; wrong_build$build <- "GRCh37"; wrong_build$reference_compatible <- TRUE
shiny::testServer(annotation_server, args = list(resources = resources, annotation_file = fixture_file,
  region = reactive(wrong_build), display_range = reactive(c(100L, 400L)), junctions = reactive(native),
  allow_confirmation = TRUE), {
  session$flushReact(); session$setInputs(confirm_grch38 = TRUE)
  stopifnot(isolate(annotation_task$status()) == "initial")
  stopifnot(grepl("differs from GENCODE v48 GRCh38", html_text(output$status), fixed = TRUE))
})
single <- base_region; single$build <- "unknown"; single$reference_compatible <- TRUE
single$reference_note <- "Synthetic matching reference length"
single_value <- reactiveVal(single)
shiny::testServer(annotation_server, args = list(resources = resources, annotation_file = fixture_file,
  region = reactive(single_value()), display_range = reactive(c(100L, 400L)), junctions = reactive(native),
  allow_confirmation = TRUE), {
  session$flushReact()
  stopifnot(isolate(annotation_task$status()) == "initial")
  session$setInputs(confirm_grch38 = TRUE)
  await_task(annotation_task, session)
  stopifnot(isolate(annotation_task$status()) == "success", isolate(effective_region())$build == "GRCh38")
  before <- length(readLines(query_log))
  x <- single; x$source_id <- "synthetic-bam-2"; single_value(x); session$flushReact()
  stopifnot(!is.null(isolate(effective_region())$unavailable), length(readLines(query_log)) == before)
  # Confirmation for one BAM never silently applies to a different BAM.
  session$setInputs(confirm_grch38 = FALSE); session$setInputs(confirm_grch38 = TRUE)
  await_task(annotation_task, session)
  stopifnot(isolate(effective_region())$build == "GRCh38", length(readLines(query_log)) == before + 1L)
  x$source_id <- "synthetic-bam-3"; x$reference_compatible <- FALSE; x$reference_note <- "Synthetic contig length mismatch"
  single_value(x); session$flushReact()
  session$setInputs(confirm_grch38 = FALSE); session$setInputs(confirm_grch38 = TRUE)
  stopifnot(grepl("length mismatch", html_text(output$status), fixed = TRUE))
  stopifnot(length(readLines(query_log)) == before + 1L)
})
unlink(fixture_dir, recursive = TRUE)
cat("ANNOTATION_UI_TESTS_OK\n")

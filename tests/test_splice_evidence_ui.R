# Synthetic-only UI and rendering tests. No BAM, VCF, FASTA or controlled data.
needed <- c("shiny", "bslib", "DT", "future", "promises", "later", "jsonlite")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("SPLICE_UI_DEPENDENCIES_MISSING: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(shiny))
source("splice_evidence_module.R")
future::plan(future::sequential)
scratch <- tempfile("splice-ui-"); dir.create(scratch)
fixture_source <- file.path(scratch, "model.R")
annotation_source <- file.path(scratch, "annotation.R")
query_log <- file.path(scratch, "annotation_queries.txt")
writeLines(character(), query_log)
fixture_build <- function(result, annotation, reference, flank = 25L) {
  v <- result$variant; cfg <- result$provenance$config; pos <- v$pos1[[1L]]
  w <- list(chrom = v$chrom[[1L]], start1 = max(cfg$start1, pos - flank),
    end1 = min(cfg$end1, pos + flank), variant_pos1 = pos, build = v$build[[1L]])
  g <- result$group_summary
  g$dna_alleles <- c("G/G", "G/A", "A/A")[match(g$genotype, c("0/0", "0/1", "1/1"))]
  g$motif_summary <- c("AG/AG", "AG/AA", "AA/AA")[match(g$genotype, c("0/0", "0/1", "1/1"))]
  status <- if (is.null(reference) || inherits(reference, "error")) "UNAVAILABLE_REFERENCE" else reference$status
  d <- result$depth[result$depth$pos1 >= w$start1 & result$depth$pos1 <= w$end1, ]
  d$coverage_status <- ifelse(is.finite(d$mean_depth), "OBSERVED", "UNAVAILABLE")
  positions <- seq.int(w$start1, w$end1)
  base <- rep("C", length(positions)); base[positions %in% c(pos - 1L, pos + 7L)] <- "A"
  base[positions %in% c(pos, pos + 8L)] <- "G"
  if (status != "REF_MATCH") base[] <- NA_character_
  bases <- data.frame(pos1 = positions, ref_base = base, alt_base = ifelse(positions == pos, "A", base), is_variant = positions == pos)
  sites <- data.frame(candidate_id = c("A0", "A8"), allele_index = 1L, site_type = "acceptor", strand = "+",
    start1 = c(pos - 1L, pos + 7L), end1 = c(pos, pos + 8L), anchor1 = c(pos, pos + 8L),
    ref_motif = "AG", alt_motif = c("AA", "AG"), canonical_expected = "AG",
    variant_overlaps = c(TRUE, FALSE), change = c("DISRUPTED_CANONICAL", "NOT_OVERLAPPING"),
    source = "RNA_OBSERVED;GENCODE", transcript_ids = "SYNTHETIC.1", distance_to_variant = c(0L, 8L))
  gs <- data.frame(genotype = c("0/0", "0/1", "1/1"), candidate_id = "A0", site_type = "acceptor", strand = "+",
    motif_alleles = c("AG/AG", "AG/AA", "AA/AA"), canonical_alleles_n = c(2L, 1L, 0L), ploidy = 2L,
    interpretation = "SYNTHETIC target-SNV motif only")
  j <- result$junctions
  j$junction_id <- paste(j$chrom, j$intron_start1, j$intron_end1, j$strand, sep = ":")
  j$distance_to_variant <- pmin(abs(j$intron_start1 - pos), abs(j$intron_end1 - pos))
  j$within_window <- j$distance_to_variant <= flank
  list(window = w, coverage = d, bases = bases, genotypes = g, sites = sites, genotype_sites = gs,
    junction_support = j, metadata = list(reference_status = status,
      focus_candidate_id = if (status == "REF_MATCH") "A0" else NA_character_, focus_allele_index = 1L,
      variant = as.list(v[1L, ]), source = "SYNTHETIC_UI_ONLY"),
    warnings = if (status != "REF_MATCH") "Synthetic reference unavailable / mismatched; motif inference held" else character())
}
dump("fixture_build", file = fixture_source)
cat("\nbuild_splice_evidence <- fixture_build\n", file = fixture_source, append = TRUE)
writeLines(c(
  'query_annotation <- function(resources, chrom, start1, end1, build) {',
  '  cat(paste(chrom, start1, end1), "\\n", file = resources$query_log, append = TRUE)',
  '  if (chrom == "chrDelay") Sys.sleep(.6)',
  '  list(exons = data.frame(chrom = chrom, start1 = 101L, end1 = 160L),',
  '    metadata = list(source = "SYNTHETIC", query = list(chrom = chrom, start1 = start1, end1 = end1)))',
  '}'
), annotation_source)
make_result <- function(chrom = "chr1") {
  groups <- data.frame(genotype = c("0/0", "0/1", "1/1"), call_status = "CALLED",
    dna_n = c(2L, 4L, 3L), selected_n = c(2L, 4L, 3L), analyzed_n = c(2L, 4L, 3L), failed_n = 0L)
  depth <- do.call(rbind, lapply(seq_len(3L), function(i) data.frame(genotype = groups$genotype[[i]], chrom = chrom,
    pos1 = 50:160, mean_depth = ifelse(50:160 > c(100L, 100L, 108L)[[i]], c(10, 7, 15)[[i]], 0), analyzed_n = groups$analyzed_n[[i]])))
  depth$mean_depth[depth$genotype == "1/1" & depth$pos1 == 90L] <- NA_real_
  j <- data.frame(genotype = rep(groups$genotype, 2L), chrom = chrom, intron_start1 = 50L,
    intron_end1 = rep(c(100L, 108L), each = 3L), strand = "+", count_sum = c(20, 12, 0, 0, 24, 45),
    count_mean = c(10, 3, 0, 0, 6, 15), analyzed_n = rep(groups$analyzed_n, 2L), source = "SAVED_RNA")
  list(variant = data.frame(chrom = chrom, pos1 = 100L, ref = "G", alt = "A", build = "GRCh38"),
    group_summary = groups, depth = depth, junctions = j,
    provenance = list(config = list(chrom = chrom, start1 = 50L, end1 = 160L), selection = list(mode = "all")))
}
reference <- list(status = "REF_MATCH")
annotation <- list(exons = data.frame(chrom = "chr1", start1 = 101L, end1 = 160L))
r <- make_result(); m <- fixture_build(r, annotation, reference)
stopifnot(identical(vapply(r$group_summary$genotype, function(g) splice_motif_label(m, g), character(1)),
  c("0/0" = "AG/AG", "0/1" = "AG/AA", "1/1" = "AA/AA")))
stopifnot(identical(splice_rna_boundaries(m)$pos1, c(100L, 108L)))
shown <- splice_junction_display(m)
stopifnot(length(unique(shown$junction_id)) == 2L, nrow(shown) == 6L)
stopifnot(shown$count_sum[shown$genotype == "1/1" & shown$intron_end1 == 100L] == 0)
stopifnot(shown$count_sum[shown$genotype == "0/0" & shown$intron_end1 == 108L] == 0)
draw_png <- file.path(scratch, "coverage.png")
grDevices::png(draw_png, width = 1800, height = 900, res = 160)
layout <- plot_splice_coverage(m, annotation)
grDevices::dev.off()
stopifnot(layout$ymax == 15, 108L %in% layout$boundary_positions, file.info(draw_png)$size > 1000)
grDevices::png(file.path(scratch, "junctions.png"), width = 1600, height = 500, res = 140)
matrix_layout <- plot_splice_junction_support(m)
grDevices::dev.off()
stopifnot(matrix_layout$max_mean == 15, length(matrix_layout$junction_ids) == 2L)
pdf_path <- file.path(scratch, "coverage.pdf")
render_warnings <- character()
withCallingHandlers({
  if (capabilities("cairo")) grDevices::cairo_pdf(pdf_path, width = 12, height = 6) else grDevices::pdf(pdf_path, width = 12, height = 6, useDingbats = FALSE)
  plot_splice_coverage(m, annotation); grDevices::dev.off()
}, warning = function(w) { render_warnings <<- c(render_warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
stopifnot(file.info(pdf_path)$size > 1000L, !length(render_warnings))
# Explicit unavailable focus cannot fall back to a fictitious motif change.
mismatch <- m; mismatch$metadata$reference_status <- "REF_MISMATCH"; mismatch$metadata$focus_candidate_id <- NA_character_
stopifnot(is.null(splice_primary_site(mismatch)), splice_motif_label(mismatch, "0/1") == "unavailable")
# Negative-strand motifs do not reverse-complement the DNA genotype labels.
negative <- m; negative$sites$strand <- "-"; negative$sites$alt_motif[[1L]] <- "TG"
negative$genotypes$dna_alleles <- c("T/T", "T/A", "A/A")
negative$genotype_sites$motif_alleles <- c("AG/AG", "AG/TG", "TG/TG")
stopifnot(splice_primary_site(negative)$strand == "-", splice_motif_label(negative, "0/1") == "AG/TG",
  negative$genotypes$dna_alleles[[2L]] == "T/A")
multi <- m; multi$metadata$variant$alt <- "A,C"; multi$sites$allele_index[[1L]] <- 2L; multi$metadata$focus_allele_index <- 2L
stopifnot(splice_alt_label(multi, splice_primary_site(multi)) == "ALT 2=C")
stopifnot(inherits(splice_evidence_ui("test"), "shiny.tag.list"))
html_text <- function(x) gsub("[[:space:]]+", " ", if (is.list(x) && !is.null(x$html)) as.character(x$html) else as.character(x))
await_model <- function(session, answer, model_task) {
  for (i in seq_len(160L)) {
    later::run_now(.05); session$flushReact()
    ready <- tryCatch({shiny::isolate(answer()); TRUE}, error = function(e) FALSE)
    if (ready && shiny::isolate(model_task$status()) == "success") {
      later::run_now(.01); session$flushReact(); return(invisible(TRUE))
    }
    if (shiny::isolate(model_task$status()) == "error") stop(shiny::isolate(model_task$result()))
  }
  stop("Synthetic splice model timed out")
}
current <- reactiveVal(r)
shiny::testServer(splice_evidence_server, args = list(result = reactive(current()), reference = reference,
  model_source = fixture_source, annotation_resources = list(query_log = query_log), annotation_source = annotation_source), {
  session$flushReact(); await_model(session, answer, model_task)
  stopifnot(isolate(model())$window$start1 == 75L, isolate(model())$window$end1 == 125L)
  stopifnot(identical(isolate(compact_result())$junctions, r$junctions))
  before <- length(readLines(query_log))
  session$setInputs(flank = "10"); await_model(session, answer, model_task)
  stopifnot(isolate(model())$window$start1 == 90L, length(readLines(query_log)) == before)
  stopifnot(grepl("reference status: REF_MATCH", html_text(output$status), fixed = TRUE))
  clipped <- r; clipped$provenance$config$start1 <- 99L; clipped$provenance$config$end1 <- 103L
  current(clipped)
  hidden <- tryCatch({isolate(answer()); FALSE}, error = function(e) inherits(e, "shiny.silent.error"))
  stopifnot(hidden)
  session$flushReact(); await_model(session, answer, model_task)
  stopifnot(isolate(model())$window$start1 == 99L, isolate(model())$window$end1 == 103L)
  stopifnot(length(unique(splice_near_junctions(isolate(model()))$junction_id)) == 1L)
})
shiny::testServer(splice_evidence_server, args = list(result = reactive(r), reference = NULL,
  model_source = fixture_source, annotation_resources = NULL, annotation_source = NULL), {
  session$flushReact(); await_model(session, answer, model_task)
  stopifnot(isolate(model())$metadata$reference_status == "UNAVAILABLE_REFERENCE")
  stopifnot(identical(isolate(model())$junction_support$count_sum, r$junctions$count_sum))
  stopifnot(any(is.finite(isolate(model())$coverage$mean_depth)), is.null(splice_primary_site(isolate(model()))))
})
# An actual delayed background annotation task must hide the previous ROI model.
future::plan(future::multisession, workers = 2L)
delayed <- reactiveVal(r)
shiny::testServer(splice_evidence_server, args = list(result = reactive(delayed()), reference = reference,
  model_source = fixture_source, annotation_resources = list(query_log = query_log), annotation_source = annotation_source), {
  session$flushReact(); await_model(session, answer, model_task)
  delayed(make_result("chrDelay")); session$flushReact()
  stopifnot(isolate(annotation_task$status()) == "running")
  hidden <- tryCatch({isolate(answer()); FALSE}, error = function(e) inherits(e, "shiny.silent.error"))
  stopifnot(hidden)
  await_model(session, answer, model_task)
  stopifnot(isolate(model())$window$chrom == "chrDelay")
})
future::plan(future::sequential)
unlink(scratch, recursive = TRUE)
cat("SPLICE_EVIDENCE_UI_TESTS_OK\n")

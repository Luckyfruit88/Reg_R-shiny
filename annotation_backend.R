# GENCODE v48 overlay for Reg_Shiny. This module does not rerun or alter RNA
# counting, genotype grouping, persistent variant jobs or their fingerprints.
.regshiny_annotation_source <- local({
  paths <- vapply(sys.frames(), function(frame) {
    path <- frame$ofile
    if (is.character(path) && length(path) == 1L) path else ""
  }, character(1))
  paths <- paths[nzchar(paths)]
  normalizePath(if (length(paths)) tail(paths, 1L) else "annotation_backend.R", mustWork = FALSE)
})

annotation_integer <- function(x, name, lower = 1L, upper = 2e9) {
  if (length(x) != 1L || is.na(x) || !is.numeric(x) || !is.finite(x) ||
      x != floor(x) || x < lower || x > upper)
    stop(name, " must be an integer from ", lower, " through ", upper, ".")
  as.integer(x)
}

annotation_file_state <- function(path) {
  s <- file.info(path)
  list(path = path, size_bytes = unname(s$size), mtime = as.numeric(s$mtime))
}

annotation_command <- function(args, script = file.path(dirname(.regshiny_annotation_source), "scripts", "build_annotation.py")) {
  python <- Sys.getenv("REGSHINY_PYTHON", "")
  if (!nzchar(python)) python <- unname(Sys.which("python3"))
  if (!nzchar(python) || !file.exists(python)) stop("Python 3 is required for read-only GENCODE annotation queries.")
  if (!file.exists(script)) stop("The GENCODE annotation query script is missing.")
  output <- tempfile("regshiny_annotation_", fileext = ".json")
  file.create(output); Sys.chmod(output, "0600")
  on.exit(unlink(output), add = TRUE)
  p <- processx::run(python, c(script, as.character(args)), stdout = output, stderr = "|",
                     timeout = 120, error_on_status = FALSE, cleanup_tree = FALSE)
  if (isTRUE(p$timeout) || is.na(p$status) || p$status != 0L)
    stop(trimws(paste("GENCODE query failed:", p$stderr)), call. = FALSE)
  if (file.info(output)$size > 128 * 1024^2)
    stop("The complete GENCODE response exceeds 128 MiB. Select a smaller interval; no annotation was truncated.")
  jsonlite::fromJSON(output, simplifyVector = TRUE)
}

load_annotation_resources <- function(db = Sys.getenv("REGSHINY_ANNOTATION_DB")) {
  if (length(db) != 1L || is.na(db) || !nzchar(db) || !file.exists(db) || file.access(db, 4) != 0L)
    stop("GENCODE v48 is not configured or its annotation index is not readable.")
  db <- normalizePath(db, mustWork = TRUE)
  state <- annotation_file_state(db)
  metadata <- annotation_command(c("metadata", "--db", db))
  if (!identical(state, annotation_file_state(db))) stop("The annotation index changed during initialization.")
  if (!identical(metadata$release, "48") || !identical(metadata$build, "GRCh38") ||
      !identical(metadata$schema_version, "regshiny_gencode_v1"))
    stop("The annotation resource must be GENCODE v48 on GRCh38.")
  structure(list(db = db, build = metadata$build, release = metadata$release,
    metadata = metadata, source_state = state,
    query_script = file.path(dirname(.regshiny_annotation_source), "scripts", "build_annotation.py")),
    class = "regshiny_annotation_resources")
}

annotation_reference_length <- function(resources, chrom) {
  if (!inherits(resources, "regshiny_annotation_resources")) stop("Load the GENCODE annotation resource first.")
  lengths <- resources$metadata$reference_contig_lengths
  if (is.null(lengths) || length(chrom) != 1L || is.na(chrom) || !chrom %in% names(lengths)) return(NA_real_)
  value <- suppressWarnings(as.numeric(lengths[[chrom]]))
  if (length(value) != 1L || !is.finite(value) || value < 1) NA_real_ else value
}

annotation_table <- function(value, schema) {
  if (is.null(value) || !length(value)) {
    columns <- lapply(schema, function(type) if (type == "integer") integer() else character())
    return(as.data.frame(columns, stringsAsFactors = FALSE))
  }
  if (!is.data.frame(value) || !all(names(schema) %in% names(value))) stop("Unexpected annotation table schema.")
  value[, names(schema), drop = FALSE]
}

query_annotation <- function(resources, chrom, start1, end1, build = "GRCh38",
                             max_transcripts = 5000L, max_records = 200000L) {
  if (!inherits(resources, "regshiny_annotation_resources")) stop("Load the GENCODE annotation resource first.")
  if (length(build) != 1L || is.na(build) || !identical(build, resources$build))
    stop("Variant/BAM and GENCODE reference assemblies differ; no liftover is performed.")
  if (length(chrom) != 1L || is.na(chrom) || !grepl("^[A-Za-z0-9_][A-Za-z0-9_.-]*$", chrom))
    stop("Use an exact, literal annotation chromosome name.")
  start1 <- annotation_integer(start1, "Annotation start")
  end1 <- annotation_integer(end1, "Annotation end", start1)
  if (end1 - start1 + 1L > 250000L) stop("The maximum annotation interval is 250,000 bp.")
  max_transcripts <- annotation_integer(max_transcripts, "Transcript limit", 1L, 50000L)
  max_records <- annotation_integer(max_records, "Annotation record limit", 1L, 2000000L)
  if (!identical(resources$source_state, annotation_file_state(resources$db)))
    stop("The GENCODE index changed; reload its resource before querying.")
  answer <- annotation_command(c("query", "--db", resources$db, "--chrom", chrom,
    "--start", start1, "--end", end1, "--build", build,
    "--max-transcripts", max_transcripts, "--max-records", max_records), resources$query_script)
  if (!identical(resources$source_state, annotation_file_state(resources$db)))
    stop("The GENCODE index changed during query; results were not published.")
  answer$genes <- annotation_table(answer$genes, c(gene_id = "character", gene_name = "character", gene_type = "character",
    chrom = "character", start1 = "integer", end1 = "integer", strand = "character", source = "character"))
  answer$transcripts <- annotation_table(answer$transcripts, c(transcript_id = "character", gene_id = "character",
    gene_name = "character", transcript_name = "character", transcript_type = "character", transcript_support_level = "character",
    tags = "character", chrom = "character", start1 = "integer", end1 = "integer", strand = "character", source = "character"))
  answer$exons <- annotation_table(answer$exons, c(transcript_id = "character", gene_id = "character", exon_id = "character",
    exon_number = "character", chrom = "character", start1 = "integer", end1 = "integer", strand = "character"))
  answer$features <- annotation_table(answer$features, c(transcript_id = "character", gene_id = "character", feature = "character",
    chrom = "character", start1 = "integer", end1 = "integer", strand = "character", frame = "character"))
  answer$introns <- annotation_table(answer$introns, c(transcript_id = "character", gene_id = "character", chrom = "character",
    start1 = "integer", end1 = "integer", strand = "character", intron_number = "integer",
    intron_start1 = "integer", intron_end1 = "integer", donor1 = "integer", acceptor1 = "integer"))
  answer$metadata$index_path <- resources$db
  answer$metadata$index_state <- resources$source_state
  answer$metadata$queried_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  answer
}

annotate_junctions <- function(junctions, annotation) {
  required <- c("chrom", "intron_start1", "intron_end1", "strand")
  if (!is.data.frame(junctions) || !all(required %in% names(junctions)))
    stop("Junction annotation requires chrom, intron_start1, intron_end1 and strand.")
  q <- annotation$metadata$query
  if (is.null(q) || !isTRUE(q$full_transcript_models) || !identical(q$truncated, FALSE))
    stop("Junction annotation requires a complete, untruncated transcript-model query.")
  out <- junctions
  for (name in c("annotation_status", "annotation_label", "matching_transcripts", "matching_genes", "matching_gene_names"))
    out[[name]] <- rep("", nrow(out))
  out$annotation_match_n <- integer(nrow(out))
  out$annotation_release <- rep("GENCODE v48", nrow(out))
  details <- list()
  introns <- annotation$introns
  key <- function(chrom, start, end) paste(chrom, start, end, sep = "|")
  intron_key <- key(introns$chrom, introns$intron_start1, introns$intron_end1)
  for (i in seq_len(nrow(out))) {
    start <- out$intron_start1[[i]]; end <- out$intron_end1[[i]]; strand <- as.character(out$strand[[i]])
    if (!is.finite(start) || !is.finite(end) || start != floor(start) || end != floor(end) || start < 1 || end < start)
      stop("Invalid 1-based inclusive junction boundaries.")
    if (out$chrom[[i]] != q$chrom || start < q$start1 || end > q$end1) {
      out$annotation_status[[i]] <- "OUTSIDE_QUERY"
      out$annotation_label[[i]] <- "Outside the complete annotation query interval"
      next
    }
    exact <- which(intron_key == key(out$chrom[[i]], start, end))
    known_strand <- !is.na(strand) && strand %in% c("+", "-")
    if (known_strand) exact <- exact[introns$strand[exact] == strand]
    if (!known_strand) {
      out$annotation_status[[i]] <- "UNKNOWN_STRAND"
      out$annotation_label[[i]] <- if (length(exact)) "Exact annotation coordinates; RNA strand unresolved" else
        "No exact coordinate match in this annotation; RNA strand unresolved"
    } else if (length(exact)) {
      out$annotation_status[[i]] <- "ANNOTATED"
      out$annotation_label[[i]] <- "Exact GENCODE v48 junction boundaries and strand"
    } else {
      out$annotation_status[[i]] <- "NOT_IN_ANNOTATION"
      out$annotation_label[[i]] <- "Not in GENCODE v48 at these exact boundaries and strand"
    }
    if (length(exact)) {
      matches <- introns[exact, , drop = FALSE]
      matches$gene_name <- annotation$genes$gene_name[match(matches$gene_id, annotation$genes$gene_id)]
      matches$observed_row <- i
      matches$observed_strand <- strand
      matches$match_type <- if (known_strand) "EXACT_BOUNDARIES_AND_STRAND" else "EXACT_BOUNDARIES_STRAND_UNRESOLVED"
      details[[length(details) + 1L]] <- matches
      join <- function(x) paste(sort(unique(x[!is.na(x) & nzchar(x)])), collapse = ";")
      out$matching_transcripts[[i]] <- join(matches$transcript_id)
      out$matching_genes[[i]] <- join(matches$gene_id)
      out$matching_gene_names[[i]] <- join(matches$gene_name)
      out$annotation_match_n[[i]] <- nrow(matches)
    }
  }
  empty <- introns[FALSE, , drop = FALSE]
  empty$gene_name <- character(); empty$observed_row <- integer(); empty$observed_strand <- character(); empty$match_type <- character()
  list(junctions = out, matches = if (length(details)) do.call(rbind, details) else empty,
       metadata = list(release = annotation$metadata$release, build = annotation$metadata$build,
                       source_sha256 = annotation$metadata$source_sha256,
                       definition = "Exact two intron boundaries plus strand; unresolved RNA strand is never a definitive unannotated call."))
}

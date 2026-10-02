# Rscript tests/test_annotations.R
# Public synthetic GTF fixture: no controlled data and no BAM analysis.
source("annotation_backend.R", encoding = "UTF-8")

annotation_expect_error <- function(expr, text = NULL) {
  e <- tryCatch({ force(expr); NULL }, error = function(e) e)
  stopifnot(inherits(e, "error"))
  if (!is.null(text)) stopifnot(grepl(text, conditionMessage(e), fixed = TRUE))
}

run_annotation_tests <- function() {
  python <- unname(Sys.which("python3"))
  stopifnot(nzchar(python))
  work <- tempfile("regshiny_annotation_test_"); dir.create(work, mode = "0700")
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  script <- normalizePath("scripts/build_annotation.py")
  gtf <- file.path(work, "gencode.v48.synthetic.gtf.gz")
  db <- file.path(work, "annotation.sqlite")
  fai <- file.path(work, "synthetic.fa.fai")
  writeLines(c("chr1\t1000\t0\t50\t51", "chrY\t1000\t0\t50\t51"), fai)
  rows <- character()
  row <- function(feature, start, end, strand, attributes, chrom = "chr1", frame = ".")
    paste(c(chrom, "SYNTHETIC", feature, start, end, ".", strand, frame, attributes), collapse = "\t")
  gene <- function(id, name, start, end, strand, chrom = "chr1")
    row("gene", start, end, strand, paste0('gene_id "', id, '"; gene_name "', name,
       '"; gene_type "protein_coding";'), chrom)
  model <- function(id, gene_id, name, starts, ends, strand = "+", chrom = "chr1", tags = character()) {
    common <- paste0('gene_id "', gene_id, '"; transcript_id "', id, '"; gene_name "', name,
      '"; transcript_name "', id, '"; transcript_type "protein_coding";')
    tx <- row("transcript", min(starts), max(ends), strand,
      paste(common, paste0('tag "', tags, '";', collapse = " ")), chrom)
    numbers <- if (strand == "+") seq_along(starts) else rev(seq_along(starts))
    exons <- vapply(seq_along(starts), function(i) row("exon", starts[[i]], ends[[i]], strand,
      paste0(common, ' exon_id "EX', id, '.', i, '"; exon_number "', numbers[[i]], '";'), chrom), character(1))
    c(tx, exons)
  }
  rows <- c("##description: synthetic test annotation of human genome (GRCh38), version 48 (Ensembl 114)",
    "##provider: SYNTHETIC FIXTURE; not biological evidence",
    gene("G1.1", "Gene;One", 100, 500, "+"),
    model("T1.2", "G1.1", "Gene;One", c(100, 200, 400), c(149, 249, 500), tags = c("basic", "Ensembl_canonical", "MANE_Select")),
    model("T2.4", "G1.1", "Gene;One", c(100, 300, 400), c(149, 349, 500)),
    model("T3.1", "G1.1", "Gene;One", c(100, 200, 400), c(149, 249, 500)),
    gene("G2.3", "Reverse", 600, 900, "-"),
    model("T4.5", "G2.3", "Reverse", c(600, 700, 800), c(649, 749, 900), "-"),
    gene("G1.1_PAR_Y", "PARgene", 100, 150, "+", "chrY"),
    model("T1.2_PAR_Y", "G1.1_PAR_Y", "PARgene", 100, 150, "+", "chrY"),
    gene("G_ALT.1", "ALT", 1, 10, "+", "chr1_alt"),
    model("T_ALT.2", "G_ALT.1", "ALT", 1, 10, "+", "chr1_alt"))
  attrs <- 'gene_id "G1.1"; transcript_id "T1.2";'
  rows <- c(rows, row("CDS", 110, 149, "+", attrs, frame = "0"),
    row("CDS", 200, 249, "+", attrs, frame = "2"), row("CDS", 400, 480, "+", attrs, frame = "0"),
    row("UTR", 100, 109, "+", attrs), row("UTR", 481, 500, "+", attrs))
  write_gtf <- function(lines, path) {
    con <- gzfile(path, "wt"); on.exit(close(con)); writeLines(lines, con)
  }
  build <- function(source, output, extra = character(), should_fail = FALSE) {
    p <- processx::run(python, c(script, "build", "--gtf", source, "--output", output,
      "--source-url", "https://example.invalid/synthetic-gencode48.gtf.gz", extra),
      timeout = 60, error_on_status = FALSE)
    if (!should_fail && p$status != 0L) stop(p$stderr)
    p
  }
  write_gtf(rows, gtf)
  built <- build(gtf, db, c("--fai", fai))
  receipt <- jsonlite::fromJSON(built$stdout)
  stopifnot(file.exists(db), receipt$metadata$counts$genes == 4L,
            receipt$metadata$counts$transcripts == 6L, receipt$metadata$counts$exons == 14L,
            receipt$metadata$counts$introns == 8L, nchar(receipt$metadata$source_sha256) == 64L,
            nchar(receipt$metadata$builder_sha256) == 64L, nchar(receipt$database_sha256) == 64L)
  resources <- load_annotation_resources(db)
  stopifnot(resources$release == "48", resources$build == "GRCh38",
    annotation_reference_length(resources, "chr1") == 1000,
    is.na(annotation_reference_length(resources, "chr1_alt")),
    is.na(annotation_reference_length(resources, "1")))
  # A wholly intronic window still returns complete transcript/exon models.
  intronic <- query_annotation(resources, "chr1", 160L, 180L)
  stopifnot(nrow(intronic$transcripts) == 3L, nrow(intronic$exons) == 9L,
            min(intronic$exons$start1) == 100L, max(intronic$exons$end1) == 500L,
            any(intronic$introns$start1 == 150L & intronic$introns$end1 == 199L),
            !any(intronic$introns$start1 == 160L | intronic$introns$end1 == 180L),
            !intronic$metadata$query$truncated,
            intronic$transcripts$tags[intronic$transcripts$transcript_id == "T1.2"] == "basic;Ensembl_canonical;MANE_Select")
  annotation <- query_annotation(resources, "chr1", 100L, 900L)
  stopifnot(nrow(annotation$features) == 5L,
            all(c("CDS", "UTR") %in% annotation$features$feature),
            annotation$genes$gene_name[annotation$genes$gene_id == "G1.1"] == "Gene;One")
  negative <- annotation$introns[annotation$introns$transcript_id == "T4.5", ]
  stopifnot(identical(as.integer(negative$intron_number), c(2L, 1L)),
            all(negative$donor1 == negative$end1), all(negative$acceptor1 == negative$start1))
  j <- data.frame(chrom = "chr1", intron_start1 = c(150L, 650L, 150L, 650L, 150L, 550L, 999L),
    intron_end1 = c(199L, 699L, 198L, 699L, 199L, 560L, 1001L), strand = c("+", "-", "+", "+", "?", "?", "+"),
    count_sum = c(1, 2, 3, 4, 5, 6, 7))
  matches <- annotate_junctions(j, annotation)
  stopifnot(identical(matches$junctions$count_sum, j$count_sum),
    identical(matches$junctions$annotation_status,
      c("ANNOTATED", "ANNOTATED", "NOT_IN_ANNOTATION", "NOT_IN_ANNOTATION", "UNKNOWN_STRAND", "UNKNOWN_STRAND", "OUTSIDE_QUERY")),
    matches$junctions$annotation_match_n[[1]] == 2L,
    matches$junctions$matching_transcripts[[1]] == "T1.2;T3.1",
    sum(matches$matches$observed_row == 1L) == 2L,
    all(matches$matches$match_type[matches$matches$observed_row == 5L] == "EXACT_BOUNDARIES_STRAND_UNRESOLVED"),
    !any(grepl("novel", matches$junctions$annotation_label, ignore.case = TRUE)))
  par <- query_annotation(resources, "chrY", 100L, 100L)
  stopifnot(par$transcripts$transcript_id == "T1.2_PAR_Y", nrow(par$introns) == 0L)
  empty <- query_annotation(resources, "chr1", 950L, 999L)
  stopifnot(nrow(empty$transcripts) == 0L, nrow(empty$exons) == 0L,
    all(c("transcript_id", "start1", "end1") %in% names(empty$exons)))
  annotation_expect_error(query_annotation(resources, "1", 100L, 200L), "not present")
  annotation_expect_error(query_annotation(resources, "chr1", 100L, 200L, build = "GRCh37"), "assemblies differ")
  annotation_expect_error(query_annotation(resources, "chr1; rm", 100L, 200L), "literal")
  annotation_expect_error(query_annotation(resources, "chr1", 0L, 200L), "integer")
  annotation_expect_error(query_annotation(resources, "chr1", 1L, 250001L), "250,000")
  annotation_expect_error(query_annotation(resources, "chr1", 100L, 900L, max_transcripts = 2L), "transcript limit")
  annotation_expect_error(query_annotation(resources, "chr1", 100L, 900L, max_records = 5L), "record limit")
  invalid <- annotation; invalid$metadata$query$truncated <- TRUE
  annotation_expect_error(annotate_junctions(j, invalid), "untruncated")
  # Read-only calls leave the derivative byte-identical and create no journals.
  hash_before <- unname(tools::md5sum(db)); Sys.chmod(db, "0440")
  query_annotation(resources, "chr1", 100L, 900L)
  stopifnot(identical(hash_before, unname(tools::md5sum(db))),
            !any(file.exists(paste0(db, c("-journal", "-wal", "-shm")))))
  old_mtime <- file.info(db)$mtime; Sys.setFileTime(db, old_mtime + 2)
  annotation_expect_error(query_annotation(resources, "chr1", 100L, 900L), "index changed")
  Sys.setFileTime(db, old_mtime)
  # Missing internal/terminal exons must not create artificial annotated introns.
  bad <- file.path(work, "bad.gtf.gz"); bad_db <- file.path(work, "bad.sqlite")
  internal <- grepl('\texon\t200\t249\t', rows) & grepl('transcript_id "T1.2";', rows, fixed = TRUE)
  write_gtf(rows[!internal], bad)
  p <- build(bad, bad_db, should_fail = TRUE)
  stopifnot(p$status != 0L, grepl("exon_number sequence", p$stderr), !file.exists(bad_db))
  terminal <- grepl('\texon\t400\t500\t', rows) & grepl('transcript_id "T1.2";', rows, fixed = TRUE)
  # Drop this transcript's now-outside CDS/UTR to reach the exon completeness check.
  features <- grepl('\t(CDS|UTR)\t', rows) & grepl('transcript_id "T1.2";', rows, fixed = TRUE)
  write_gtf(rows[!terminal & !features], bad)
  p <- build(bad, bad_db, should_fail = TRUE)
  stopifnot(p$status != 0L, grepl("full transcript boundaries", p$stderr), !file.exists(bad_db))
  write_gtf(sub("version 48", "version 47", rows, fixed = TRUE), bad)
  p <- build(bad, bad_db, should_fail = TRUE)
  stopifnot(p$status != 0L, grepl("version 48", p$stderr), !file.exists(bad_db))
  writeLines("chr1\t499\t0\t50\t51", fai)
  p <- build(gtf, bad_db, c("--fai", fai), should_fail = TRUE)
  stopifnot(p$status != 0L, grepl("FAI contig length", p$stderr), !file.exists(bad_db))
  p <- build(gtf, db, should_fail = TRUE)
  stopifnot(p$status != 0L, grepl("already exists", p$stderr))
  cat("GENCODE index/query/strand/exon-completeness/one-to-many/read-only/provenance tests: PASS\n")
  invisible(TRUE)
}

run_annotation_tests()

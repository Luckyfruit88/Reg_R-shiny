# Rscript tests/test_backend.R [--unit]
# backend.R: CIGAR/coordinate parsing, BAM contig choice and RNA base counting.
# --unit skips the native samtools / RegTools integration section.
unit_only <- "--unit" %in% commandArgs(trailingOnly = TRUE)
source("backend.R", encoding = "UTF-8")
cat("Core CIGAR / coordinate tests\n")
a <- parse_cigar("20M100N30M", 181L)
stopifnot(identical(as.integer(a$introns[1, ]), c(200L, 300L)))
stopifnot(identical(as.integer(a$blocks[1, ]), c(180L, 200L)))
stopifnot(nrow(parse_cigar("20M100D30M", 181L)$introns) == 0L)
stopifnot(nrow(parse_cigar("20M100N30M100N20M", 181L)$introns) == 2L)
b <- parse_cigar("5S10=2X8M3I100N20M2D5M", 181L)
stopifnot(identical(as.integer(b$introns[1, ]), c(200L, 300L)))
# Same-fragment paired ends both infer the same transcript strand.
stopifnot(all(infer_regtools_strand(c(99L, 147L), c(NA, NA), "RF") == "-"))
stopifnot(all(infer_regtools_strand(c(99L, 147L), c(NA, NA), "FR") == "+"))
stopifnot(infer_regtools_strand(0L, NA_character_, "XS") == "?")

s <- tempfile(fileext = ".sam"); make_demo_sam(s)
lines <- readLines(s); records <- lines[!startsWith(lines, "@")]
cfg <- default_demo_config()
# Filter flags and MAPQ as samtools would before parse_sam_records().
fields <- strsplit(records, "\t", fixed = TRUE)
keep <- vapply(fields, function(x) as.integer(x[[5]]) >= 20L &&
                 bitwAnd(as.integer(x[[2]]), 2820L) == 0L, logical(1))
p <- parse_sam_records(records[keep], cfg)
stopifnot(nrow(p$reads) == 4L, p$span_only_excluded == 1L, nrow(p$events) == 3L)
unlink(s)
cat("Core tests: PASS\n")


x <- parse_bam_contigs(c("@HD\tVN:1.6", "@SQ\tSN:chr2\tLN:242193529", "@SQ\tSN:chr1\tLN:248956422", "@SQ\tSN:chrX\tLN:156040895"))
stopifnot(identical(x$chrom, c("chr2","chr1","chrX")), choose_bam_contig(x) == "chr1",
          choose_bam_contig(x, "chrX") == "chrX", choose_bam_contig(x, "chrDemo") == "chr1")
y <- parse_bam_contigs(c("@SQ\tSN:2\tLN:200", "@SQ\tSN:1\tLN:100"))
stopifnot(choose_bam_contig(y) == "1")
z <- parse_bam_contigs("@SQ\tSN:scaffoldA\tLN:123")
stopifnot(choose_bam_contig(z) == "scaffoldA")
bad <- list(character(), "@SQ\tSN:chr1", "@SQ\tSN:chr1\tLN:no",
            c("@SQ\tSN:chr1\tLN:100", "@SQ\tSN:chr1\tLN:100"))
stopifnot(all(vapply(bad, function(h) inherits(try(parse_bam_contigs(h), silent=TRUE), "try-error"), logical(1))))
bam <- Sys.getenv("REGTOOLS_TEST_BAM")
if (nzchar(bam)) {
  actual <- read_bam_contigs(bam)
  stopifnot(all(c("chr1","chr22","chrX","chrY") %in% actual$chrom),
            !"chrDemo" %in% actual$chrom, choose_bam_contig(actual, "chrDemo") == "chr1")
  cat("Real BAM header choices: PASS; contigs=", nrow(actual), "\n", sep="")
}
cat("Contig selection tests: PASS\n")

rec <- function(cigar, seq, qual = paste(rep("I", nchar(seq)), collapse = ""), flag = 0L) {
  paste(c("synthetic", flag, "chrDemo", 100, 60, cigar, "*", 0, 0, seq, qual), collapse = "\t")
}
# Soft clips and insertions consume query, N and D consume reference only.
x <- count_rna_bases(c(rec("2S3M2I2M", "TTAAAGGCT"),
                       rec("3M2D2M", "AAACC"),
                       rec("3M2N2M", "AAACC"),
                       rec("5M", "AAAGT", flag = 16L),
                       rec("5M", "AAAAT", "III!I"),
                       rec("5M", "AAAAA", "*")), 103L, 20L)
stopifnot(x$count[x$base == "C"] == 1L, x$count[x$base == "G"] == 1L,
          sum(x$count) == 2L, attr(x, "excluded")[["deleted_D"]] == 1L,
          attr(x, "excluded")[["skipped_N"]] == 1L,
          attr(x, "excluded")[["low_baseq"]] == 1L,
          attr(x, "excluded")[["missing_quality"]] == 1L)
stopifnot(sum(count_rna_bases(character(), 103L)$count) == 0L)
stopifnot(sum(count_rna_bases(rec("5M", "AAAAA"), 110L)$count) == 0L)
cat("RNA base CIGAR, reverse-strand, BQ and missingness tests: PASS\n")

if (unit_only) {
  cat("Native integration tests: skipped (--unit)\n")
} else {
  cat("Native samtools / RegTools integration tests\n")
  r <- analyze_bam(cfg)
  m <- summarize_target(r, 201L, 300L, 5L)
  stopifnot(r$native_audit_passed, m$denominator_reads == 4L,
            m$junction_reads == 3L, m$junction_exact_count == 2L,
            m$junction_near_count == 1L, m$junction_per_100_reads == 50,
            abs(m$mean_read_depth - 200/350) < 1e-8,
            sum(r$depth$depth) == 200L, r$span_only_excluded == 1L)
  stopifnot(summarize_target(r, 201L, 300L, 0L)$junction_near_count == 0L)
  stopifnot(summarize_target(r, 201L, 300L, 1L)$junction_near_count == 0L)
  stopifnot(summarize_target(r, 201L, 300L, 2L)$junction_near_count == 1L)
  stopifnot(inherits(try(summarize_target(r, 101L, 300L, 5L), silent = TRUE), "try-error"))
  # No retained reads: no-junction is zero, but a zero-denominator rate is NA.
  cfg$mapq <- 61L
  z <- analyze_bam(cfg); mz <- summarize_target(z, 201L, 300L, 5L)
  stopifnot(mz$denominator_reads == 0L, mz$junction_exact_count == 0L,
            is.na(mz$junction_per_100_reads), sum(z$depth$depth) == 0L)
  cat("Native integration tests: PASS\n")
}

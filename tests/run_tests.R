# Run from the project directory: Rscript tests/run_tests.R
# Core tests need R; integration tests additionally need all environment.yml dependencies.
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

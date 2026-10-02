source("backend.R")
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

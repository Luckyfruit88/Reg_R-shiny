source("backend.R", encoding = "UTF-8")
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

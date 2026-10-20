# From project root: Rscript tests/test_variants.R [--unit]
# Synthetic fixtures contain no study data. Native checks require bcftools,
# samtools and RegTools; --unit explicitly runs only parser/aggregation checks,
# Model 1 gate tests and the Model 1 null-calibration simulations.
source("backend.R", encoding = "UTF-8")
source("variant_backend.R", encoding = "UTF-8")

expect_error <- function(expr, text = NULL) {
  e <- tryCatch({ force(expr); NULL }, error = function(e) e)
  stopifnot(inherits(e, "error"))
  if (!is.null(text)) stopifnot(grepl(text, conditionMessage(e), fixed = TRUE))
}

# Builds a Model 1 input: one row per sample, and per-sample exact reads for
# each junction key (zeros are left out, as in the sparse result table).
model1_fixture <- function(genotype, retained, counts) {
  ids <- sprintf("s%05d", seq_along(genotype))
  keys <- names(counts)
  sj <- do.call(rbind, lapply(keys, function(k) {
    hit <- counts[[k]] > 0
    data.frame(vcf_sample = ids[hit], genotype = genotype[hit], key = rep(k, sum(hit)),
               score = counts[[k]][hit], stringsAsFactors = FALSE)
  }))
  pos <- seq_along(keys) * 100L
  list(variant = data.frame(alt = "T"),
       samples = data.frame(vcf_sample = ids, genotype = genotype, status = "SUCCESS",
                            retained_reads = retained, stringsAsFactors = FALSE),
       junctions = data.frame(key = keys, chrom = "c", intron_start1 = pos + 1L, intron_end1 = pos + 50L,
                              strand = "+", stringsAsFactors = FALSE),
       sample_junctions = sj)
}

run_model1_gate_tests <- function() {
  g <- variant_association_gates()
  stopifnot(g$min_denominator_reads == 20, g$min_nonzero_samples == 10L, g$min_total_exact_reads == 20,
            g$min_mac == 5L, g$min_group_samples == 5L, g$min_groups == 2L)
  geno <- rep(c("0/0", "0/1", "1/1"), c(60L, 30L, 10L))
  dose <- rep(c(0, 1, 2), c(60L, 30L, 10L))
  status_of <- function(r, ...) variant_genotype_association(r, ...)$status

  # Singleton carrier with a private junction: one heterozygote, 7 reads.
  private <- numeric(100); private[61] <- 7
  a <- variant_genotype_association(model1_fixture(geno, rep(100, 100), list("c|0|1|+" = private)))
  stopifnot(a$status == "LOW_SUPPORT", is.na(a$p_value), is.na(a$q_value_bh), a$n_nonzero == 1L,
            a$total_exact_reads == 7)

  # LOW_SUPPORT: 9 supporting samples, or 15 supporting samples with 15 reads in total.
  nine <- c(rep(5, 9), numeric(91)); fifteen <- c(rep(1, 15), numeric(85))
  stopifnot(status_of(model1_fixture(geno, rep(100, 100), list("c|0|1|+" = nine))) == "LOW_SUPPORT",
            status_of(model1_fixture(geno, rep(100, 100), list("c|0|1|+" = fifteen))) == "LOW_SUPPORT")

  # LOW_MAC: mac 2, and mac 6 carried by one 1/1 group of 3 (only one group reaches 5).
  rare <- rep(c("0/0", "0/1"), c(98L, 2L))
  hom_only <- rep(c("0/0", "1/1"), c(97L, 3L))
  busy <- rep(3, 100)
  stopifnot(status_of(model1_fixture(rare, rep(100, 100), list("c|0|1|+" = busy))) == "LOW_MAC",
            status_of(model1_fixture(hom_only, rep(100, 100), list("c|0|1|+" = busy))) == "LOW_MAC",
            variant_genotype_association(model1_fixture(hom_only, rep(100, 100), list("c|0|1|+" = busy)))$mac == 6L)

  # TOO_FEW_SAMPLES precedes every other gate.
  stopifnot(status_of(model1_fixture(c("0/0", "0/1"), c(100, 100), list("c|0|1|+" = c(5, 9)))) == "TOO_FEW_SAMPLES")

  # CONSTANT_OUTCOME: every sample has 1 read per 100.
  stopifnot(status_of(model1_fixture(geno, rep(100, 100), list("c|0|1|+" = rep(1, 100)))) == "CONSTANT_OUTCOME")

  # Rates constant within dosage groups and exactly linear (1, 6, 11 per 900):
  # floating-point residue must not be reported as an OK fit with p ~ 0.
  linear <- c(1, 6, 11)[dose + 1]
  pf <- variant_genotype_association(model1_fixture(geno, rep(900, 100), list("c|0|1|+" = linear)))
  stopifnot(pf$status == "PERFECT_FIT", is.na(pf$se), is.na(pf$p_value), is.na(pf$q_value_bh), is.finite(pf$beta))
  # Exactly linear per-sample rates with unequal denominators.
  den <- rep(c(100, 300), 50)
  exact <- (2 + 3 * dose) * den / 100
  stopifnot(status_of(model1_fixture(geno, den, list("c|0|1|+" = exact))) == "PERFECT_FIT")
  # Group-constant but not linear rates (1, 6, 7) keep residual variance and are fitted.
  stopifnot(status_of(model1_fixture(geno, rep(900, 100), list("c|0|1|+" = c(1, 6, 7)[dose + 1]))) == "OK")

  # Per-sample denominator gate: 10 samples with 10 reads are excluded from every pair.
  low_den <- rep(100, 100); low_den[1:10] <- 10
  set.seed(11); noisy <- stats::rpois(100, 5)
  ld <- variant_genotype_association(model1_fixture(geno, low_den, list("c|0|1|+" = noisy)))
  acc <- attr(ld, "sample_accounting")
  stopifnot(ld$n == 90L, ld$n_dosage0 == 50L, acc$below_min_denominator == 10L, acc$analysed == 90L,
            abs(ld$mean_denominator_dosage0 - 100) < 1e-12)

  # BH runs over OK rows only; gated rows keep NA q-values.
  set.seed(12)
  mixed <- list("c|0|1|+" = stats::rpois(100, 4) + 2 * dose, "c|1|2|+" = stats::rpois(100, 4),
                "c|2|3|+" = c(rep(2, 5), numeric(95)))
  mx <- variant_genotype_association(model1_fixture(geno, rep(100, 100), mixed))
  ok <- mx$status == "OK"
  stopifnot(sum(ok) == 2L, all(mx$status[!ok] == "LOW_SUPPORT"), all(is.na(mx$q_value_bh[!ok])),
            all(abs(mx$q_value_bh[ok] - stats::p.adjust(mx$p_value[ok], "BH")) < 1e-15))

  # A strong effect whose p-value underflows double precision keeps its magnitude.
  set.seed(13)
  strong <- stats::rpois(2000, 50) + 400 * rep(c(0, 1, 2), c(1000, 800, 200))
  st <- variant_genotype_association(model1_fixture(rep(c("0/0", "0/1", "1/1"), c(1000L, 800L, 200L)),
                                                    rep(1000, 2000), list("c|0|1|+" = strong)))
  stopifnot(st$status == "OK", st$log10_p_value < -300, st$p_value < 1e-300,
            variant_format_p(st$p_value) == "< 1e-300")
  stopifnot(identical(variant_format_p(c(0, 1e-310, 0.5, NA, 2e-5)),
                      c("< 1e-300", "< 1e-300", "0.5", NA, "2e-05")))
  cat("Model 1 gates, perfect-fit tolerance and p-value display: PASS\n")
}

# Acceptance simulations for the gates (deterministic seeds).
run_model1_calibration <- function() {
  set.seed(20261020)
  n <- 4700L; maf <- 0.05
  geno <- sample(c("0/0", "0/1", "1/1"), n, TRUE, c((1 - maf)^2, 2 * maf * (1 - maf), maf^2))
  retained <- stats::rpois(n, 1000) + 1
  # Null 1: 1000 junctions, each seen in one random sample with one read.
  carriers <- sample(n, 1000L, TRUE)
  singles <- stats::setNames(lapply(carriers, function(i) { v <- numeric(n); v[i] <- 1; v }),
                             sprintf("c|%d|%d|+", seq_len(1000L), seq_len(1000L) + 1L))
  s1 <- variant_genotype_association(model1_fixture(geno, retained, singles))
  false_q <- sum(s1$q_value_bh < 0.05, na.rm = TRUE)
  # Null 2: 300 well-supported junctions whose expected rate (2 per 100
  # denominator reads) does not depend on genotype. Counts scale with each
  # sample's denominator; otherwise the shared denominator vector would itself
  # carry a chance genotype association into every pair.
  m <- 1000L; maf2 <- 0.3
  geno2 <- sample(c("0/0", "0/1", "1/1"), m, TRUE, c((1 - maf2)^2, 2 * maf2 * (1 - maf2), maf2^2))
  den2 <- stats::rpois(m, 1000) + 1
  supported <- stats::setNames(lapply(seq_len(300L), function(i) stats::rpois(m, 0.02 * den2)),
                               sprintf("c|%d|%d|+", seq_len(300L), seq_len(300L) + 1L))
  s2 <- variant_genotype_association(model1_fixture(geno2, den2, supported))
  ks <- suppressWarnings(stats::ks.test(s2$p_value[s2$status == "OK"], "punif"))$p.value
  cat(sprintf("Singleton null: %d of %d pairs with q < 0.05 (statuses: %s)\n", false_q, nrow(s1),
              paste(names(table(s1$status)), table(s1$status), collapse = ", ")))
  cat(sprintf("Well-supported null: %d OK pairs, KS p = %.3f\n", sum(s2$status == "OK"), ks))
  stopifnot(false_q == 0L, sum(s2$status == "OK") == 300L, ks > 0.05)
  cat("Model 1 null calibration: PASS\n")
}

run_variant_tests <- function(native = TRUE) {
  q <- parse_variant_query("chr1:101:A:C,G")
  stopifnot(q$chrom == "chr1", q$pos1 == 101L, q$ref == "A", q$alt == "C,G")
  for (bad in c("chr1:0", "chr1:1.1", "chr1:1;touch /tmp/bad", "../chr1:1", "chr1:1:A:C;bad"))
    expect_error(parse_variant_query(bad))
  gt <- lapply(c("0/0", "1|0", "2/1", "1", "./.", "0/.", ".|1", "0/1/2"), canonical_dna_gt, n_alt = 2)
  stopifnot(gt[[2]]$genotype == "0/1", gt[[2]]$phased,
            gt[[3]]$genotype == "1/2", gt[[4]]$ploidy == 1L,
            gt[[5]]$call_status == "NO_CALL", gt[[6]]$genotype == "PARTIAL_CALL:0/.",
            gt[[7]]$genotype == "PARTIAL_CALL:1/.", gt[[8]]$ploidy == 3L)
  expect_error(canonical_dna_gt("0/2", 1), "exceeds")
  expect_error(canonical_dna_gt("1x0", 1))
  expect_error(canonical_dna_gt("999999999999999999999", 1))

  # Independent arithmetic: a successful zero is in the mean denominator; a
  # failed BAM is not. Entirely failed and unlinked groups remain NA.
  cfg <- default_demo_config(); cfg$start1 <- 190L; cfg$end1 <- 191L
  calls <- data.frame(vcf_sample = c("a", "b", "c", "d", "e"),
    genotype = c("0/0", "0/0", "0/0", "1/1", "NO_CALL:./."),
    call_status = c("CALLED", "CALLED", "CALLED", "CALLED", "NO_CALL"), linked = c(TRUE, TRUE, TRUE, TRUE, FALSE))
  selected <- calls[1:4, ]; selected$status <- c("SUCCESS", "SUCCESS", "FAILED", "FAILED")
  junction <- data.frame(key = "chrDemo|200|300|+", chrom = "chrDemo", intron_start1 = 201L,
                        intron_end1 = 300L, strand = "+", score = 4L)
  r1 <- list(depth = data.frame(depth = c(2L, 4L)), junctions = junction,
             rna_bases = data.frame(base = c("A", "C", "G", "T", "N"), count = c(2L, 0L, 0L, 0L, 0L)))
  attr(r1$rna_bases, "excluded") <- c(skipped_N = 0L, deleted_D = 0L, low_baseq = 0L,
                                     missing_sequence = 0L, missing_quality = 0L)
  r2 <- r1; r2$depth$depth <- c(0L, 0L); r2$junctions <- junction[FALSE, ]; r2$rna_bases$count[] <- 0L
  agg <- variant_aggregate(NULL, calls, selected, list(r1, r2, NULL, NULL), cfg, TRUE)
  s <- agg$group_summary
  stopifnot(s$analyzed_n[s$genotype == "0/0"] == 2L, s$failed_n[s$genotype == "0/0"] == 1L,
            all(agg$depth$mean_depth[agg$depth$genotype == "0/0"] == c(1, 2)),
            all(is.na(agg$depth$mean_depth[agg$depth$genotype == "1/1"])),
            agg$junctions$count_mean[agg$junctions$genotype == "0/0"] == 2,
            is.na(agg$junctions$count_sum[agg$junctions$genotype == "1/1"]),
            all(is.na(agg$rna_bases$count[agg$rna_bases$genotype == "NO_CALL:./."])))
  # Sparse per-sample rows: only the successful sample with support appears; the
  # successful zero (b) and failed BAMs (c, d) have no rows.
  sj <- agg$sample_junctions
  stopifnot(identical(names(sj), c("vcf_sample", "genotype", "key", "score")), nrow(sj) == 1L,
            sj$vcf_sample == "a", sj$genotype == "0/0", sj$key == junction$key, sj$score == 4,
            sum(sj$score) == agg$junctions$count_sum[agg$junctions$genotype == "0/0"])
  tampered <- sj; tampered$score <- 5
  expect_error(variant_check_sample_junctions(tampered, agg$junctions))
  expect_error(variant_check_sample_junctions(sj[FALSE, ], agg$junctions))
  stopifnot(isTRUE(variant_check_sample_junctions(variant_empty_sample_junctions(),
                                                  agg$junctions[FALSE, ])))
  # Model 1 matches lm() on the independently computed per-sample rates. Strands
  # of one boundary pair are summed; a zero denominator and a failed BAM are excluded.
  stopifnot(identical(variant_genotype_dosage(c("0/0", "0/1", "1/1"), 1L), c(0, 1, 2)))
  stopifnot(all(is.na(variant_genotype_dosage(c("0", "1", "0/1/1", "PARTIAL_CALL:0/.", "NO_CALL:./."), 1L))),
            all(is.na(variant_genotype_dosage(c("0/0", "0/1", "1/2"), 2L))))
  assoc_result <- list(variant = data.frame(alt = "C"),
    samples = data.frame(vcf_sample = letters[1:10], genotype = rep(c("0/0", "0/1", "1/1", "1/1"), c(3, 3, 3, 1)),
      status = c(rep("SUCCESS", 9), "FAILED"), retained_reads = c(100, 200, 100, 100, 100, 50, 100, 0, 200, NA),
      stringsAsFactors = FALSE),
    junctions = data.frame(key = c("c|9|20|+", "c|9|20|-", "c|29|40|+"), chrom = "c",
      intron_start1 = c(10L, 10L, 30L), intron_end1 = c(20L, 20L, 40L), strand = c("+", "-", "+"),
      stringsAsFactors = FALSE),
    sample_junctions = data.frame(vcf_sample = c("a", "b", "c", "c", "d", "e", "f", "g", "i", "a"),
      genotype = "", key = c(rep("c|9|20|+", 2), "c|9|20|-", "c|9|20|+", "c|9|20|+", "c|9|20|-",
        "c|9|20|+", "c|9|20|+", "c|9|20|+", "c|29|40|+"),
      score = c(10, 30, 5, 5, 20, 25, 15, 40, 90, 3), stringsAsFactors = FALSE))
  # Gates relaxed to their minimum so the small fixture checks the arithmetic.
  open_gates <- variant_association_gates(min_denominator_reads = 1, min_nonzero_samples = 1L,
    min_total_exact_reads = 1, min_mac = 1L, min_group_samples = 1L)
  assoc <- variant_genotype_association(assoc_result, gates = open_gates)
  fit <- summary(lm(rate ~ dosage, data.frame(rate = c(10, 15, 10, 20, 25, 30, 40, 45), dosage = c(0, 0, 0, 1, 1, 1, 2, 2))))$coefficients
  p1 <- assoc[assoc$intron_start1 == 10L, ]
  stopifnot(nrow(assoc) == 2L, p1$status == "OK", p1$n == 8L, p1$strands == "+,-",
            c(p1$n_dosage0, p1$n_dosage1, p1$n_dosage2) == c(3L, 3L, 2L),
            p1$mac == 7L, p1$n_nonzero == 8L, p1$total_exact_reads == 240,
            abs(p1$mean_denominator_dosage2 - 150) < 1e-12,
            abs(p1$mean_rate_dosage2 - 42.5) < 1e-12,
            abs(p1$beta - fit["dosage", "Estimate"]) < 1e-10, abs(p1$se - fit["dosage", "Std. Error"]) < 1e-10,
            abs(p1$p_value - fit["dosage", "Pr(>|t|)"]) < 1e-12,
            abs(p1$log10_p_value - log10(fit["dosage", "Pr(>|t|)"])) < 1e-10,
            all(is.finite(assoc$q_value_bh)), assoc$p_value[[1L]] <= assoc$p_value[[2L]])
  # The same fixture under the default gates: no dosage group reaches 5 samples.
  stopifnot(all(variant_genotype_association(assoc_result)$status == "LOW_MAC"),
            all(is.na(variant_genotype_association(assoc_result)$p_value)))
  single <- assoc_result; single$samples$genotype <- "0/1"
  stopifnot(all(variant_genotype_association(single)$status == "SINGLE_DOSAGE"),
            all(is.na(variant_genotype_association(single)$p_value)))
  legacy <- assoc_result; legacy$sample_junctions <- NULL
  expect_error(variant_genotype_association(legacy), "rerun")
  run_model1_gate_tests()
  run_model1_calibration()
  unsupported <- variant_aggregate(NULL, calls, selected, list(r1, r2, NULL, NULL), cfg, FALSE)
  stopifnot(all(is.na(unsupported$rna_bases$count)), all(unsupported$rna_bases$evidence_status == "UNSUPPORTED_NON_SNV"))
  missing_base <- r2
  attr(missing_base$rna_bases, "excluded")[["missing_quality"]] <- 2L
  only_missing <- variant_aggregate(NULL, calls, selected, list(missing_base, r2, NULL, NULL), cfg, TRUE)
  mb <- only_missing$rna_bases[only_missing$rna_bases$genotype == "0/0", ]
  stopifnot(all(mb$count == 0), all(mb$excluded_missing_quality == 2L),
            all(mb$evidence_status == "NO_VALID_BASES_WITH_MISSING_DATA"), all(mb$valid_base_count == 0))
  cat("Variant parser, native GT preservation and missingness arithmetic: PASS\n")
  if (!native) return(invisible(TRUE))

  executables <- Sys.which(c("bcftools", "samtools", "regtools"))
  if (any(!nzchar(executables))) stop("Native tests require bcftools, samtools and regtools; use --unit only for explicitly partial verification.")
  work <- tempfile("regshiny_variant_test_"); dir.create(work, mode = "0700")
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  run <- function(tool, args) {
    p <- processx::run(unname(executables[[tool]]), as.character(args), timeout = 120, error_on_status = FALSE)
    if (p$status != 0L) stop(tool, ": ", p$stderr)
    p$stdout
  }
  sam <- file.path(work, "demo.sam"); make_demo_sam(sam)
  bam <- file.path(work, "template.bam")
  run("samtools", c("sort", "-o", bam, sam)); run("samtools", c("index", bam))
  ids <- paste0("S", 1:8)
  bams <- file.path(work, paste0(ids, ".bam"))
  for (i in 1:5) {
    stopifnot(file.copy(bam, bams[[i]]), file.copy(paste0(bam, ".bai"), paste0(bams[[i]], ".bai")))
  }
  # A valid, empty BAM represents observable zero coverage.
  empty_sam <- file.path(work, "empty.sam")
  writeLines(c("@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chrDemo\tLN:1000"), empty_sam)
  run("samtools", c("sort", "-o", bams[[6]], empty_sam)); run("samtools", c("index", bams[[6]]))
  # S7 intentionally missing; S8 has a wrong contig length.
  wrong_sam <- file.path(work, "wrong.sam")
  writeLines(c("@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chrDemo\tLN:1001"), wrong_sam)
  run("samtools", c("sort", "-o", bams[[8]], wrong_sam)); run("samtools", c("index", bams[[8]]))
  raw_vcf <- file.path(work, "synthetic.vcf"); vcf <- paste0(raw_vcf, ".gz")
  record <- function(pos, ref, alt, gt, id = ".") paste(c("chrDemo", pos, id, ref, alt, ".", "PASS", ".", "GT", gt), collapse = "\t")
  writeLines(c("##fileformat=VCFv4.2", "##reference=synthetic_fixture",
    "##contig=<ID=chrDemo,length=1000>", "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", ids), collapse = "\t"),
    record(190, "A", "C", c("0/0", "0|1", "1/1", "./.", "1", "0/1", "0/1", "1/1"), "two_alleles"),
    record(190, "A", "G", rep("0/0", 8), "same_position_other_allele"),
    record(191, "A", "C,G", c("1/2", "2|1", "0/2", "0/.", "./.", "0/0", "1/1", "2")),
    record(192, "A", "C", c(rep("0/0", 6), "0/1", "0/1")),
    record(193, "A", "AT", rep("0/0", 8)),
    record(194, "A", "C", rep("./.", 8))), raw_vcf)
  run("bcftools", c("view", "-Oz", "-o", vcf, raw_vcf)); run("bcftools", c("index", "-t", vcf))
  registry <- file.path(work, "vcfs.tsv"); manifest <- file.path(work, "samples.tsv")
  vr <- data.frame(chrom = "chrDemo", vcf = vcf, build = "synthetic_v1")
  sm <- data.frame(vcf_sample = ids, bam = bams, build = "synthetic_v1")
  write_tsv <- function(x, path) utils::write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE)
  write_tsv(vr, registry); write_tsv(sm, manifest)
  resources <- load_variant_resources(registry, manifest)
  variants <- lookup_variants(resources, "chrDemo:190")
  stopifnot(nrow(variants) == 2L, setequal(variants$alt, c("C", "G")),
            all(variants$contig_length == 1000), nrow(lookup_variants(resources, "chrDemo:999")) == 0L)
  expect_error(lookup_variants(resources, "Demo:190"), "not configured")
  expect_error(lookup_variants(resources, "chrDemo:1001"), "exceeds")
  variant <- lookup_variants(resources, "chrDemo:190:A:C")
  cfg <- default_demo_config(); cfg$demo <- FALSE
  expect_error(analyze_variant(resources, sampling_mode = "preview", variants, cfg), "Select exactly one")
  result <- analyze_variant(resources, sampling_mode = "preview", variant, cfg, max_per_group = 3L)
  g <- result$group_summary
  stopifnot(nrow(g) == 5L, result$genotypes$raw_gt[result$genotypes$vcf_sample == "S2"] == "0|1",
            result$genotypes$genotype[result$genotypes$vcf_sample == "S2"] == "0/1",
            g$analyzed_n[g$genotype == "0/1"] == 2L, g$failed_n[g$genotype == "0/1"] == 1L,
            g$failed_n[g$genotype == "1/1"] == 1L,
            grepl("contig lengths differ", result$samples$error[result$samples$vcf_sample == "S8"], fixed = TRUE),
            result$samples$mean_depth[result$samples$vcf_sample == "S6"] == 0,
            is.na(result$samples$mean_depth[result$samples$vcf_sample == "S7"]))
  d <- result$depth[result$depth$genotype == "0/1", ]
  stopifnot(abs(mean(d$mean_depth) - (200 / 350 / 2)) < 1e-10)
  j <- result$junctions[result$junctions$genotype == "0/1" & result$junctions$intron_start1 == 201L, ]
  stopifnot(nrow(j) == 1L, j$count_sum == 2, j$count_mean == 1,
            result$rna_bases$count[result$rna_bases$genotype == "0/1" & result$rna_bases$base == "A"] == 3)
  # Bounded deterministic subset is explicit, with no failure replacement.
  capped <- analyze_variant(resources, sampling_mode = "preview", variant, cfg, max_per_group = 1L)
  stopifnot(capped$samples$vcf_sample[capped$samples$genotype == "0/1"] == "S2",
            capped$group_summary$available_n[capped$group_summary$genotype == "0/1"] == 3L,
            capped$group_summary$selected_n[capped$group_summary$genotype == "0/1"] == 1L)
  multi <- lookup_variants(resources, "chrDemo:191:A:C,G")
  gm <- variant_get_genotypes(resources, multi)
  stopifnot(all(gm$genotype[gm$vcf_sample %in% c("S1", "S2")] == "1/2"),
            gm$raw_gt[gm$vcf_sample == "S2"] == "2|1", gm$call_status[gm$vcf_sample == "S4"] == "PARTIAL_CALL")
  failed <- analyze_variant(resources, sampling_mode = "preview", lookup_variants(resources, "chrDemo:192"), cfg, max_per_group = 2L)
  stopifnot(all(is.na(failed$depth$mean_depth[failed$depth$genotype == "0/1"])),
            all(is.na(failed$junctions$count_sum[failed$junctions$genotype == "0/1"])),
            all(is.na(failed$rna_bases$count[failed$rna_bases$genotype == "0/1"])))
  indel <- analyze_variant(resources, sampling_mode = "preview", lookup_variants(resources, "chrDemo:193"), cfg, max_per_group = 1L)
  stopifnot(all(is.na(indel$rna_bases$count)), indel$provenance$rna_evidence_status == "UNSUPPORTED_NON_SNV")
  missing_variant <- lookup_variants(resources, "chrDemo:194")
  missing_calls <- variant_get_genotypes(resources, missing_variant)
  stopifnot(all(missing_calls$call_status == "NO_CALL"))
  expect_error(analyze_variant(resources, sampling_mode = "preview", missing_variant, cfg, max_per_group = 1L), "No called DNA genotype")
  # Manifest validation is independent of expensive native counting.
  sm_bad <- sm; sm_bad$build[[1L]] <- "other_build"; write_tsv(sm_bad, manifest)
  expect_error(load_variant_resources(registry, manifest), "same known reference")
  expect_error(lookup_variants(resources, "chrDemo:190"), "manifest changed")
  sm_bad <- sm; sm_bad$vcf_sample[[2L]] <- "S1"; write_tsv(sm_bad, manifest)
  expect_error(load_variant_resources(registry, manifest), "unique")
  sm_bad <- sm; sm_bad$bam[[2L]] <- sm_bad$bam[[1L]]; write_tsv(sm_bad, manifest)
  expect_error(load_variant_resources(registry, manifest), "multiple VCF samples")
  write_tsv(sm, manifest)
  cat("Native indexed VCF + RNA BAM genotype/junction/depth/base integration: PASS\n")
  invisible(TRUE)
}

run_variant_tests(native = !"--unit" %in% commandArgs(trailingOnly = TRUE))

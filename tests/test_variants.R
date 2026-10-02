# From project root: Rscript tests/test_variants.R [--unit]
# Synthetic fixtures contain no study data. Native checks require bcftools,
# samtools and RegTools; --unit explicitly runs only parser/aggregation checks.
source("backend.R", encoding = "UTF-8")
source("variant_backend.R", encoding = "UTF-8")

expect_error <- function(expr, text = NULL) {
  e <- tryCatch({ force(expr); NULL }, error = function(e) e)
  stopifnot(inherits(e, "error"))
  if (!is.null(text)) stopifnot(grepl(text, conditionMessage(e), fixed = TRUE))
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

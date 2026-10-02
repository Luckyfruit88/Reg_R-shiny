# Rscript tests/test_all_samples.R [--native]
# Deterministic mock orchestration tests exercise full-size selection, streaming,
# atomic progress and interruption/resume without requiring study data or tools.
# --native additionally processes 36 synthetic BAMs and an indexed native VCF.
source("backend.R", encoding = "UTF-8")
source("variant_backend.R", encoding = "UTF-8")

all_expect_error <- function(expr, text = NULL) {
  e <- tryCatch({ force(expr); NULL }, error = function(e) e)
  stopifnot(inherits(e, "error"))
  if (!is.null(text)) stopifnot(grepl(text, conditionMessage(e), fixed = TRUE))
}

all_equal_table <- function(a, b, columns) {
  a <- a[, columns, drop = FALSE]; b <- b[, columns, drop = FALSE]
  key <- intersect(c("genotype", "chrom", "pos1", "intron_start1", "intron_end1", "strand", "base"), columns)
  if (nrow(a)) a <- a[do.call(order, a[key]), , drop = FALSE]
  if (nrow(b)) b <- b[do.call(order, b[key]), , drop = FALSE]
  rownames(a) <- rownames(b) <- NULL
  stopifnot(isTRUE(all.equal(a, b, check.attributes = FALSE, tolerance = 1e-12)))
}

run_all_mock_tests <- function() {
  work <- tempfile("regshiny_all_mock_"); dir.create(work, mode = "0700")
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  replaced <- c("lookup_variants", "variant_vcf_info", "variant_get_genotypes", "variant_run", "read_bam_contigs", "analyze_bam")
  originals <- mget(replaced, envir = .GlobalEnv)
  on.exit(list2env(originals, envir = .GlobalEnv), add = TRUE)
  old_slots <- Sys.getenv("NSLOTS", unset = NA_character_); Sys.unsetenv("NSLOTS")
  on.exit(if (is.na(old_slots)) Sys.unsetenv("NSLOTS") else Sys.setenv(NSLOTS = old_slots), add = TRUE)
  ids <- sprintf("S%02d", 1:42)
  bams <- file.path(work, paste0(ids, ".bam"))
  stopifnot(all(file.create(bams)))
  vcf <- file.path(work, "synthetic.vcf.gz"); file.create(vcf, paste0(vcf, ".tbi"))
  manifest <- file.path(work, "samples.tsv"); registry <- file.path(work, "vcfs.tsv")
  write_tsv <- function(x, path) utils::write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE)
  write_tsv(data.frame(chrom = "chrDemo", vcf = vcf, build = "synthetic_v1"), registry)
  write_tsv(data.frame(vcf_sample = ids[1:41], bam = bams[1:41], build = "synthetic_v1"), manifest)
  resources <- load_variant_resources(registry, manifest)
  variant <- data.frame(record_id = "synthetic_v1:chrDemo:190:A:C", chrom = "chrDemo", pos1 = 190L,
    id = ".", ref = "A", alt = "C", qual = ".", filter = "PASS", build = "synthetic_v1", vcf = vcf, contig_length = 1000)
  calls <- data.frame(vcf_sample = ids, raw_gt = c(rep(c("0/0", "0/1", "1/1"), each = 13), "./.", "0/.", "1/1"),
    genotype = c(rep(c("0/0", "0/1", "1/1"), each = 13), "NO_CALL:./.", "PARTIAL_CALL:0/.", "1/1"),
    call_status = c(rep("CALLED", 39), "NO_CALL", "PARTIAL_CALL", "CALLED"),
    ploidy = 2L, phased = FALSE, linked = c(rep(TRUE, 41), FALSE), bam = bams)
  calls$bam[[42]] <- NA_character_
  ctx <- new.env(); ctx$attempts <- 0L; ctx$interrupt_at <- Inf; ctx$mutate_vcf <- FALSE
  fixture <- function(cfg) {
    id <- as.integer(sub("S([0-9]+)\\.bam$", "\\1", basename(cfg$bam)))
    if (id == 3L || id %in% 27:39) stop("Synthetic unreadable BAM")
    zero <- id == 2L
    j <- data.frame(key = c("chrDemo|200|300|+", "chrDemo|202|302|+"), chrom = "chrDemo",
      intron_start1 = c(201L, 203L), intron_end1 = c(300L, 302L), strand = "+", score = c(2L, 1L))
    if (zero) j <- j[FALSE, , drop = FALSE]
    b <- data.frame(base = c("A", "C", "G", "T", "N"), count = c(if (zero) 0L else 2L, 0L, 0L, 0L, 0L))
    attr(b, "excluded") <- c(skipped_N = 0L, deleted_D = 0L, low_baseq = 0L, missing_sequence = 0L, missing_quality = 0L)
    list(config = cfg, reads = data.frame(read_id = seq_len(if (zero) 0L else 4L)), events = data.frame(large = 1:10),
      depth = data.frame(chrom = cfg$chrom, pos1 = seq.int(cfg$start1, cfg$end1), depth = if (zero) 0 else 2),
      junctions = j, rna_bases = b, native_audit_passed = TRUE, candidate_alignments = if (zero) 0L else 4L,
      span_only_excluded = 0L, warnings = character(), log = "synthetic command",
      versions = list(samtools = "synthetic", regtools = "synthetic", regtools_extract_help = "omit this repeated help"),
      source_metadata = list(size_bytes = 0), completed_at = "synthetic fixed time")
  }
  assign("lookup_variants", function(resources, query) { variant_check_resources(resources); variant }, .GlobalEnv)
  assign("variant_vcf_info", function(path, chrom) list(state = variant_file_state(c(path, paste0(path, ".tbi"))),
    contig_length = 1000, reference_header = "synthetic", warnings = character(), log = "synthetic header"), .GlobalEnv)
  assign("variant_get_genotypes", function(resources, variant) calls, .GlobalEnv)
  assign("variant_run", function(...) list(text = "synthetic bcftools", log = "synthetic", stderr = ""), .GlobalEnv)
  assign("read_bam_contigs", function(bam) data.frame(chrom = "chrDemo", length_bp = 1000), .GlobalEnv)
  assign("analyze_bam", function(cfg) {
    ctx$attempts <- ctx$attempts + 1L
    if (ctx$attempts == ctx$interrupt_at)
      stop(structure(list(message = "Synthetic abrupt interruption"), class = c("interrupt", "condition")))
    if (ctx$mutate_vcf) { cat("changed", file = vcf, append = TRUE); ctx$mutate_vcf <- FALSE }
    fixture(cfg)
  }, .GlobalEnv)
  cfg <- default_demo_config(); cfg$demo <- FALSE
  progress <- file.path(work, "progress.json"); checkpoints <- file.path(work, "checkpoints")
  full <- analyze_variant(resources, variant, cfg, progress_file = progress, checkpoint_dir = checkpoints)
  stopifnot(nrow(full$samples) == 39L, ctx$attempts == 39L,
            all(full$group_summary$selected_n[full$group_summary$call_status == "CALLED"] == 13L),
            all(full$group_summary$analyzed_n[full$group_summary$call_status != "CALLED"] == 0L),
            sum(full$samples$status == "SUCCESS") == 25L,
            sum(full$samples$status == "FAILED") == 14L,
            full$group_summary$analyzed_n[full$group_summary$genotype == "0/0"] == 12L,
            all(is.na(full$depth$mean_depth[full$depth$genotype == "1/1"])),
            all(is.na(full$junctions$count_sum[full$junctions$genotype == "1/1"])),
            abs(full$depth$mean_depth[full$depth$genotype == "0/0"][[1]] - 22/12) < 1e-12,
            all(vapply(Filter(Negate(is.null), full$sample_results), function(z)
              !any(c("reads", "events", "depth", "junctions", "log") %in% names(z)), logical(1))))
  p <- jsonlite::read_json(progress, simplifyVector = TRUE)
  stopifnot(p$state == "COMPLETE", p$total == 39L, p$completed == 39L, p$success == 25L,
            p$failed == 14L, p$incomplete_gt_n == 2L, p$missing_link_n == 1L,
            sum(p$group_summary$completed) == 39L)
  # Independent reference aggregation holds full fixtures only inside this test.
  reference <- lapply(full$samples$bam, function(path) {
    c <- cfg; c$bam <- path
    tryCatch(fixture(c), error = function(e) NULL)
  })
  expected <- variant_aggregate(variant, calls, full$samples, reference, cfg, TRUE)
  for (field in c("depth", "junctions", "rna_bases")) all_equal_table(full[[field]], expected[[field]], names(expected[[field]]))
  # Preview remains explicitly capped and unchanged; full mode ignores the cap.
  preview <- analyze_variant(resources, variant, cfg, sampling_mode = "preview", max_per_group = 3L)
  stopifnot(nrow(preview$samples) == 9L)
  legacy_preview <- analyze_variant(resources, variant, cfg, max_per_group = 3L)
  stopifnot(nrow(legacy_preview$samples) == 9L,
            legacy_preview$provenance$selection$mode == "preview")
  explicit_all <- analyze_variant(resources, variant, cfg, sampling_mode = "all", max_per_group = 3L)
  stopifnot(nrow(explicit_all$samples) == 39L,
            explicit_all$provenance$selection$mode == "all")
  all_expect_error(analyze_variant(resources, variant, cfg, sampling_mode = "preview", max_per_group = 11L), "integer")
  before <- ctx$attempts
  resumed <- analyze_variant(resources, variant, cfg, checkpoint_dir = checkpoints, progress_file = progress, max_per_group = 1L)
  stopifnot(ctx$attempts == before, resumed$provenance$checkpoint$reused == 39L,
            all(resumed$samples$checkpoint_reused), nrow(resumed$samples) == 39L)
  all_equal_table(full$depth, resumed$depth, names(full$depth))
  # One corrupt contribution is recomputed, without duplicating other samples.
  checkpoint_files <- list.files(checkpoints, pattern = "^sample_", full.names = TRUE)
  writeLines("corrupted synthetic checkpoint", checkpoint_files[[1]])
  repaired <- analyze_variant(resources, variant, cfg, checkpoint_dir = checkpoints)
  stopifnot(ctx$attempts == before + 1L, repaired$provenance$checkpoint$reused == 38L,
            repaired$provenance$checkpoint$invalidated == 1L)
  all_equal_table(full$depth, repaired$depth, names(full$depth))
  # Changed config, BAM stat and code identity all reject stale checkpoints.
  changed <- cfg; changed$mapq <- cfg$mapq + 1L
  all_expect_error(analyze_variant(resources, variant, changed, checkpoint_dir = checkpoints), "fingerprint mismatch")
  old_mtime <- file.info(bams[[1]])$mtime
  Sys.setFileTime(bams[[1]], old_mtime + 5)
  all_expect_error(analyze_variant(resources, variant, cfg, checkpoint_dir = checkpoints), "fingerprint mismatch")
  Sys.setFileTime(bams[[1]], old_mtime)
  code_identity <- variant_code_identity
  assign("variant_code_identity", function() list(changed = "synthetic changed code"), .GlobalEnv)
  all_expect_error(analyze_variant(resources, variant, cfg, checkpoint_dir = checkpoints), "fingerprint mismatch")
  assign("variant_code_identity", code_identity, .GlobalEnv)
  # Interruption after two committed samples resumes from those two only.
  interrupted_dir <- file.path(work, "interrupted")
  ctx$interrupt_at <- ctx$attempts + 3L
  interrupted <- tryCatch({ analyze_variant(resources, variant, cfg, checkpoint_dir = interrupted_dir, progress_file = progress); FALSE },
                          interrupt = function(e) TRUE)
  stopifnot(interrupted, length(list.files(interrupted_dir, pattern = "^sample_")) == 2L,
            jsonlite::read_json(progress)$state == "FAILED")
  ctx$interrupt_at <- Inf; before <- ctx$attempts
  recovered <- analyze_variant(resources, variant, cfg, checkpoint_dir = interrupted_dir, progress_file = progress)
  stopifnot(ctx$attempts - before == 37L, recovered$provenance$checkpoint$reused == 2L)
  all_equal_table(full$depth, recovered$depth, names(full$depth))
  # Explicit resume=FALSE reruns all selected samples, including failed ones.
  before <- ctx$attempts
  retried <- analyze_variant(resources, variant, cfg, checkpoint_dir = interrupted_dir, resume = FALSE)
  stopifnot(ctx$attempts - before == 39L, retried$provenance$checkpoint$reused == 0L)
  # Worker admission and parallel results are deterministic and equivalent.
  Sys.setenv(NSLOTS = "2")
  all_expect_error(analyze_variant(resources, variant, cfg, workers = 2L), "NSLOTS")
  Sys.setenv(NSLOTS = "3")
  parallel_result <- analyze_variant(resources, variant, cfg, workers = 2L)
  stopifnot(parallel_result$provenance$execution$workers == 2L)
  for (field in c("depth", "junctions", "rna_bases")) all_equal_table(full[[field]], parallel_result[[field]], names(full[[field]]))
  # The final validation may fail even after all samples were processed. It
  # must never publish COMPLETE for a source-changing run.
  ctx$mutate_vcf <- TRUE
  all_expect_error(analyze_variant(resources, variant, cfg, progress_file = progress), "VCF or its index changed")
  stopifnot(jsonlite::read_json(progress)$state == "FAILED")
  cat("Full-cohort streaming, 39-sample selection, parallel bounds, interruption/checkpoint/source validation: PASS\n")
  invisible(TRUE)
}

run_all_native_tests <- function() {
  executables <- Sys.which(c("bcftools", "samtools", "regtools"))
  if (any(!nzchar(executables))) stop("--native requires bcftools, samtools and RegTools.")
  work <- tempfile("regshiny_all_native_"); dir.create(work, mode = "0700")
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  run <- function(tool, args) {
    p <- processx::run(unname(executables[[tool]]), as.character(args), error_on_status = FALSE, timeout = 120)
    if (p$status != 0) stop(p$stderr)
  }
  sam <- file.path(work, "template.sam"); make_demo_sam(sam)
  bam <- file.path(work, "template.bam")
  run("samtools", c("sort", "-o", bam, sam)); run("samtools", c("index", bam))
  ids <- sprintf("N%02d", 1:36); bams <- file.path(work, paste0(ids, ".bam"))
  for (path in bams) stopifnot(file.copy(bam, path), file.copy(paste0(bam, ".bai"), paste0(path, ".bai")))
  raw <- file.path(work, "synthetic.vcf"); vcf <- paste0(raw, ".gz")
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=chrDemo,length=1000>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", ids), collapse = "\t"),
    paste(c("chrDemo", 190, ".", "A", "C", ".", "PASS", ".", "GT", rep(c("0/0", "0/1", "1/1"), each = 12)), collapse = "\t")), raw)
  run("bcftools", c("view", "-Oz", "-o", vcf, raw)); run("bcftools", c("index", "-t", vcf))
  registry <- file.path(work, "vcfs.tsv"); manifest <- file.path(work, "samples.tsv")
  utils::write.table(data.frame(chrom = "chrDemo", vcf = vcf, build = "synthetic_v1"), registry, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(data.frame(vcf_sample = ids, bam = bams, build = "synthetic_v1"), manifest, sep = "\t", quote = FALSE, row.names = FALSE)
  resources <- load_variant_resources(registry, manifest); variant <- lookup_variants(resources, "chrDemo:190:A:C")
  cfg <- default_demo_config(); cfg$demo <- FALSE
  full <- analyze_variant(resources, variant, cfg, checkpoint_dir = file.path(work, "checkpoints"))
  stopifnot(nrow(full$samples) == 36L, all(full$samples$status == "SUCCESS"),
            all(full$group_summary$selected_n == 12L), all(full$group_summary$analyzed_n == 12L),
            all(full$samples$native_audit_passed), all(full$samples$rna_valid_bases == 3L))
  exact <- full$junctions[full$junctions$intron_start1 == 201L, ]
  stopifnot(all(exact$count_sum == 24L), all(exact$count_mean == 2L))
  preview <- analyze_variant(resources, variant, cfg, sampling_mode = "preview", max_per_group = 3L)
  all_equal_table(full$depth, preview$depth, c("genotype", "chrom", "pos1", "mean_depth"))
  all_equal_table(full$junctions, preview$junctions, c("genotype", "chrom", "intron_start1", "intron_end1", "strand", "count_mean"))
  resumed <- analyze_variant(resources, variant, cfg, checkpoint_dir = file.path(work, "checkpoints"))
  stopifnot(resumed$provenance$checkpoint$reused == 36L)
  cat("Native full cohort: 36 BAMs, 12 per genotype, junction/depth/base audit and resume: PASS\n")
}

run_all_mock_tests()
if ("--native" %in% commandArgs(trailingOnly = TRUE)) run_all_native_tests()

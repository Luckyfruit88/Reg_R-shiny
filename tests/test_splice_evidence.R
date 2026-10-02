# Rscript tests/test_splice_evidence.R [--native]
# Synthetic reference and saved-result fixtures; no BAM analysis or study data.
# --native exercises real samtools faidx, otherwise only that leaf call is mocked.
source("splice_evidence_backend.R")

splice_expect_error <- function(expr, text = NULL) {
  e <- tryCatch({ force(expr); NULL }, error = function(e) e)
  stopifnot(inherits(e, "error"))
  if (!is.null(text)) stopifnot(grepl(text, conditionMessage(e), fixed = TRUE))
}

run_splice_tests <- function(native = FALSE) {
  work <- tempfile("regshiny_splice_test_"); dir.create(work, mode = "0700")
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  sequence <- paste(rep("C", 220L), collapse = "")
  substr(sequence, 28L, 29L) <- "NG"
  substr(sequence, 60L, 61L) <- "GT"
  substr(sequence, 99L, 100L) <- "AG"
  substr(sequence, 107L, 108L) <- "AG"
  substr(sequence, 149L, 150L) <- "AC"
  substr(sequence, 169L, 170L) <- "CT"
  substr(sequence, 199L, 200L) <- "AA"
  fasta <- file.path(work, "synthetic.fa")
  writeLines(c(">chr1", sequence), fasta)
  writeLines("chr1\t220\t6\t220\t221", paste0(fasta, ".fai"))
  fixture_mtime <- as.POSIXct(floor(as.numeric(Sys.time())), origin = "1970-01-01", tz = "UTC")
  Sys.setFileTime(fasta, fixture_mtime); Sys.setFileTime(paste0(fasta, ".fai"), fixture_mtime)
  original_fetch <- splice_run_faidx
  on.exit(assign("splice_run_faidx", original_fetch, .GlobalEnv), add = TRUE)
  if (!native) assign("splice_run_faidx", function(fasta, region) {
    # Assert the query receives only a private alias/private index.
    stopifnot(nzchar(Sys.readlink(fasta)), file.exists(paste0(fasta, ".fai")))
    pos <- as.integer(strsplit(sub("^chr1:", "", region), "-", fixed = TRUE)[[1L]])
    paste0(">", region, "\n", substr(sequence, pos[[1L]], pos[[2L]]), "\n")
  }, .GlobalEnv)
  if (native && !nzchar(Sys.which("samtools"))) stop("--native requires samtools.")
  reference <- load_splice_reference(fasta, receipt = "")
  before <- splice_file_state(c(fasta, paste0(fasta, ".fai")))
  source_hash <- splice_hash_file(fasta); fai_hash <- splice_hash_file(paste0(fasta, ".fai"))
  receipt <- file.path(work, "receipt.json")
  jsonlite::write_json(list(fasta_path = fasta, size_bytes = before$size_bytes[[1L]],
    mtime_epoch = before$mtime_epoch[[1L]], sha256 = source_hash$value, fai_sha256 = fai_hash$value),
    receipt, auto_unbox = TRUE, digits = NA)
  with_receipt <- load_splice_reference(fasta, receipt)
  stopifnot(with_receipt$prior_receipt$sha256 == source_hash$value)
  fixture <- function(position = 100L, ref = "G", alt = "A", start = 70L, end = 125L) {
    groups <- data.frame(genotype = c("0/0", "0/1", "1/1", "NO_CALL:./."),
      call_status = c(rep("CALLED", 3L), "NO_CALL"), dna_n = c(2L, 2L, 2L, 1L),
      available_n = c(2L, 2L, 2L, 1L), selected_n = c(2L, 2L, 2L, 0L),
      analyzed_n = c(2L, 2L, 2L, 0L), failed_n = 0L)
    depth <- do.call(rbind, lapply(seq_len(nrow(groups)), function(i)
      data.frame(genotype = groups$genotype[[i]], chrom = "chr1", pos1 = seq.int(start, end),
        mean_depth = if (i == 4L) NA_real_ else if (i == 3L) 0 else 12 / i, analyzed_n = groups$analyzed_n[[i]])))
    junctions <- do.call(rbind, lapply(seq_len(nrow(groups)), function(i)
      data.frame(genotype = groups$genotype[[i]], chrom = "chr1", intron_start1 = 60L,
        intron_end1 = c(100L, 108L), strand = "+", count_sum = if (i == 4L) NA_real_ else
          list(c(20, 0), c(10, 16), c(0, 32))[[i]],
        count_mean = if (i == 4L) NA_real_ else list(c(10, 0), c(5, 8), c(0, 16))[[i]],
        analyzed_n = groups$analyzed_n[[i]])))
    list(variant = data.frame(chrom = "chr1", pos1 = position, ref = ref, alt = alt, build = "GRCh38", contig_length = 220),
      group_summary = groups, depth = depth, junctions = junctions,
      provenance = list(config = list(chrom = "chr1", start1 = start, end1 = end), selection = list(mode = "all")))
  }
  annotation <- list(
    introns = data.frame(chrom = "chr1", intron_start1 = c(60L, 60L), intron_end1 = c(100L, 108L),
      strand = "+", transcript_id = c("TX1.2", "TX2.3")),
    transcripts = data.frame(chrom = "chr1", start1 = 1L, end1 = 220L, strand = "+",
      transcript_id = c("TX1.2", "TX2.3")),
    metadata = list(build = "GRCh38", release = "48",
      query = list(chrom = "chr1", start1 = 1L, end1 = 220L, full_transcript_models = TRUE, truncated = FALSE)))
  result <- fixture(); original_result <- serialize(result, NULL)
  model <- build_splice_evidence(result, annotation, reference)
  id <- "chr1:99:100:+:acceptor"
  target <- model$sites[model$sites$candidate_id == id, ]
  stopifnot(model$metadata$reference_status == "REF_MATCH", model$metadata$focus_candidate_id == id,
    nrow(target) == 1L, target$ref_motif == "AG", target$alt_motif == "AA",
    target$change == "DISRUPTED_CANONICAL", target$variant_overlaps, !target$sequence_only,
    identical(model$genotypes$motif_summary[1:3], c("AG/AG", "AG/AA", "AA/AA")),
    identical(model$genotypes$dna_alleles[1:3], c("G/G", "G/A", "A/A")),
    all(is.na(model$coverage$mean_depth[model$coverage$genotype == "NO_CALL:./."])),
    all(model$coverage$coverage_status[model$coverage$genotype == "1/1"] == "OBSERVED_ZERO"),
    identical(serialize(result, NULL), original_result))
  nearby <- model$sites[model$sites$candidate_id == "chr1:107:108:+:acceptor", ]
  stopifnot(nrow(nearby) == 1L, nearby$ref_motif == "AG", nearby$alt_motif == "AG",
    nearby$change == "NOT_OVERLAPPING", !nearby$variant_overlaps,
    all(model$junction_support$count_sum == result$junctions$count_sum, na.rm = TRUE),
    all(model$junction_support$count_mean == result$junctions$count_mean, na.rm = TRUE),
    nrow(model$junction_support) == nrow(result$junctions))
  tight <- build_splice_evidence(result, annotation, reference, flank = 10L)
  wide <- build_splice_evidence(result, annotation, reference, flank = 100L)
  stopifnot(tight$window$start1 == 90L, wide$window$start1 == 70L, wide$window$end1 == 125L,
    wide$window$clipped, !identical(tight$metadata$sequence_query$query_hash, wide$metadata$sequence_query$query_hash))
  missing_depth <- result; missing_depth$depth <- missing_depth$depth[
    !(missing_depth$depth$genotype == "0/1" & missing_depth$depth$pos1 == 105L), ]
  gap <- build_splice_evidence(missing_depth, annotation, reference)
  stopifnot(gap$coverage$coverage_status[gap$coverage$genotype == "0/1" & gap$coverage$pos1 == 105L] == "UNMEASURED")
  mismatch <- fixture(ref = "T")
  bad <- build_splice_evidence(mismatch, annotation, reference)
  stopifnot(bad$metadata$reference_status == "REF_MISMATCH", all(is.na(bad$genotypes$motif_summary)),
    all(startsWith(bad$sites$change, "UNAVAILABLE")), all(is.na(bad$bases$alt_base)),
    all(!bad$sequences$available), all(is.na(bad$sequences$sequence)), bad$bases$ref_base[bad$bases$is_variant] == "G",
    identical(bad$coverage$mean_depth, model$coverage$mean_depth))
  indel <- build_splice_evidence(fixture(alt = "GA"), annotation, reference)
  stopifnot(indel$metadata$reference_status == "UNAVAILABLE_NON_SNV_OR_AMBIGUOUS_ALLELE",
    all(is.na(indel$genotype_sites$motif_alleles)))
  absent <- build_splice_evidence(result, NULL, NULL)
  stopifnot(absent$metadata$reference_status == "UNAVAILABLE_REFERENCE",
    all(is.na(absent$bases$ref_base)), identical(absent$coverage$mean_depth, model$coverage$mean_depth))
  only_rna <- build_splice_evidence(result, NULL, reference)
  stopifnot(only_rna$sites$ref_motif[only_rna$sites$candidate_id == id] == "AG",
    all(!grepl("SEQUENCE_ONLY", only_rna$sites$source, fixed = TRUE)))
  unknown <- result; unknown$junctions$strand <- "?"
  unknown_model <- build_splice_evidence(unknown, NULL, reference)
  stopifnot(all(unknown_model$sites$change == "UNAVAILABLE_UNKNOWN_STRAND"),
    all(is.na(unknown_model$genotype_sites$motif_alleles)), is.na(unknown_model$metadata$focus_candidate_id))
  away <- build_splice_evidence(fixture(position = 104L, ref = "C", alt = "T"), annotation, reference)
  stopifnot(all(away$sites$change[!away$sites$sequence_only] == "NOT_OVERLAPPING"),
    is.na(away$metadata$focus_candidate_id), all(is.na(away$genotypes$motif_summary)))
  # Negative-strand donor AC -> AT is GT -> AT after reverse complement.
  reverse <- fixture(150L, "C", "T", 125L, 175L)
  reverse$junctions$intron_start1 <- 110L; reverse$junctions$intron_end1 <- 150L; reverse$junctions$strand <- "-"
  neg <- build_splice_evidence(reverse, NULL, reference)
  donor <- neg$sites[neg$sites$site_type == "donor", ]
  stopifnot(donor$ref_motif == "GT", donor$alt_motif == "AT", donor$change == "DISRUPTED_CANONICAL",
    identical(neg$genotypes$motif_summary[1:3], c("GT/GT", "GT/AT", "AT/AT")))
  reverse <- fixture(169L, "C", "T", 145L, 195L)
  reverse$junctions$intron_start1 <- 169L; reverse$junctions$intron_end1 <- 190L; reverse$junctions$strand <- "-"
  neg_acceptor <- build_splice_evidence(reverse, NULL, reference)
  acceptor <- neg_acceptor$sites[neg_acceptor$sites$site_type == "acceptor", ]
  stopifnot(acceptor$ref_motif == "AG", acceptor$alt_motif == "AA", acceptor$change == "DISRUPTED_CANONICAL")
  # Creation requires a supported transcript strand and remains sequence-only.
  creation <- fixture(200L, "A", "G", 180L, 220L); creation$junctions <- creation$junctions[FALSE, ]
  made <- build_splice_evidence(creation, annotation, reference)
  gains <- made$sites[made$sites$change == "CREATED_CANONICAL", ]
  stopifnot(nrow(gains) >= 1L, all(gains$sequence_only), all(gains$junction_ids == ""),
    any(gains$ref_motif == "AA" & gains$alt_motif == "AG"),
    nrow(build_splice_evidence(creation, NULL, reference)$sites) == 0L)
  creation_multi <- creation; creation_multi$variant$alt <- "G,T"
  creation_multi$group_summary$genotype[[3L]] <- "2/2"
  creation_multi$depth$genotype[creation_multi$depth$genotype == "1/1"] <- "2/2"
  created_multi <- build_splice_evidence(creation_multi, annotation, reference)
  candidate <- created_multi$sites[created_multi$sites$candidate_id == "chr1:199:200:+:acceptor", ]
  stopifnot(nrow(candidate) == 2L, setequal(candidate$alt_motif, c("AG", "AT")),
    created_multi$genotype_sites$motif_alleles[created_multi$genotype_sites$genotype == "2/2" &
      created_multi$genotype_sites$candidate_id == "chr1:199:200:+:acceptor"] == "AT/AT")
  multi <- fixture(alt = "A,T")
  multi$group_summary$genotype[[3L]] <- "1/2"
  multi$depth$genotype[multi$depth$genotype == "1/1"] <- "1/2"
  multi$junctions$genotype[multi$junctions$genotype == "1/1"] <- "1/2"
  mult <- build_splice_evidence(multi, annotation, reference)
  stopifnot(nrow(mult$sites[mult$sites$candidate_id == id, ]) == 2L,
    mult$genotypes$motif_summary[[3L]] == "AA/AT", mult$genotypes$dna_alleles[[3L]] == "A/T")
  haploid <- result
  haploid$group_summary$genotype[[3L]] <- "1"
  haploid$depth$genotype[haploid$depth$genotype == "1/1"] <- "1"
  haploid$junctions$genotype[haploid$junctions$genotype == "1/1"] <- "1"
  one_allele <- build_splice_evidence(haploid, annotation, reference)
  stopifnot(one_allele$genotypes$ploidy[[3L]] == 1L, one_allele$genotypes$dna_alleles[[3L]] == "A",
    one_allele$genotypes$motif_summary[[3L]] == "AA")
  ambiguous <- fixture(29L, "G", "A", 10L, 50L)
  ambiguous$junctions$intron_start1 <- 10L; ambiguous$junctions$intron_end1 <- 29L
  uncertain <- build_splice_evidence(ambiguous, NULL, reference)
  stopifnot(uncertain$metadata$reference_status == "REF_MATCH",
    all(uncertain$sites$change[uncertain$sites$site_type == "acceptor"] == "UNAVAILABLE_AMBIGUOUS_REFERENCE"),
    all(is.na(uncertain$genotypes$motif_summary)), is.na(uncertain$metadata$focus_candidate_id))
  wrong_build <- result; wrong_build$variant$build <- "GRCh37"
  blocked <- build_splice_evidence(wrong_build, annotation, reference)
  stopifnot(blocked$metadata$reference_status == "REFERENCE_BUILD_MISMATCH", !blocked$metadata$annotation_available)
  wrong_length <- reference; wrong_length$contigs[["chr1"]] <- 219
  blocked <- build_splice_evidence(result, annotation, wrong_length)
  stopifnot(blocked$metadata$reference_status == "REFERENCE_LENGTH_MISMATCH")
  splice_expect_error(build_splice_evidence(result, annotation, reference, 101L), "integer")
  outside <- result; outside$provenance$config$end1 <- 99L
  splice_expect_error(build_splice_evidence(outside, annotation, reference), "outside")
  old <- file.info(paste0(fasta, ".fai"))$mtime
  Sys.setFileTime(paste0(fasta, ".fai"), old + 10)
  changed <- build_splice_evidence(result, annotation, reference)
  stopifnot(changed$metadata$reference_status == "UNAVAILABLE_REFERENCE")
  Sys.setFileTime(paste0(fasta, ".fai"), old)
  stopifnot(identical(before, splice_file_state(c(fasta, paste0(fasta, ".fai")))),
    identical(source_hash, splice_hash_file(fasta)), identical(fai_hash, splice_hash_file(paste0(fasta, ".fai"))))
  cat(if (native) "Native faidx + " else "Mock faidx + ",
    "saved-result splice motifs, strand, nearby junctions, multiallelic, missingness, ROI and provenance: PASS\n", sep = "")
}

run_splice_tests(native = "--native" %in% commandArgs(trailingOnly = TRUE))

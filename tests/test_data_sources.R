# Rscript tests/test_data_sources.R [--native]
# Source/profile tests use only synthetic IDs and reference headers.
# --native uses real samtools/bcftools; otherwise leaf header commands are mocked.
source("data_sources.R")

source_expect_error <- function(expr, text = NULL) {
  error <- tryCatch({ force(expr); NULL }, error = function(e) e)
  stopifnot(inherits(error, "error"))
  if (!is.null(text)) stopifnot(grepl(text, conditionMessage(error), fixed = TRUE))
}

run_data_source_tests <- function(native = FALSE) {
  app <- normalizePath(getwd())
  work <- tempfile("regshiny-source-test-"); dir.create(work, mode = "0700")
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  state <- file.path(work, "state"); dir.create(state, mode = "0700")
  env_names <- c("PATH", "REGSHINY_VCF_MANIFEST", "REGSHINY_SAMPLE_MANIFEST", "REGSHINY_JOB_ROOT",
                 "REGSHINY_ANNOTATION_DB", "REGSHINY_REFERENCE_FASTA", "REGSHINY_REFERENCE_RECEIPT")
  environment_before <- Sys.getenv(env_names, unset = NA_character_)
  run_original <- data_source_run
  on.exit(assign("data_source_run", run_original, .GlobalEnv), add = TRUE)
  fakebin <- file.path(work, "bin"); dir.create(fakebin)
  child_env <- character()
  if (!native) {
    python <- data_source_python()
    fake <- paste0("#!", python, "\nimport sys\nfrom pathlib import Path\n",
      "if Path(sys.argv[0]).name == 'samtools' and '-c' in sys.argv:\n",
      " p=Path(sys.argv[-2]+'.bai')\n",
      " if p.read_bytes()!=b'synthetic index\\n': sys.exit('invalid BAM index')\n",
      " print(0)\n",
      "elif Path(sys.argv[0]).name == 'bcftools' and sys.argv[1]=='query':\n",
      " p=Path(sys.argv[-1]+'.tbi')\n",
      " if p.read_bytes()!=b'synthetic index\\n': sys.exit('invalid VCF index')\n",
      "else: sys.stdout.write(Path(sys.argv[-1]).read_text())\n")
    for (tool in c("bcftools", "samtools")) {
      writeLines(fake, file.path(fakebin, tool))
      Sys.chmod(file.path(fakebin, tool), "0700")
    }
    child_env <- c(PATH = paste(fakebin, Sys.getenv("PATH"), sep = .Platform$path.sep))
    assign("data_source_run", function(script, args, log_path, timeout = 3600) {
      output <- paste0(log_path, ".out")
      p <- processx::run(data_source_python(), c(script, as.character(args)), env = child_env,
        stdout = output, stderr = log_path, timeout = timeout, error_on_status = FALSE)
      if (p$status != 0L) stop(paste(readLines(log_path, warn = FALSE), collapse = "\n"))
      invisible(output)
    }, .GlobalEnv)
  } else stopifnot(all(nzchar(Sys.which(c("samtools", "bcftools")))))
  native_run <- function(tool, args) {
    p <- processx::run(unname(Sys.which(tool)), args, timeout = 120, error_on_status = FALSE)
    if (p$status != 0L) stop(p$stderr)
  }
  make_bam <- function(path, sample = "DNA_A", extra_rg = NULL, lengths = c(chr1 = 1000L)) {
    lines <- c("@HD\tVN:1.6\tSO:coordinate",
      paste0("@SQ\tSN:", names(lengths), "\tLN:", lengths),
      paste0("@RG\tID:rg1\tSM:", sample), extra_rg)
    if (native) {
      sam <- paste0(path, ".sam"); writeLines(lines, sam)
      native_run("samtools", c("view", "-b", "-o", path, sam))
      native_run("samtools", c("index", path)); unlink(sam)
    } else {
      writeLines(lines, path); writeLines("synthetic index", paste0(path, ".bai"))
    }
    normalizePath(path)
  }
  make_vcf <- function(path, samples = c("DNA_A", "DNA_B"), lengths = c(chr1 = 1000L)) {
    header <- c("##fileformat=VCFv4.2",
      paste0("##contig=<ID=", names(lengths), ",length=", lengths, ">"),
      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
      paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples), collapse = "\t"))
    if (native) {
      raw <- sub("\\.gz$", "", path); writeLines(header, raw)
      native_run("bcftools", c("view", "-Oz", "-o", path, raw))
      native_run("bcftools", c("index", "-t", path)); unlink(raw)
    } else {
      writeLines(header, path); writeLines("synthetic index", paste0(path, ".tbi"))
    }
    normalizePath(path)
  }
  a_dir <- file.path(work, "person-a"); b_dir <- file.path(work, "person-b")
  dir.create(a_dir); dir.create(b_dir)
  bam_a <- make_bam(file.path(a_dir, "Aligned.bam"), "RNA_A")
  bam_b <- make_bam(file.path(b_dir, "Aligned.bam"), "RNA_B")
  vcf <- make_vcf(file.path(work, "synthetic.vcf.gz"))
  mapping <- file.path(work, "mapping.tsv")
  write_mapping <- function(samples, paths, path = mapping) utils::write.table(
    data.frame(vcf_sample = samples, bam = paths), path, sep = "\t", quote = FALSE, row.names = FALSE)
  write_mapping(c("DNA_A", "DNA_B"), c(bam_a, bam_b))
  spec <- custom_source_spec(label = "Synthetic explicit mapping", vcf = vcf, mapping = mapping,
    annotation_gtf = "", reference_fasta = "")
  progress <- file.path(work, "progress.json")
  bundle <- prepare_data_profile(spec, app, state, progress_file = progress)
  stopifnot(inherits(bundle$variant_resources, "regshiny_variant_resources"),
    nrow(bundle$variant_resources$samples) == 2L, length(bundle$bam_choices) == 2L,
    !anyDuplicated(names(bundle$bam_choices)), all(grepl("DNA_", names(bundle$bam_choices), fixed = TRUE)),
    grepl("/profiles/custom-", bundle$profile_dir, fixed = TRUE),
    dirname(bundle$job_root) == bundle$profile_dir,
    jsonlite::read_json(progress)$state == "READY",
    inherits(bundle$annotation_resources, "error"), inherits(bundle$splice_reference, "error"))
  for (path in c(bundle$profile_dir, bundle$job_root))
    stopifnot(bitwAnd(as.integer(file.info(path)$mode), 63L) == 0L)
  before <- unname(tools::md5sum(bundle$config$sample_manifest))
  again <- prepare_data_profile(spec, app, state)
  stopifnot(identical(bundle$id, again$id), identical(before, unname(tools::md5sum(bundle$config$sample_manifest))))
  listed <- list_data_profiles(app, state)
  stopifnot(nrow(listed) == 1L, listed$id == bundle$id, listed$profile_dir == bundle$profile_dir)
  # Distinct profiles have distinct job histories and cannot be loaded as another clone.
  other_spec <- custom_source_spec(label = "Second dataset", bam = bam_a, annotation_gtf = "", reference_fasta = "")
  other <- prepare_data_profile(other_spec, app, state)
  stopifnot(other$id != bundle$id, other$job_root != bundle$job_root,
            inherits(other$variant_resources, "error"), length(other$bam_choices) == 1L)
  another_clone <- file.path(work, "another-clone"); dir.create(another_clone)
  source_expect_error(load_data_profile(bundle$profile_dir, another_clone, state), "clone/account")
  source_expect_error(load_data_profile(file.path(bundle$profile_dir, "..", basename(bundle$profile_dir)), app, state), "Select a saved profile")
  # Advisory locks reject a concurrent holder and recover after its process exits.
  lock_path <- file.path(work, "lock")
  lock <- data_source_lock(lock_path, app)
  source_expect_error(data_source_lock(lock_path, app), "already being prepared")
  lock$kill()
  unlocked <- data_source_lock(lock_path, app); unlocked$kill()
  symlink_state <- file.path(work, "state-link"); file.symlink(state, symlink_state)
  source_expect_error(data_source_state_root(app, symlink_state), "symbolic link")
  ignore_repo <- file.path(work, "ignore-fixture-repository")
  processx::run("git", c("init", "--quiet", ignore_repo))
  protected <- file.path(ignore_repo, "nondefault-state")
  dir.create(protected, mode = "0700")
  on.exit(unlink(protected, recursive = TRUE), add = TRUE)
  protected_uid <- data_source_state_root(app, protected)
  private_json <- file.path(protected_uid, "synthetic-private-audit.json")
  writeLines('{"synthetic":true}', private_json)
  ignored <- processx::run("git", c("-C", ignore_repo, "check-ignore", "-v", private_json), error_on_status = FALSE)
  stopifnot(ignored$status == 0L, grepl(".gitignore:1:*", ignored$stdout, fixed = TRUE))
  conflicting <- file.path(work, "existing-ignore"); dir.create(conflicting, mode = "0700")
  writeLines("*.tmp", file.path(conflicting, ".gitignore"))
  source_expect_error(data_source_state_root(app, conflicting), "was not changed")
  stopifnot(identical(readLines(file.path(conflicting, ".gitignore")), "*.tmp"))
  # Explicit mappings are authoritative; incompatible RG SM never overrides them.
  stopifnot(setequal(bundle$variant_resources$samples$vcf_sample, c("DNA_A", "DNA_B")))
  rg_bam <- make_bam(file.path(work, "chosen-name-unrelated-to-ID.bam"), "DNA_A")
  rg_spec <- custom_source_spec(vcf = vcf, bam = rg_bam, mapping_mode = "read_group",
    annotation_gtf = "", reference_fasta = "")
  rg <- prepare_data_profile(rg_spec, app, state)
  stopifnot(rg$variant_resources$samples$vcf_sample == "DNA_A")
  wrong_bam <- make_bam(file.path(work, "DNA_A.bam"), "UNMATCHED_RG")
  wrong_spec <- custom_source_spec(vcf = vcf, bam = wrong_bam, mapping_mode = "read_group",
    annotation_gtf = "", reference_fasta = "")
  source_expect_error(prepare_data_profile(wrong_spec, app, state, progress_file = progress), "does not exactly match")
  stopifnot(jsonlite::read_json(progress)$state == "FAILED")
  mixed <- make_bam(file.path(work, "mixed.bam"), "DNA_A", "@RG\tID:rg2\tSM:DNA_B")
  mixed_spec <- custom_source_spec(vcf = vcf, bam = mixed, mapping_mode = "read_group",
    annotation_gtf = "", reference_fasta = "")
  source_expect_error(prepare_data_profile(mixed_spec, app, state), "mixed")
  # Duplicate sample mappings and symlink/hardlink aliases are never collapsed.
  bad_map <- file.path(work, "duplicate-map.tsv")
  write_mapping(c("DNA_A", "DNA_A"), c(bam_a, bam_b), bad_map)
  bad_spec <- spec; bad_spec$mapping <- bad_map
  source_expect_error(prepare_data_profile(bad_spec, app, state), "one-to-one")
  alias <- file.path(work, "alias.bam"); file.symlink(bam_a, alias)
  write_mapping(c("DNA_A", "DNA_B"), c(bam_a, alias), bad_map)
  source_expect_error(prepare_data_profile(bad_spec, app, state), "duplicate physical")
  unlink(alias); file.link(bam_a, alias); file.copy(paste0(bam_a, ".bai"), paste0(alias, ".bai"))
  write_mapping(c("DNA_A", "DNA_B"), c(bam_a, alias), bad_map)
  source_expect_error(prepare_data_profile(bad_spec, app, state), "hard-link")
  # Incompatible VCF/BAM assembly lengths fail before a bundle becomes active.
  mismatch <- make_bam(file.path(work, "mismatch.bam"), "DNA_A", lengths = c(chr1 = 1001L))
  wrong_spec$bam <- mismatch
  source_expect_error(prepare_data_profile(wrong_spec, app, state), "lengths differ")
  source_expect_error(prepare_data_profile(custom_source_spec(vcf = vcf, bam = bam_a,
    annotation_gtf = "", reference_fasta = ""), app, state), "requires an explicit mapping")
  bad_index_bam <- make_bam(file.path(work, "bad-index.bam"), "DNA_A")
  writeBin(as.raw(1:9), paste0(bad_index_bam, ".bai"))
  index_spec <- custom_source_spec(vcf = vcf, bam = bad_index_bam, mapping_mode = "read_group",
    annotation_gtf = "", reference_fasta = "")
  source_expect_error(prepare_data_profile(index_spec, app, state), "index")
  bad_index_vcf <- make_vcf(file.path(work, "bad-index.vcf.gz"))
  writeBin(as.raw(1:9), paste0(bad_index_vcf, ".tbi"))
  index_spec$vcf <- bad_index_vcf; index_spec$bam <- rg_bam
  source_expect_error(prepare_data_profile(index_spec, app, state), "index")
  # Changed source identity and changed private manifests are rejected on reload.
  original_manifest <- readBin(rg$config$sample_manifest, "raw", n = file.info(rg$config$sample_manifest)$size)
  cat("\n", file = rg$config$sample_manifest, append = TRUE)
  source_expect_error(load_data_profile(rg$profile_dir, app, state), "manifest changed")
  writeBin(original_manifest, rg$config$sample_manifest)
  old_time <- file.info(bam_a)$mtime; Sys.setFileTime(bam_a, old_time + 10)
  source_expect_error(load_data_profile(bundle$profile_dir, app, state), "source changed")
  Sys.setFileTime(bam_a, old_time)
  # FHS setup uses original source maps, without any prior Reg_Shiny deployment.
  fhs <- file.path(work, "fhs"); dir.create(fhs)
  fhs_bams <- file.path(fhs, "bams"); fhs_vcfs <- file.path(fhs, "vcfs")
  dir.create(fhs_bams); dir.create(fhs_vcfs)
  lengths <- c(chr1 = 248956422L, chr17 = 83257441L, chr21 = 46709983L, chrX = 156040895L)
  make_bam(file.path(fhs_bams, "RNA_001.accepted_hits.merged.markeddups.recal.bam"), "RNA_001", lengths = lengths)
  for (chrom in c(paste0("chr", 1:22), "chrX")) {
    size <- if (chrom %in% names(lengths)) lengths[[chrom]] else 1000L
    make_vcf(file.path(fhs_vcfs, paste0(chrom, ".vcf.gz")), samples = c("DNA_A", "DNA_B"), lengths = setNames(size, chrom))
  }
  rna_map <- file.path(fhs, "original-rna.csv"); bam_map <- file.path(fhs, "original-bams.csv")
  utils::write.csv(data.frame(NWGC_ID = c("RNA_001", "RNA_002"), framid = c("DNA_A", "DNA_B"), Batch_ID = "SYNTHETIC"), rna_map, row.names = FALSE)
  utils::write.csv(data.frame(NWGC_ID = c("RNA_001", "RNA_002"), directory = fhs_bams), bam_map, row.names = FALSE)
  fhs_spec <- default_fhs_spec(rna_map = rna_map, bam_map = bam_map,
    vcf_dir = fhs_vcfs, annotation_gtf = "", reference_fasta = "")
  fhs_bundle <- prepare_data_profile(fhs_spec, app, state)
  stopifnot(nrow(fhs_bundle$variant_resources$vcf_registry) == 23L,
    nrow(fhs_bundle$variant_resources$samples) == 1L, fhs_bundle$variant_resources$samples$vcf_sample == "DNA_A",
    !grepl("Jinjie_Jay/Reg_Shiny", fhs_bundle$config$sample_manifest, fixed = TRUE))
  added_bam <- make_bam(file.path(fhs_bams, "RNA_002.accepted_hits.merged.markeddups.recal.bam"), "RNA_002", lengths = lengths)
  source_expect_error(load_data_profile(fhs_bundle$profile_dir, app, state), "source changed")
  expanded <- prepare_data_profile(fhs_spec, app, state)
  stopifnot(expanded$id != fhs_bundle$id, nrow(expanded$variant_resources$samples) == 2L,
    expanded$job_root != fhs_bundle$job_root)
  Sys.setFileTime(paste0(added_bam, ".bai"), file.info(paste0(added_bam, ".bai"))$mtime + 5)
  source_expect_error(load_data_profile(expanded$profile_dir, app, state), "source changed")
  stopifnot(identical(environment_before, Sys.getenv(env_names, unset = NA_character_)))
  cat(if (native) "Native" else "Mock-header",
    " source profiles: explicit/RG mapping, FHS raw sources, identity, permissions, locks and isolated history: PASS\n")
}

run_data_source_tests(native = "--native" %in% commandArgs(trailingOnly = TRUE))

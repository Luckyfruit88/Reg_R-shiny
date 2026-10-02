# Private per-account, per-clone source profiles. No session-global environment
# variables are changed. Input data/indexes remain read-only and outside Git.
data_source_fields <- function() list(kind = "custom", label = "Custom data", build = "GRCh38",
  rna_map = "", bam_map = "", vcf_dir = "", vcf = "", vcf_registry = "", mapping = "",
  bam_dir = "", bam = "", mapping_mode = "explicit", annotation_gtf = "", annotation_db = "",
  reference_fasta = "", reference_receipt = "")

default_fhs_spec <- function(
    label = "FHS WGS + RNA-seq", build = "GRCh38",
    rna_map = "/restricted/projectnb/mtdna-alcohol/Jian_Yang/DATA/FHS_data/RNA_seq/IDs/RNA_seq_ID_summary_4batchs.csv",
    bam_map = "/restricted/projectnb/mtdna-alcohol/Jian_Yang/DATA/FHS_data/RNA_seq/IDs/NWGC_bam_files_location.csv",
    vcf_dir = "/restricted/projectnb/sequencing/topMed/data/freeze_10a/passgt.minDP10.fhsids",
    annotation_gtf = "/restricted/projectnb/mtdna-alcohol/Jian_Yang/DATA/Gene/gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz",
    annotation_db = "", reference_fasta = "/share/pkg.8/spliceai/1.3.1/install/share/examples/hg38.fa",
    reference_receipt = "") {
  spec <- data_source_fields()
  values <- as.list(environment()); values$spec <- NULL
  for (name in intersect(names(values), names(spec))) spec[[name]] <- values[[name]]
  spec$kind <- "fhs"
  spec
}

custom_source_spec <- function(label = "Custom data", build = "GRCh38",
    vcf = "", vcf_registry = "", mapping = "", bam_dir = "", bam = "",
    annotation_db = "", annotation_gtf = default_fhs_spec()$annotation_gtf,
    reference_fasta = default_fhs_spec()$reference_fasta, reference_receipt = "",
    mapping_mode = c("explicit", "read_group")) {
  mapping_mode <- match.arg(mapping_mode)
  spec <- data_source_fields(); values <- as.list(environment()); values$spec <- NULL
  for (name in intersect(names(values), names(spec))) spec[[name]] <- values[[name]]
  spec$kind <- "custom"
  spec
}

startup_profile_spec <- function(app_dir = getwd(), config_file = NULL) {
  if (is.null(config_file) || !nzchar(config_file)) return(default_fhs_spec())
  path <- path.expand(config_file)
  if (!startsWith(path, "/")) path <- file.path(app_dir, path)
  if (!file.exists(path)) stop("The explicitly selected profile JSON configuration does not exist.")
  spec <- jsonlite::read_json(path, simplifyVector = TRUE)
  data_source_normalize_spec(spec, app_dir)
}

data_source_uid <- function() {
  value <- trimws(processx::run("id", "-u", timeout = 5)$stdout)
  if (!grepl("^[0-9]+$", value)) stop("Cannot identify the current Unix account.")
  value
}

data_source_private_dir <- function(path, create = FALSE) {
  is_link <- function(path) { value <- Sys.readlink(path); length(value) == 1L && !is.na(value) && nzchar(value) }
  if (is_link(path)) stop("Private profile directories cannot be symbolic links.")
  if (create && !dir.exists(path) && !dir.create(path, recursive = FALSE, mode = "0700"))
    stop("Cannot create a private profile directory.")
  if (!dir.exists(path) || is_link(path)) stop("Private profile directory is unavailable.")
  info <- file.info(path)
  if (is.null(info$uid) || as.character(info$uid) != data_source_uid()) stop("Profile directory belongs to another account.")
  if (bitwAnd(as.integer(info$mode), 63L) != 0L) stop("Profile directory must be private (0700).")
  normalizePath(path, mustWork = TRUE)
}

data_source_state_root <- function(app_dir = getwd(), state_root = NULL) {
  app_dir <- normalizePath(app_dir, mustWork = TRUE)
  base <- if (is.null(state_root)) file.path(app_dir, ".reg_shiny") else path.expand(state_root)
  if (!startsWith(base, "/")) stop("State root must be an absolute persistent path.")
  link <- Sys.readlink(base)
  if (length(link) == 1L && !is.na(link) && nzchar(link)) stop("State root cannot be a symbolic link.")
  if (!dir.exists(base) && !dir.create(base, recursive = TRUE, mode = "0700"))
    stop("The clone/state location is not writable; choose a writable persistent state_root.")
  base <- normalizePath(base, mustWork = TRUE)
  mode <- as.integer(file.info(base)$mode)
  if (bitwAnd(mode, 18L) != 0L && bitwAnd(mode, 512L) == 0L)
    stop("A group/world-writable state root must have the sticky bit; choose a private persistent directory instead.")
  # Protect generated JSON/audits even when state_root is a nondefault directory
  # inside a clone. Never replace or append to an existing user .gitignore.
  ignore <- file.path(base, ".gitignore")
  if (!file.exists(ignore)) {
    temp <- tempfile(".ignore-", tmpdir = base)
    writeLines("*", temp); Sys.chmod(temp, "0644")
    file.link(temp, ignore)
    unlink(temp)
  }
  link <- Sys.readlink(ignore)
  if ((length(link) == 1L && !is.na(link) && nzchar(link)) || !file.exists(ignore) || file.access(ignore, 4) != 0L)
    stop("State-root .gitignore must be a readable regular file.")
  rules <- readLines(ignore, warn = FALSE)
  rules <- rules[!grepl("^[[:space:]]*(#|$)", rules)]
  if (!length(rules) || !identical(tail(rules, 1L), "*"))
    stop("Existing state-root .gitignore must end with an unconditional '*' rule. It was not changed; choose another private state_root.")
  uid <- data_source_private_dir(file.path(base, paste0("uid-", data_source_uid())), create = TRUE)
  data_source_private_dir(file.path(uid, "profiles"), create = TRUE)
  data_source_private_dir(file.path(uid, "resources"), create = TRUE)
  uid
}

data_source_normalize_spec <- function(spec, app_dir) {
  if (!is.list(spec) || is.null(names(spec)) || anyDuplicated(names(spec)) ||
      any(!names(spec) %in% names(data_source_fields()))) stop("Invalid or unknown source-profile configuration fields.")
  base <- if (isTRUE(spec$kind == "fhs")) default_fhs_spec() else custom_source_spec()
  for (name in names(spec)) base[[name]] <- spec[[name]]
  if (any(vapply(base, function(x) !is.character(x) || length(x) != 1L || is.na(x) || grepl("[\r\n\t]", x), logical(1))))
    stop("Profile values must be single text values without line breaks or tabs.")
  if (!base$kind %in% c("fhs", "custom") || !base$mapping_mode %in% c("explicit", "read_group"))
    stop("Unknown profile type or sample-mapping mode.")
  if (!nzchar(base$label) || nchar(base$label) > 120L || !grepl("^[A-Za-z0-9_.-]+$", base$build) ||
      tolower(base$build) %in% c("unknown", "na", "none")) stop("Provide a profile label and an explicit reference assembly.")
  if (base$kind == "fhs" && base$build != "GRCh38") stop("The default FHS data require GRCh38.")
  path_fields <- setdiff(names(base), c("kind", "label", "build", "mapping_mode"))
  for (name in path_fields) if (nzchar(base[[name]])) {
    path <- path.expand(base[[name]])
    if (!startsWith(path, "/")) path <- file.path(app_dir, path)
    base[[name]] <- normalizePath(path, mustWork = FALSE)
  }
  if (base$kind == "custom" && nzchar(base$vcf) && nzchar(base$vcf_registry))
    stop("Choose one indexed VCF/BCF file or a VCF registry.")
  if (nzchar(base$bam) && nzchar(base$bam_dir)) stop("Choose one BAM file or a BAM directory.")
  base
}

data_source_hash <- function(value) {
  path <- tempfile("profile-hash-")
  on.exit(unlink(path), add = TRUE)
  saveRDS(value, path, version = 2L, compress = FALSE)
  unname(tools::md5sum(path))
}

data_source_paths <- function(spec) {
  paths <- unique(unlist(spec[c("rna_map", "bam_map", "vcf", "vcf_registry", "mapping",
    "annotation_gtf", "annotation_db", "reference_fasta", "reference_receipt")], use.names = FALSE))
  paths <- paths[nzchar(paths)]
  if (spec$kind == "fhs") {
    paths <- c(paths, file.path(spec$vcf_dir, paste0(c(paste0("chr", 1:22), "chrX"), ".vcf.gz")))
  } else if (nzchar(spec$vcf_registry) && file.exists(spec$vcf_registry)) {
    x <- utils::read.delim(spec$vcf_registry, sep = "\t", quote = "", comment.char = "", colClasses = "character")
    if (!all(c("chrom", "vcf", "build") %in% names(x))) stop("VCF registry needs chrom, vcf and build columns.")
    paths <- c(paths, x$vcf)
  }
  vcfs <- paths[grepl("\\.(vcf\\.(gz|bgz)|bcf)$", paths, ignore.case = TRUE)]
  if (length(vcfs)) paths <- c(paths, paste0(vcfs, ".csi"), paste0(vcfs, ".tbi"))
  bams <- if (nzchar(spec$bam)) spec$bam else character()
  if (spec$kind == "fhs" && file.exists(spec$rna_map) && file.exists(spec$bam_map)) {
    rna <- utils::read.csv(spec$rna_map, colClasses = "character", check.names = FALSE, na.strings = character())
    locations <- utils::read.csv(spec$bam_map, colClasses = "character", check.names = FALSE, na.strings = character())
    if (!all(c("NWGC_ID", "framid", "Batch_ID") %in% names(rna)) ||
        !all(c("NWGC_ID", "directory") %in% names(locations))) stop("Original FHS source maps have an unexpected schema.")
    locations <- locations[locations$NWGC_ID %in% rna$NWGC_ID, , drop = FALSE]
    # Include every source-map candidate, even currently missing/inaccessible
    # BAMs and indexes. Adding one can change the admitted sample universe.
    bams <- c(bams, file.path(locations$directory,
      paste0(locations$NWGC_ID, ".accepted_hits.merged.markeddups.recal.bam")))
  }
  if (nzchar(spec$bam_dir) && dir.exists(spec$bam_dir))
    bams <- c(bams, list.files(spec$bam_dir, pattern = "\\.bam$", ignore.case = TRUE, full.names = TRUE))
  if (nzchar(spec$mapping) && file.exists(spec$mapping)) {
    x <- utils::read.delim(spec$mapping, sep = "\t", quote = "", comment.char = "", colClasses = "character")
    if (!all(c("vcf_sample", "bam") %in% names(x))) stop("Sample mapping needs vcf_sample and bam columns.")
    bams <- c(bams, x$bam)
  }
  bams <- unique(bams)
  if (length(bams)) paths <- c(paths, bams, paste0(bams, ".bai"), paste0(bams, ".csi"),
    sub("\\.bam$", ".bai", bams, ignore.case = TRUE), sub("\\.bam$", ".csi", bams, ignore.case = TRUE))
  if (nzchar(spec$reference_fasta)) paths <- c(paths, paste0(spec$reference_fasta, ".fai"))
  sort(unique(paths[nzchar(paths)]))
}

data_source_identity <- function(spec) {
  paths <- data_source_paths(spec)
  if (any(!startsWith(paths, "/"))) stop("Registry/mapping source paths must be absolute.")
  info <- file.info(paths)
  # Small authoritative maps, annotation GTF and FAI are content bound. Large
  # BAM/VCF/FASTA files use size/mtime plus native runtime identity checks.
  hash_paths <- paths %in% unlist(spec[c("rna_map", "bam_map", "vcf_registry", "mapping",
    "annotation_gtf", "reference_receipt")]) | endsWith(paths, ".fai")
  digest <- rep(NA_character_, length(paths))
  ok <- hash_paths & file.exists(paths) & !info$isdir & file.access(paths, 4) == 0
  digest[ok] <- unname(tools::md5sum(paths[ok]))
  data.frame(path = paths, resolved_path = vapply(paths, normalizePath, character(1), mustWork = FALSE),
    exists = file.exists(paths), readable = file.access(paths, 4) == 0L, size_bytes = info$size,
    mtime_epoch = as.numeric(info$mtime), ctime_epoch = as.numeric(info$ctime), md5 = digest,
    stringsAsFactors = FALSE, row.names = NULL)
}

data_source_atomic_rds <- function(value, path) {
  temp <- tempfile(".profile-", tmpdir = dirname(path)); on.exit(unlink(temp), add = TRUE)
  saveRDS(value, temp, version = 2L); Sys.chmod(temp, "0600")
  if (!file.rename(temp, path)) stop("Could not publish private profile state.")
}

data_source_progress <- function(path, state, stage, message, id = NULL) {
  if (is.null(path)) return(invisible(NULL))
  if (!dir.exists(dirname(path))) stop("Preparation progress directory does not exist.")
  temp <- tempfile(".progress-", tmpdir = dirname(path)); on.exit(unlink(temp), add = TRUE)
  jsonlite::write_json(list(state = state, stage = stage, message = message, id = id,
    updated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)), temp, auto_unbox = TRUE, na = "null")
  Sys.chmod(temp, "0600")
  if (!file.rename(temp, path)) stop("Could not publish preparation progress.")
}

data_source_python <- function() {
  path <- unname(Sys.which("python3"))
  if (!nzchar(path)) stop("Python 3 is required for profile preparation.")
  path
}

data_source_lock <- function(path, app_dir) {
  process <- processx::process$new(data_source_python(),
    c(file.path(app_dir, "scripts", "profile_lock.py"), "--path", path),
    stdin = "|", stdout = "|", stderr = "|", cleanup = TRUE)
  for (i in 1:50) {
    process$poll_io(100)
    lines <- process$read_output_lines()
    if ("LOCKED" %in% lines) return(process)
    if ("BUSY" %in% lines || !process$is_alive()) break
  }
  process$kill()
  stop("This profile/resource is already being prepared in another session, or its lock is unavailable. Retry after that preparation finishes; another data profile can still be selected.")
}

data_source_run <- function(script, args, log_path, timeout = 3600) {
  output <- paste0(log_path, ".out")
  p <- processx::run(data_source_python(), c(script, as.character(args)),
    stdout = output, stderr = log_path, timeout = timeout, error_on_status = FALSE)
  for (path in c(log_path, output)) if (file.exists(path)) Sys.chmod(path, "0600")
  if (p$status != 0L || isTRUE(p$timeout)) {
    message <- if (file.exists(log_path)) paste(tail(readLines(log_path, warn = FALSE), 8L), collapse = "\n") else ""
    stop("Source preparation failed: ", message, call. = FALSE)
  }
  invisible(output)
}

data_source_api <- function(app_dir) {
  e <- new.env(parent = globalenv())
  for (name in c("backend.R", "variant_backend.R", "annotation_backend.R", "splice_evidence_backend.R"))
    sys.source(file.path(app_dir, name), envir = e)
  e$.regshiny_variant_source <- file.path(app_dir, "variant_backend.R")
  e$.regshiny_annotation_source <- file.path(app_dir, "annotation_backend.R")
  e
}

data_source_profile_path <- function(profile_dir, app_dir, state_root = NULL) {
  root <- data_source_state_root(app_dir, state_root)
  if (length(profile_dir) != 1L || !startsWith(profile_dir, "/") ||
      dirname(profile_dir) != file.path(root, "profiles") || !grepl("^(fhs|custom)-[0-9a-f]{20}$", basename(profile_dir)))
    stop("Select a saved profile from this clone and the current account.")
  data_source_private_dir(profile_dir)
}

load_data_profile <- function(profile_dir, app_dir = getwd(), state_root = NULL) {
  app_dir <- normalizePath(app_dir, mustWork = TRUE)
  profile_dir <- data_source_profile_path(profile_dir, app_dir, state_root)
  path <- file.path(profile_dir, "profile.rds")
  config <- readRDS(path)
  if (!is.list(config) || !identical(config$schema, "regshiny-source-profile-v1") ||
      !identical(config$id, basename(profile_dir)) || !identical(config$owner_uid, data_source_uid()) ||
      !identical(config$app_dir, app_dir)) stop("Saved profile identity does not match this clone/account.")
  if (!identical(config$source_identity, data_source_identity(config$spec)))
    stop("A configured source changed. Prepare a new profile before loading these data; old job history remains in its original profile.")
  if (!identical(unname(tools::md5sum(names(config$file_md5))), unname(config$file_md5)))
    stop("A private source manifest changed or is missing; prepare a fresh profile.")
  e <- data_source_api(app_dir)
  variant_resources <- if (isTRUE(config$has_vcf))
    e$load_variant_resources(config$vcf_manifest, config$sample_manifest) else simpleError("This profile has BAM data only; no DNA VCF/mapping was configured.")
  if (!inherits(variant_resources, "error")) {
    vcfs <- unique(variant_resources$vcf_registry$vcf)
    if (any(!file.exists(vcfs)) || any(file.access(vcfs, 4) != 0L)) stop("A configured DNA VCF/BCF is no longer readable.")
    indexed <- vapply(vcfs, function(vcf) {
      paths <- c(paste0(vcf, ".csi"), paste0(vcf, ".tbi"))
      any(file.exists(paths) & file.access(paths, 4) == 0L)
    }, logical(1))
    if (!all(indexed)) stop("A configured DNA index is no longer readable.")
  }
  bams <- utils::read.delim(config$bam_choices, sep = "\t", quote = "", colClasses = "character")$bam
  readable <- file.exists(bams) & file.access(bams, 4) == 0L
  if (!all(readable)) stop("A selected BAM is no longer readable; restore source access or prepare a new profile.")
  indexed <- vapply(bams, function(bam) {
    paths <- c(paste0(bam, ".bai"), paste0(bam, ".csi"),
      sub("\\.bam$", ".bai", bam, ignore.case = TRUE), sub("\\.bam$", ".csi", bam, ignore.case = TRUE))
    any(file.exists(paths) & file.access(paths, 4) == 0L)
  }, logical(1))
  if (!all(indexed)) stop("A selected BAM index is no longer available.")
  issues <- config$issues
  annotation_resources <- if (config$spec$build != "GRCh38") simpleError("GENCODE v48 requires GRCh38.") else
    tryCatch(e$load_annotation_resources(config$annotation_db), error = function(error) error)
  splice_reference <- if (config$spec$build != "GRCh38") simpleError("The configured splice reference overlay requires GRCh38.") else
    tryCatch(e$load_splice_reference(config$spec$reference_fasta, receipt = config$spec$reference_receipt), error = function(error) error)
  if (!inherits(splice_reference, "error") && !inherits(variant_resources, "error")) {
    registry <- variant_resources$vcf_registry
    if ("length_bp" %in% names(registry)) {
      expected <- unname(splice_reference$contigs[registry$chrom])
      actual <- suppressWarnings(as.numeric(registry$length_bp))
      if (any(!is.na(expected) & !is.na(actual) & expected != actual))
        stop("Configured VCF and reference FASTA chromosome lengths differ; the data profile is incompatible.")
    }
  }
  for (resource in list(annotation_resources, splice_reference))
    if (inherits(resource, "error")) issues <- c(issues, conditionMessage(resource))
  bam_roots <- unique(dirname(bams))
  labels <- basename(bams)
  duplicate_labels <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  if (any(duplicate_labels)) {
    mapped <- if (!inherits(variant_resources, "error"))
      variant_resources$samples$vcf_sample[match(bams, variant_resources$samples$bam)] else rep(NA_character_, length(bams))
    labels[duplicate_labels] <- ifelse(!is.na(mapped[duplicate_labels]),
      paste(mapped[duplicate_labels], labels[duplicate_labels], sep = " | "), bams[duplicate_labels])
  }
  list(id = config$id, label = config$spec$label, variant_resources = variant_resources,
    bam_root = if (length(bam_roots) == 1L) bam_roots[[1L]] else NULL,
    bam_choices = stats::setNames(bams, labels),
    annotation_resources = annotation_resources, splice_reference = splice_reference,
    job_root = data_source_private_dir(file.path(profile_dir, "jobs")), profile_dir = profile_dir,
    config = config, issues = unique(issues))
}

list_data_profiles <- function(app_dir = getwd(), state_root = NULL) {
  root <- data_source_state_root(app_dir, state_root)
  paths <- list.dirs(file.path(root, "profiles"), recursive = FALSE, full.names = TRUE)
  answer <- lapply(paths, function(path) tryCatch({
    data_source_profile_path(path, normalizePath(app_dir), state_root)
    config <- readRDS(file.path(path, "profile.rds"))
    if (!identical(config$schema, "regshiny-source-profile-v1") || !identical(config$app_dir, normalizePath(app_dir)) ||
        !identical(config$owner_uid, data_source_uid())) return(NULL)
    data.frame(profile_dir = path, id = config$id, label = config$spec$label,
      kind = config$spec$kind, created_at = config$created_at, stringsAsFactors = FALSE)
  }, error = function(e) NULL))
  answer <- Filter(Negate(is.null), answer)
  if (length(answer)) do.call(rbind, answer) else data.frame(profile_dir = character(), id = character(),
    label = character(), kind = character(), created_at = character(), stringsAsFactors = FALSE)
}

prepare_data_profile <- function(spec, app_dir = getwd(), state_root = NULL, progress_file = NULL,
                                 build_annotation = TRUE) {
  app_dir <- normalizePath(app_dir, mustWork = TRUE)
  root <- data_source_state_root(app_dir, state_root)
  done <- FALSE; stage <- "VALIDATING"; id <- NULL
  publish <- function(message) data_source_progress(progress_file, "PREPARING", stage, message, id)
  on.exit(if (!done) try(data_source_progress(progress_file, "FAILED", stage,
    "Preparation did not complete. Existing active data and saved jobs are unchanged; correct the inputs or select another profile.", id), silent = TRUE), add = TRUE)
  publish("Validating source configuration and private state location.")
  spec <- data_source_normalize_spec(spec, app_dir)
  identity <- data_source_identity(spec)
  id <- paste0(spec$kind, "-", substr(data_source_hash(list(schema = "regshiny-source-profile-v1",
    app_dir = app_dir, spec = spec, identity = identity)), 1L, 20L))
  profile_dir <- data_source_private_dir(file.path(root, "profiles", id), create = TRUE)
  if (is.null(progress_file)) progress_file <- file.path(profile_dir, "progress.json")
  lock <- data_source_lock(file.path(profile_dir, ".prepare.lock"), app_dir)
  on.exit(try(lock$kill(), silent = TRUE), add = TRUE)
  if (file.exists(file.path(profile_dir, "profile.rds"))) {
    bundle <- load_data_profile(profile_dir, app_dir, state_root)
    repair_annotation <- isTRUE(build_annotation) && inherits(bundle$annotation_resources, "error") &&
      spec$build == "GRCh38" && !nzchar(spec$annotation_db) &&
      nzchar(spec$annotation_gtf) && file.exists(spec$annotation_gtf)
    if (!repair_annotation) {
      data_source_progress(progress_file, "READY", "READY", "Validated cached private data profile.", id)
      done <- TRUE
      return(bundle)
    }
    old_db <- bundle$config$annotation_db
    if (nzchar(old_db) && file.exists(old_db) && startsWith(old_db, paste0(root, "/resources/")))
      file.rename(old_db, tempfile("unreadable-annotation-", tmpdir = dirname(old_db), fileext = ".sqlite"))
  }
  manifest_dir <- file.path(profile_dir, "manifests")
  if (dir.exists(manifest_dir)) {
    commit <- tryCatch(readRDS(file.path(manifest_dir, "COMMIT.rds")), error = function(e) NULL)
    valid <- !is.null(commit) && identical(commit$id, id) && identical(commit$source_identity, identity) &&
      identical(unname(tools::md5sum(file.path(manifest_dir, names(commit$file_md5)))), unname(commit$file_md5))
    if (!valid && !file.rename(manifest_dir, tempfile("rejected-manifests-", tmpdir = profile_dir)))
      stop("An incomplete private manifest attempt cannot be safely archived.")
  }
  if (!dir.exists(manifest_dir)) {
    stage <- "PREPARING_DATA"
    publish(if (spec$kind == "fhs") "Preparing FHS manifests from the original RNA-ID/BAM-location maps and 23 indexed WGS files." else
      "Checking selected BAM/VCF headers, existing indexes and exact sample mapping.")
    attempt <- tempfile("import-", tmpdir = profile_dir); dir.create(attempt, mode = "0700")
    if (spec$kind == "fhs") {
      required <- c(spec$rna_map, spec$bam_map)
      if (any(!file.exists(required)) || any(file.access(required, 4) != 0L) || !dir.exists(spec$vcf_dir))
        stop("Default FHS sources are not accessible under this account. Select custom data or restore SCC project permissions.")
      data_source_run(file.path(app_dir, "scripts", "build_manifests.py"),
        c("--rna-map", spec$rna_map, "--bam-map", spec$bam_map, "--vcf-dir", spec$vcf_dir,
          "--output-dir", attempt, "--build", spec$build), file.path(profile_dir, "prepare-data.stderr"))
      samples <- utils::read.delim(file.path(attempt, "sample_manifest.tsv"), sep = "\t", quote = "", colClasses = "character")
      utils::write.table(data.frame(bam = samples$bam), file.path(attempt, "bam_choices.tsv"),
                         sep = "\t", quote = FALSE, row.names = FALSE)
      Sys.chmod(file.path(attempt, "bam_choices.tsv"), "0600")
    } else {
      args <- c("--build", spec$build, "--output-dir", attempt, "--mapping-mode", spec$mapping_mode)
      flags <- c(vcf = "--vcf", vcf_registry = "--vcf-registry", mapping = "--mapping", bam = "--bam", bam_dir = "--bam-dir")
      for (field in names(flags)) if (nzchar(spec[[field]])) args <- c(args, flags[[field]], spec[[field]])
      data_source_run(file.path(app_dir, "scripts", "import_data_sources.py"), args,
        file.path(profile_dir, "prepare-data.stderr"))
    }
    if (!identical(identity, data_source_identity(spec))) stop("A source changed during preparation; retry against a stable source snapshot.")
    prepared_files <- list.files(attempt, pattern = "\\.(tsv|json)$", full.names = TRUE)
    hashes <- tools::md5sum(prepared_files); names(hashes) <- basename(prepared_files)
    data_source_atomic_rds(list(id = id, source_identity = identity, file_md5 = hashes), file.path(attempt, "COMMIT.rds"))
    if (!file.rename(attempt, manifest_dir)) stop("Could not atomically publish private source manifests.")
  }
  annotation_db <- spec$annotation_db; issues <- character()
  if (!nzchar(annotation_db) && spec$build == "GRCh38" && nzchar(spec$annotation_gtf)) {
    annotation_db <- tryCatch({
      if (!file.exists(spec$annotation_gtf) || file.access(spec$annotation_gtf, 4) != 0L)
        stop("GENCODE source GTF is unavailable. RNA data remain usable.")
      key <- data_source_hash(list(gtf = spec$annotation_gtf, gtf_md5 = unname(tools::md5sum(spec$annotation_gtf)),
        fai = if (nzchar(spec$reference_fasta) && file.exists(paste0(spec$reference_fasta, ".fai")))
          unname(tools::md5sum(paste0(spec$reference_fasta, ".fai"))) else NULL,
        builder = unname(tools::md5sum(file.path(app_dir, "scripts", "build_annotation.py")))))
      resource_dir <- data_source_private_dir(file.path(root, "resources", paste0("gencode48-", key)), create = TRUE)
      db <- file.path(resource_dir, "annotation.sqlite")
      if (!file.exists(db)) {
        if (!isTRUE(build_annotation)) stop("GENCODE index has not been prepared for this clone.")
        stage <- "PREPARING_ANNOTATION"; publish("Building this account's private GENCODE v48 index; original GTF and RNA data are read-only.")
        annotation_lock <- data_source_lock(file.path(resource_dir, ".build.lock"), app_dir)
        on.exit(try(annotation_lock$kill(), silent = TRUE), add = TRUE)
        if (!file.exists(db)) {
          args <- c("build", "--gtf", spec$annotation_gtf, "--output", db, "--release", "48", "--build", "GRCh38",
            "--source-url", "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz")
          fai <- paste0(spec$reference_fasta, ".fai")
          if (nzchar(spec$reference_fasta) && file.exists(fai)) args <- c(args, "--fai", fai)
          data_source_run(file.path(app_dir, "scripts", "build_annotation.py"), args,
            file.path(resource_dir, "build.stderr"))
          Sys.chmod(db, "0600")
        }
        annotation_lock$kill()
      }
      db
    }, error = function(e) { issues <<- c(issues, conditionMessage(e)); "" })
  }
  data_source_private_dir(file.path(profile_dir, "jobs"), create = TRUE)
  files <- list.files(manifest_dir, pattern = "\\.(tsv|json|rds)$", full.names = TRUE)
  if (!length(files)) stop("Prepared manifests are missing.")
  config <- list(schema = "regshiny-source-profile-v1", id = id, app_dir = app_dir,
    owner_uid = data_source_uid(), spec = spec, source_identity = identity,
    has_vcf = file.exists(file.path(manifest_dir, "vcf_manifest.tsv")),
    vcf_manifest = file.path(manifest_dir, "vcf_manifest.tsv"),
    sample_manifest = file.path(manifest_dir, "sample_manifest.tsv"),
    bam_choices = file.path(manifest_dir, "bam_choices.tsv"), annotation_db = annotation_db,
    file_md5 = tools::md5sum(files), issues = issues,
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
  if (!identical(identity, data_source_identity(spec))) stop("A source changed before profile publication.")
  data_source_atomic_rds(config, file.path(profile_dir, "profile.rds"))
  stage <- "LOADING"; publish("Validating the prepared manifests and loading optional reference overlays.")
  bundle <- tryCatch(load_data_profile(profile_dir, app_dir, state_root), error = function(e) {
    file.rename(file.path(profile_dir, "profile.rds"), tempfile("failed-profile-", tmpdir = profile_dir, fileext = ".rds"))
    stop(e)
  })
  data_source_progress(progress_file, "READY", "READY", "Private data profile is ready.", id)
  done <- TRUE
  bundle
}

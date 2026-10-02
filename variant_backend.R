# Whole-genome extension for Reg_Shiny. Source backend.R first.
# DNA GT comes exclusively from the selected VCF record. RNA observations are
# descriptive alignment evidence, never substituted for a missing DNA genotype.
.regshiny_variant_source <- tryCatch(normalizePath(sys.frame(1L)$ofile, mustWork = TRUE),
                                    error = function(e) normalizePath("variant_backend.R", mustWork = FALSE))

variant_lines <- function(x) {
  if (!nzchar(x)) character() else strsplit(sub("\n$", "", x), "\n", fixed = TRUE)[[1L]]
}

variant_file_state <- function(paths) {
  s <- file.info(paths)
  data.frame(path = paths, size = s$size, mtime = as.numeric(s$mtime),
             stringsAsFactors = FALSE, row.names = NULL)
}

variant_run <- function(args, timeout = 60) {
  tool <- unname(Sys.which("bcftools"))
  if (!nzchar(tool)) stop("bcftools is not available in PATH.", call. = FALSE)
  p <- processx::run(tool, as.character(args), timeout = timeout,
                    error_on_status = FALSE, cleanup_tree = TRUE)
  if (isTRUE(p$timeout) || is.na(p$status) || p$status != 0L)
    stop("bcftools failed: ", p$stderr, call. = FALSE)
  if (nchar(p$stdout, type = "bytes") > 16 * 1024^2)
    stop("The locus response exceeds 16 MiB; this record cannot be loaded interactively.")
  list(text = p$stdout, log = paste(c("bcftools", vapply(as.character(args), shQuote, character(1))), collapse = " "),
       stderr = p$stderr)
}

variant_read_manifest <- function(path, required, label) {
  if (length(path) != 1L || is.na(path) || !nzchar(path) || !file.exists(path))
    stop(label, " is not configured or does not exist.", call. = FALSE)
  path <- normalizePath(path, mustWork = TRUE)
  if (file.info(path)$size > 16 * 1024^2) stop(label, " exceeds 16 MiB.")
  x <- utils::read.delim(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
                        colClasses = "character", check.names = FALSE, na.strings = character())
  if (anyDuplicated(names(x)) || !all(required %in% names(x)) || !nrow(x))
    stop(label, " requires nonempty columns: ", paste(required, collapse = ", "))
  if (any(vapply(x[required], function(z) anyNA(z) || any(!nzchar(z)) || any(grepl("[\r\n]", z)), logical(1))))
    stop(label, " has empty or invalid required values.")
  list(data = x, path = path, md5 = unname(tools::md5sum(path)))
}

load_variant_resources <- function(vcf_manifest = Sys.getenv("REGSHINY_VCF_MANIFEST"),
                                   sample_manifest = Sys.getenv("REGSHINY_SAMPLE_MANIFEST")) {
  v <- variant_read_manifest(vcf_manifest, c("chrom", "vcf", "build"), "REGSHINY_VCF_MANIFEST")
  s <- variant_read_manifest(sample_manifest, c("vcf_sample", "bam", "build"), "REGSHINY_SAMPLE_MANIFEST")
  vr <- v$data; sm <- s$data
  if (any(!grepl("^[A-Za-z0-9_][A-Za-z0-9_.-]*$", vr$chrom)) || anyDuplicated(vr$chrom))
    stop("The VCF manifest requires one unique, literal chromosome name per row.")
  builds <- unique(c(vr$build, sm$build))
  if (length(builds) != 1L || !grepl("^[A-Za-z0-9_.-]+$", builds) ||
      tolower(builds) %in% c("unknown", "na", "none"))
    stop("VCF and RNA manifests must declare the same known reference assembly; no liftover is performed.")
  if (anyDuplicated(sm$vcf_sample) || any(grepl("[\t\r\n ]", sm$vcf_sample)))
    stop("Sample IDs must be unique, exact VCF sample names without whitespace.")
  if (any(!startsWith(vr$vcf, "/")) || any(!startsWith(sm$bam, "/")))
    stop("Manifest data paths must be absolute, server-configured paths.")
  # Canonicalize existing paths to catch multiple rows referencing the same BAM.
  sm$bam <- vapply(sm$bam, normalizePath, character(1), mustWork = FALSE)
  vr$vcf <- vapply(vr$vcf, normalizePath, character(1), mustWork = FALSE)
  if (anyDuplicated(sm$bam)) stop("A BAM cannot be assigned to multiple VCF samples.")
  if (any(!grepl("\\.bam$", sm$bam, ignore.case = TRUE))) stop("RNA manifest paths must name BAM files.")
  if (any(!grepl("\\.(vcf\\.(gz|bgz)|bcf)$", vr$vcf, ignore.case = TRUE)))
    stop("VCF manifest paths must name indexed bgzip VCF or BCF files.")
  provenance <- list(vcf_manifest = v$path, sample_manifest = s$path,
    vcf_manifest_md5 = v$md5, sample_manifest_md5 = s$md5,
    build = builds, build_validation = "Explicit matching manifest assembly; contig lengths checked when declared in VCF. Reference sequence identity is not independently proven.",
    loaded_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
  structure(list(vcf_registry = vr, samples = sm, build = builds, provenance = provenance),
            class = "regshiny_variant_resources")
}

variant_check_resources <- function(resources) {
  if (!inherits(resources, "regshiny_variant_resources")) stop("Load trusted variant resources first.")
  p <- resources$provenance
  if (!identical(unname(tools::md5sum(p$vcf_manifest)), p$vcf_manifest_md5) ||
      !identical(unname(tools::md5sum(p$sample_manifest)), p$sample_manifest_md5))
    stop("A resource manifest changed; reload the application before analysis.")
  invisible(TRUE)
}

parse_variant_query <- function(query) {
  if (length(query) != 1L || is.na(query) || nchar(query) > 2000L)
    stop("Enter CHROM:POS, optionally CHROM:POS:REF:ALT.")
  fields <- strsplit(trimws(query), ":", fixed = TRUE)[[1L]]
  if (!length(fields) %in% c(2L, 4L) ||
      !grepl("^[A-Za-z0-9_][A-Za-z0-9_.-]*$", fields[[1L]]) ||
      !grepl("^[0-9]+$", fields[[2L]]))
    stop("Enter a literal chromosome and 1-based position: CHROM:POS or CHROM:POS:REF:ALT.")
  pos <- int_scalar(as.numeric(fields[[2L]]), "Variant position", 1)
  ref <- alt <- NULL
  if (length(fields) == 4L) {
    ref <- toupper(fields[[3L]]); alt <- toupper(fields[[4L]])
    if (!grepl("^[ACGTN]+$", ref) ||
        !grepl("^([ACGTN]+|\\*|<[A-Z0-9_.-]+>)(,([ACGTN]+|\\*|<[A-Z0-9_.-]+>))*$", alt))
      stop("Optional REF and ALT must match the VCF exactly; for breakends use CHROM:POS and select a record.")
  }
  list(chrom = fields[[1L]], pos1 = pos, ref = ref, alt = alt)
}

variant_vcf_info <- function(path, chrom) {
  if (!file.exists(path) || file.access(path, 4) != 0) stop("The configured VCF is not readable: ", path)
  indices <- c(paste0(path, ".csi"), paste0(path, ".tbi"))
  indices <- indices[file.exists(indices)]
  if (!length(indices)) stop("An existing CSI or TBI index is required next to the VCF: ", path)
  state <- variant_file_state(c(path, indices))
  h <- variant_run(c("view", "-h", path))
  header <- variant_lines(h$text)
  contig_line <- header[startsWith(header, paste0("##contig=<ID=", chrom, ",")) |
                          header == paste0("##contig=<ID=", chrom, ">")]
  if (length(contig_line) > 1L) stop("Duplicate selected contig declarations in the VCF header.")
  len <- NA_real_
  if (length(contig_line) && grepl("[,<]length=[0-9]+[,>]", contig_line))
    len <- as.numeric(sub(".*[,<]length=([0-9]+)[,>].*", "\\1", contig_line))
  warnings <- character()
  if (is.na(len)) warnings <- "VCF contig length is undeclared; assembly relies on the explicit manifests and exact contig names."
  if (any(state$mtime[-1L] < state$mtime[[1L]]))
    warnings <- c(warnings, "A VCF index is older than its data file; the administrator must verify that the index matches.")
  list(state = state, contig_length = len, reference_header = header[startsWith(header, "##reference=")],
       warnings = warnings, log = h$log)
}

lookup_variants <- function(resources, query) {
  variant_check_resources(resources)
  q <- parse_variant_query(query)
  reg <- resources$vcf_registry[resources$vcf_registry$chrom == q$chrom, , drop = FALSE]
  if (nrow(reg) != 1L) stop("Chromosome is not configured: ", q$chrom, ". No automatic chr-prefix conversion is performed.")
  info <- variant_vcf_info(reg$vcf, q$chrom)
  if (!is.na(info$contig_length) && q$pos1 > info$contig_length) stop("Position exceeds the declared VCF contig length.")
  region <- paste0(q$chrom, ":", q$pos1, "-", q$pos1)
  p <- variant_run(c("query", "--regions-overlap", "0", "-r", region,
                     "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n", reg$vcf))
  lines <- variant_lines(p$text)
  x <- data.frame(chrom = character(), pos1 = integer(), id = character(), ref = character(),
                  alt = character(), qual = character(), filter = character(), stringsAsFactors = FALSE)
  if (length(lines)) {
    fields <- strsplit(lines, "\t", fixed = TRUE)
    if (any(lengths(fields) != 7L)) stop("Unexpected VCF query output schema.")
    x <- as.data.frame(do.call(rbind, fields), stringsAsFactors = FALSE)
    names(x) <- c("chrom", "pos1", "id", "ref", "alt", "qual", "filter")
    x$pos1 <- as.integer(x$pos1)
    x <- x[x$chrom == q$chrom & x$pos1 == q$pos1, , drop = FALSE]
    if (!is.null(q$ref)) x <- x[x$ref == q$ref & x$alt == q$alt, , drop = FALSE]
  }
  x$build <- rep(resources$build, nrow(x)); x$vcf <- rep(reg$vcf, nrow(x))
  x$contig_length <- rep(info$contig_length, nrow(x))
  x$record_id <- if (nrow(x)) paste(x$build, x$chrom, x$pos1, x$ref, x$alt, sep = ":") else character()
  if (anyDuplicated(x$record_id)) stop("Duplicate VCF records share the same CHROM/POS/REF/ALT. Resolve duplicate records before analysis.")
  if (!identical(info$state, variant_file_state(info$state$path))) stop("The VCF or its index changed during lookup.")
  attr(x, "source_state") <- info$state
  attr(x, "reference_header") <- info$reference_header
  attr(x, "warnings") <- info$warnings
  attr(x, "log") <- c(info$log, p$log, p$stderr)
  x
}

canonical_dna_gt <- function(gt, n_alt) {
  # Phase is preserved in raw_gt; descriptive group labels ignore phase/order.
  if (length(gt) != 1L || is.na(gt) || !grepl("^(\\.|[0-9]+)([/|](\\.|[0-9]+))*$", gt))
    stop("Invalid DNA GT encountered; malformed genotypes cannot be grouped.")
  alleles <- strsplit(gt, "[/|]")[[1L]]
  numeric_alleles <- suppressWarnings(as.integer(alleles[alleles != "."]))
  if (anyNA(numeric_alleles) || any(numeric_alleles > n_alt)) stop("GT allele index exceeds the selected record's ALT alleles.")
  ordered <- c(as.character(sort(numeric_alleles)), rep(".", sum(alleles == ".")))
  status <- if (all(alleles == ".")) "NO_CALL" else if (any(alleles == ".")) "PARTIAL_CALL" else "CALLED"
  label <- paste(ordered, collapse = "/")
  if (status != "CALLED") label <- paste(status, label, sep = ":")
  list(genotype = label, call_status = status, ploidy = length(alleles), phased = grepl("|", gt, fixed = TRUE))
}

variant_get_genotypes <- function(resources, variant) {
  p <- variant_run(c("query", "-l", variant$vcf))
  samples <- variant_lines(p$text)
  if (!length(samples) || anyDuplicated(samples)) stop("The VCF must contain unique sample columns and DNA GT calls.")
  region <- paste0(variant$chrom, ":", variant$pos1, "-", variant$pos1)
  g <- variant_run(c("query", "--regions-overlap", "0", "-r", region,
                     "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", variant$vcf))
  fields <- strsplit(variant_lines(g$text), "\t", fixed = TRUE)
  hits <- Filter(function(f) length(f) >= 4L && identical(f[1:4],
    c(variant$chrom, as.character(variant$pos1), variant$ref, variant$alt)), fields)
  if (length(hits) != 1L || length(hits[[1L]]) != length(samples) + 4L)
    stop("The selected CHROM/POS/REF/ALT does not resolve to one complete genotype record.")
  gt <- hits[[1L]][-(1:4)]
  n_alt <- length(strsplit(variant$alt, ",", fixed = TRUE)[[1L]])
  parsed <- lapply(gt, canonical_dna_gt, n_alt = n_alt)
  x <- data.frame(vcf_sample = samples, raw_gt = gt,
    genotype = vapply(parsed, `[[`, character(1), "genotype"),
    call_status = vapply(parsed, `[[`, character(1), "call_status"),
    ploidy = vapply(parsed, `[[`, integer(1), "ploidy"),
    phased = vapply(parsed, `[[`, logical(1), "phased"), stringsAsFactors = FALSE)
  m <- match(x$vcf_sample, resources$samples$vcf_sample)
  x$linked <- !is.na(m); x$bam <- resources$samples$bam[m]
  attr(x, "log") <- c(p$log, g$log, p$stderr, g$stderr)
  x
}

variant_aggregate <- function(variant, genotypes, samples, sample_results, cfg, rna_supported) {
  groups <- sort(unique(genotypes$genotype))
  summary <- list(); depths <- list(); junctions <- list(); bases <- list()
  success <- samples$status == "SUCCESS"
  global_j <- lapply(sample_results[success], function(z) z$junctions[, c("key", "chrom", "intron_start1", "intron_end1", "strand"), drop = FALSE])
  global_j <- if (length(global_j)) unique(do.call(rbind, global_j)) else data.frame()
  for (i in seq_along(groups)) {
    group <- groups[[i]]
    indices <- which(samples$genotype == group & success)
    ok <- sample_results[indices]; n <- length(ok)
    summary[[i]] <- data.frame(genotype = group,
      call_status = genotypes$call_status[match(group, genotypes$genotype)],
      dna_n = sum(genotypes$genotype == group),
      available_n = sum(genotypes$genotype == group & genotypes$linked),
      linked_n = sum(genotypes$genotype == group & genotypes$linked),
      selected_n = sum(samples$genotype == group), analyzed_n = n,
      failed_n = sum(samples$genotype == group & !success),
      selection_mode = "deterministic_first_n_by_sample_id", stringsAsFactors = FALSE)
    d <- data.frame(genotype = group, chrom = cfg$chrom, pos1 = seq.int(cfg$start1, cfg$end1),
      mean_depth = if (n) Reduce(`+`, lapply(ok, function(z) z$depth$depth)) / n else NA_real_, analyzed_n = n)
    depths[[i]] <- d
    if (nrow(global_j)) {
      j <- global_j
      counts <- lapply(ok, function(z) { m <- match(j$key, z$junctions$key); v <- z$junctions$score[m]; v[is.na(m)] <- 0; v })
      j$genotype <- group; j$count_sum <- if (n) Reduce(`+`, counts) else NA_real_
      j$count_mean <- if (n) j$count_sum / n else NA_real_; j$analyzed_n <- n
      junctions[[i]] <- j
    }
    b <- data.frame(genotype = group, base = c("A", "C", "G", "T", "N"), count = NA_real_,
                    analyzed_n = if (rna_supported) n else 0L,
                    evidence_status = if (!rna_supported) "UNSUPPORTED_NON_SNV" else if (!n) "NO_SUCCESSFUL_RNA_ANALYSIS" else "OBSERVED")
    exclusion_names <- c("skipped_N", "deleted_D", "low_baseq", "missing_sequence", "missing_quality")
    for (name in exclusion_names) b[[paste0("excluded_", name)]] <- NA_real_
    b$valid_base_count <- b$samples_with_valid_bases <- NA_real_
    if (n && rna_supported) b$count <- Reduce(`+`, lapply(ok, function(z) {
      if (is.null(z$rna_bases)) stop("RNA base evidence is missing from a successful analysis.")
      m <- match(b$base, z$rna_bases$base)
      if (anyNA(m)) stop("RNA base evidence has an unexpected schema.")
      as.numeric(z$rna_bases$count[m])
    }))
    if (n && rna_supported) {
      exclusions <- lapply(ok, function(z) {
        e <- attr(z$rna_bases, "excluded")
        if (is.null(e) || !all(exclusion_names %in% names(e)))
          stop("RNA base exclusion audit is missing from a successful analysis.")
        e[exclusion_names]
      })
      totals <- Reduce(`+`, exclusions)
      for (name in exclusion_names) b[[paste0("excluded_", name)]] <- unname(totals[[name]])
      b$valid_base_count <- sum(b$count)
      b$samples_with_valid_bases <- sum(vapply(ok, function(z) sum(z$rna_bases$count) > 0, logical(1)))
      missing <- sum(totals[c("missing_sequence", "missing_quality")]) > 0
      b$evidence_status <- if (sum(b$count) == 0) {
        if (missing) "NO_VALID_BASES_WITH_MISSING_DATA" else "NO_VALID_BASES"
      } else if (missing) "OBSERVED_WITH_MISSING_DATA" else "OBSERVED"
    }
    bases[[i]] <- b
  }
  empty_j <- data.frame(key = character(), chrom = character(), intron_start1 = integer(), intron_end1 = integer(),
    strand = character(), genotype = character(), count_sum = numeric(), count_mean = numeric(), analyzed_n = integer())
  list(group_summary = do.call(rbind, summary), depth = do.call(rbind, depths),
       junctions = if (length(junctions)) do.call(rbind, junctions) else empty_j,
       rna_bases = do.call(rbind, bases))
}

analyze_variant_preview <- function(resources, variant, cfg, max_per_group = 3L) {
  variant_check_resources(resources)
  if (!is.data.frame(variant) || nrow(variant) != 1L ||
      !all(c("record_id", "chrom", "pos1", "ref", "alt", "build", "vcf") %in% names(variant)))
    stop("Select exactly one CHROM/POS/REF/ALT record from the current lookup.")
  max_per_group <- int_scalar(max_per_group, "Maximum RNA samples per genotype", 1, 10)
  cfg <- validate_config(cfg)
  if (cfg$demo || cfg$chrom != variant$chrom || variant$pos1 < cfg$start1 || variant$pos1 > cfg$end1)
    stop("The real RNA interval must contain the selected variant on exactly the same chromosome.")
  current <- lookup_variants(resources, paste0(variant$chrom, ":", variant$pos1))
  current <- current[current$record_id == variant$record_id, , drop = FALSE]
  if (nrow(current) != 1L || !identical(as.character(current$vcf), as.character(variant$vcf)) ||
      !identical(as.character(current$build), as.character(variant$build)))
    stop("The selected variant or resource identity changed; repeat lookup.")
  # Read the current record again, but bind all sample analyses to this snapshot.
  variant <- current
  info <- variant_vcf_info(variant$vcf, variant$chrom)
  calls <- variant_get_genotypes(resources, variant)
  selected <- calls[calls$linked & calls$call_status == "CALLED", , drop = FALSE]
  if (!nrow(selected)) stop("No called DNA genotype has an explicitly linked RNA BAM. Missing and partial DNA calls are not used for RNA genotype comparison.")
  selected <- selected[order(selected$genotype, selected$vcf_sample, method = "radix"), , drop = FALSE]
  selected <- selected[ave(seq_len(nrow(selected)), selected$genotype, FUN = seq_along) <= max_per_group, , drop = FALSE]
  if (nrow(selected) > 30L) stop("The selected genotype groups exceed 30 RNA samples. Reduce samples per genotype.")
  rownames(selected) <- NULL
  selected$status <- "FAILED"; selected$error <- NA_character_
  selected$retained_reads <- selected$junction_count <- selected$mean_depth <- NA_real_
  selected$rna_valid_bases <- NA_real_
  for (name in c("skipped_N", "deleted_D", "low_baseq", "missing_sequence", "missing_quality"))
    selected[[paste0("rna_excluded_", name)]] <- NA_real_
  selected$native_audit_passed <- NA
  alleles <- c(variant$ref, strsplit(variant$alt, ",", fixed = TRUE)[[1L]])
  rna_supported <- all(grepl("^[ACGTN]$", alleles))
  results <- vector("list", nrow(selected)); logs <- attr(calls, "log")
  warnings <- info$warnings
  if (!rna_supported) warnings <- c(warnings, "RNA base comparison is available only for SNV records. Indels and symbolic alleles retain DNA GT, junction and depth evidence; RNA allele support is unavailable.")
  if (any(!calls$linked)) warnings <- c(warnings, paste(sum(!calls$linked), "VCF samples have no configured RNA BAM and are not analyzed."))
  if (any(calls$call_status != "CALLED")) warnings <- c(warnings, "Missing and partial DNA calls are shown in the audit and are excluded from RNA genotype comparison.")
  for (i in seq_len(nrow(selected))) {
    sample_cfg <- cfg; sample_cfg$bam <- selected$bam[[i]]; sample_cfg$demo <- FALSE
    sample_cfg$variant_position <- if (rna_supported) variant$pos1 else NULL
    result <- tryCatch({
      contigs <- read_bam_contigs(sample_cfg$bam)
      m <- match(variant$chrom, contigs$chrom)
      if (is.na(m)) stop("The BAM has no exact matching chromosome name.")
      if (!is.na(info$contig_length) && contigs$length_bp[[m]] != info$contig_length)
        stop("VCF and BAM contig lengths differ; reference assembly compatibility failed.")
      analyze_bam(sample_cfg)
    }, error = function(e) e)
    if (inherits(result, "error")) {
      selected$error[[i]] <- conditionMessage(result)
    } else {
      results[[i]] <- result
      selected$status[[i]] <- "SUCCESS"; selected$retained_reads[[i]] <- nrow(result$reads)
      selected$junction_count[[i]] <- nrow(result$junctions)
      selected$mean_depth[[i]] <- mean(result$depth$depth)
      selected$native_audit_passed[[i]] <- isTRUE(result$native_audit_passed)
      if (rna_supported) {
        selected$rna_valid_bases[[i]] <- sum(result$rna_bases$count)
        e <- attr(result$rna_bases, "excluded")
        for (name in names(e)) selected[[paste0("rna_excluded_", name)]][[i]] <- unname(e[[name]])
      }
      logs <- c(logs, paste("SAMPLE", selected$vcf_sample[[i]]), result$log)
      warnings <- c(warnings, result$warnings)
    }
  }
  if (!identical(info$state, variant_file_state(info$state$path)))
    stop("The VCF or its index changed during analysis; results were not published.")
  variant_check_resources(resources)
  aggregated <- variant_aggregate(variant, calls, selected, results, cfg, rna_supported)
  if (any(aggregated$group_summary$available_n > aggregated$group_summary$selected_n & aggregated$group_summary$call_status == "CALLED"))
    warnings <- c(warnings, "RNA results describe a deterministic bounded subset of linked samples, not the complete genotype cohort or a random sample.")
  if (any(selected$status != "SUCCESS")) warnings <- c(warnings, "One or more selected RNA analyses failed; their evidence remains NA and they are excluded from successful-sample means.")
  version <- variant_run("--version")
  provenance <- c(resources$provenance, list(record_id = variant$record_id,
    vcf_state = info$state, vcf_reference_header = info$reference_header,
    selection = list(mode = "preview", order = "genotype, then exact VCF sample ID (radix sort)", max_per_group = max_per_group,
                     max_total = 30L, replacement_on_failure = FALSE),
    rna_evidence_status = if (rna_supported) "SNV_BASE_COUNTS" else "UNSUPPORTED_NON_SNV",
    count_unit = "retained primary alignment record; paired ends separate",
    aggregate_definition = "Unweighted mean across successful selected BAMs, including observed zeros; failed samples excluded and counted.",
    bcftools_version = version$text, config = cfg,
    completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE), log = c(info$log, logs)))
  c(list(variant = variant, genotypes = calls, samples = selected, sample_results = results,
         provenance = provenance, warnings = unique(warnings)), aggregated)
}

# Full-cohort mode keeps only one BAM's detailed result in memory. Each committed
# checkpoint is an independent contribution, replayed once into a fresh stream
# after restart. No cumulative counter/sum pair can be partially committed.
variant_hash_object <- function(x) {
  path <- tempfile("regshiny_hash_")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path, version = 2L, compress = FALSE); Sys.chmod(path, "0600")
  unname(tools::md5sum(path))
}

variant_atomic_rds <- function(object, path) {
  tmp <- tempfile(".pending_", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp, version = 2L); Sys.chmod(tmp, "0600")
  if (!file.rename(tmp, path)) stop("Could not atomically commit checkpoint: ", path)
  invisible(path)
}

variant_write_progress <- function(progress, path) {
  if (is.null(path)) return(invisible(NULL))
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required for persistent progress.")
  if (length(path) != 1L || !nzchar(path) || !dir.exists(dirname(path)))
    stop("The progress file's parent directory must exist.")
  tmp <- tempfile(".progress_", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  jsonlite::write_json(progress, tmp, auto_unbox = TRUE, pretty = FALSE, na = "null")
  Sys.chmod(tmp, "0600")
  if (!file.rename(tmp, path)) stop("Could not atomically publish progress.")
  invisible(path)
}

variant_bam_state <- function(bam) {
  indices <- unique(c(paste0(bam, ".bai"), sub("\\.bam$", ".bai", bam, ignore.case = TRUE),
                      paste0(bam, ".csi"), sub("\\.bam$", ".csi", bam, ignore.case = TRUE)))
  # Missing paths are retained as NA, so a repaired missing BAM/index invalidates
  # its failed checkpoint while unchanged failures remain explicitly recorded.
  variant_file_state(c(bam, indices))
}

variant_code_identity <- function() {
  paths <- c(.regshiny_variant_source, file.path(dirname(.regshiny_variant_source), "backend.R"))
  if (any(!file.exists(paths))) stop("Cannot bind checkpoint identity to backend source files.")
  list(files = data.frame(file = basename(paths), md5 = unname(tools::md5sum(paths))),
       R = R.version.string,
       processx = as.character(utils::packageVersion("processx")),
       data_table = as.character(utils::packageVersion("data.table")))
}

variant_stream_init <- function(calls, cfg, rna_supported) {
  s <- new.env(parent = emptyenv())
  s$cfg <- cfg; s$rna_supported <- rna_supported
  s$groups <- sort(unique(calls$genotype))
  s$totals <- setNames(lapply(s$groups, function(group) list(
    n = 0L, depth = numeric(cfg$end1 - cfg$start1 + 1L), junctions = numeric(),
    bases = setNames(numeric(5L), c("A", "C", "G", "T", "N")),
    exclusions = setNames(numeric(5L), c("skipped_N", "deleted_D", "low_baseq", "missing_sequence", "missing_quality")),
    samples_with_bases = 0L)), s$groups)
  s$junctions <- data.frame(key = character(), chrom = character(), intron_start1 = integer(),
                           intron_end1 = integer(), strand = character(), stringsAsFactors = FALSE)
  s
}

variant_stream_add <- function(stream, genotype, contribution) {
  if (!identical(contribution$status, "SUCCESS")) return(invisible(NULL))
  x <- stream$totals[[genotype]]
  if (is.null(x) || length(contribution$depth) != length(x$depth) ||
      anyNA(contribution$depth) || any(!is.finite(contribution$depth)) || any(contribution$depth < 0))
    stop("Checkpoint depth or genotype does not match the analysis contract.")
  x$depth <- x$depth + contribution$depth; x$n <- x$n + 1L
  j <- contribution$junctions
  if (nrow(j)) {
    if (anyDuplicated(j$key) || anyNA(j$score) || any(j$score < 0)) stop("Invalid checkpoint junction counts.")
    added <- !j$key %in% names(x$junctions)
    x$junctions <- c(x$junctions, setNames(numeric(sum(added)), j$key[added]))
    x$junctions[j$key] <- x$junctions[j$key] + j$score
    new_global <- !j$key %in% stream$junctions$key
    if (any(new_global)) stream$junctions <- rbind(stream$junctions,
      j[new_global, c("key", "chrom", "intron_start1", "intron_end1", "strand"), drop = FALSE])
  }
  if (stream$rna_supported) {
    b <- contribution$audit$rna_bases; e <- attr(b, "excluded")
    m <- match(names(x$bases), b$base)
    if (anyNA(m) || is.null(e) || !all(names(x$exclusions) %in% names(e)) ||
        anyNA(b$count[m]) || any(b$count[m] < 0)) stop("Invalid checkpoint RNA base evidence.")
    x$bases <- x$bases + b$count[m]
    x$exclusions <- x$exclusions + e[names(x$exclusions)]
    x$samples_with_bases <- x$samples_with_bases + as.integer(sum(b$count) > 0)
  }
  stream$totals[[genotype]] <- x
  invisible(NULL)
}

variant_stream_finish <- function(stream, calls, samples) {
  summary <- depths <- junctions <- bases <- vector("list", length(stream$groups))
  global_j <- stream$junctions
  if (nrow(global_j)) global_j <- global_j[order(global_j$chrom, global_j$intron_start1,
                                                global_j$intron_end1, global_j$strand), , drop = FALSE]
  for (i in seq_along(stream$groups)) {
    group <- stream$groups[[i]]; x <- stream$totals[[group]]; n <- x$n
    summary[[i]] <- data.frame(genotype = group,
      call_status = calls$call_status[match(group, calls$genotype)], dna_n = sum(calls$genotype == group),
      available_n = sum(calls$genotype == group & calls$linked), linked_n = sum(calls$genotype == group & calls$linked),
      selected_n = sum(samples$genotype == group), analyzed_n = n,
      failed_n = sum(samples$genotype == group & samples$status != "SUCCESS"),
      selection_mode = "all_linked_called_samples", stringsAsFactors = FALSE)
    depths[[i]] <- data.frame(genotype = group, chrom = stream$cfg$chrom,
      pos1 = seq.int(stream$cfg$start1, stream$cfg$end1),
      mean_depth = if (n) x$depth / n else NA_real_, analyzed_n = n)
    if (nrow(global_j)) {
      j <- global_j; m <- match(j$key, names(x$junctions)); counts <- x$junctions[m]
      counts[is.na(m)] <- 0
      j$genotype <- group; j$count_sum <- if (n) counts else NA_real_
      j$count_mean <- if (n) counts / n else NA_real_; j$analyzed_n <- n
      junctions[[i]] <- j
    }
    supported <- stream$rna_supported && n > 0L
    missing <- sum(x$exclusions[c("missing_sequence", "missing_quality")]) > 0
    status <- if (!stream$rna_supported) "UNSUPPORTED_NON_SNV" else if (!n) "NO_SUCCESSFUL_RNA_ANALYSIS" else
      if (sum(x$bases) == 0) { if (missing) "NO_VALID_BASES_WITH_MISSING_DATA" else "NO_VALID_BASES" } else
      if (missing) "OBSERVED_WITH_MISSING_DATA" else "OBSERVED"
    b <- data.frame(genotype = group, base = names(x$bases), count = if (supported) unname(x$bases) else NA_real_,
                    analyzed_n = if (stream$rna_supported) n else 0L, evidence_status = status)
    for (name in names(x$exclusions)) b[[paste0("excluded_", name)]] <- if (supported) unname(x$exclusions[[name]]) else NA_real_
    b$valid_base_count <- if (supported) sum(x$bases) else NA_real_
    b$samples_with_valid_bases <- if (supported) x$samples_with_bases else NA_real_
    bases[[i]] <- b
  }
  empty_j <- data.frame(key = character(), chrom = character(), intron_start1 = integer(), intron_end1 = integer(),
    strand = character(), genotype = character(), count_sum = numeric(), count_mean = numeric(), analyzed_n = integer())
  junctions <- Filter(Negate(is.null), junctions)
  list(group_summary = do.call(rbind, summary), depth = do.call(rbind, depths),
       junctions = if (length(junctions)) do.call(rbind, junctions) else empty_j,
       rna_bases = do.call(rbind, bases))
}

variant_compact_contribution <- function(result, source_state) {
  if (!isTRUE(result$native_audit_passed)) stop("Native junction/CIGAR audit did not pass.")
  audit_fields <- c("config", "candidate_alignments", "span_only_excluded", "native_audit_passed",
                    "versions", "source_metadata", "completed_at", "rna_bases", "warnings")
  audit <- result[intersect(audit_fields, names(result))]
  audit$versions$regtools_extract_help <- NULL
  audit$retained_reads <- nrow(result$reads)
  audit$junction_count <- nrow(result$junctions)
  audit$mean_depth <- mean(result$depth$depth)
  list(status = "SUCCESS", error = NA_character_, source_state = source_state, audit = audit,
       depth = as.numeric(result$depth$depth),
       junctions = result$junctions[, c("key", "chrom", "intron_start1", "intron_end1", "strand", "score"), drop = FALSE])
}

analyze_variant_all <- function(resources, variant, cfg, progress_file = NULL,
                                checkpoint_dir = NULL, resume = TRUE, workers = 1L) {
  progress <- list(state = "RUNNING", total = 0L, completed = 0L, success = 0L, failed = 0L,
                   resumed = 0L, invalidated = 0L, updated_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
  publish_progress <- function() {
    progress$updated_at <<- format(Sys.time(), tz = "UTC", usetz = TRUE)
    variant_write_progress(progress, progress_file)
  }
  complete <- FALSE
  on.exit(if (!complete) { progress$state <- "FAILED"; try(publish_progress(), silent = TRUE) }, add = TRUE)
  workers <- int_scalar(workers, "Concurrent RNA workers", 1, 16)
  slots <- suppressWarnings(as.integer(Sys.getenv("NSLOTS", "")))
  if (!is.na(slots) && (workers > slots || (workers > 1L && workers + 1L > slots)))
    stop("Concurrent RNA workers exceed the allocated NSLOTS after reserving one coordinator slot.")
  if (workers > 1L && .Platform$OS.type != "unix") stop("Parallel RNA workers require Unix; use workers=1 on this platform.")
  variant_check_resources(resources)
  if (!is.data.frame(variant) || nrow(variant) != 1L ||
      !all(c("record_id", "chrom", "pos1", "ref", "alt", "build", "vcf") %in% names(variant)))
    stop("Select exactly one CHROM/POS/REF/ALT record from the current lookup.")
  cfg <- validate_config(cfg)
  if (cfg$demo || cfg$chrom != variant$chrom || variant$pos1 < cfg$start1 || variant$pos1 > cfg$end1)
    stop("The real RNA interval must contain the selected variant on exactly the same chromosome.")
  current <- lookup_variants(resources, paste0(variant$chrom, ":", variant$pos1))
  current <- current[current$record_id == variant$record_id, , drop = FALSE]
  if (nrow(current) != 1L || !identical(as.character(current$vcf), as.character(variant$vcf)) ||
      !identical(as.character(current$build), as.character(variant$build)))
    stop("The selected variant or resource identity changed; repeat lookup.")
  variant <- current
  info <- variant_vcf_info(variant$vcf, variant$chrom)
  calls <- variant_get_genotypes(resources, variant)
  selected <- calls[calls$linked & calls$call_status == "CALLED", , drop = FALSE]
  if (!nrow(selected)) stop("No called DNA genotype has an explicitly linked RNA BAM. Missing and partial DNA calls are not used for RNA genotype comparison.")
  selected <- selected[order(selected$genotype, selected$vcf_sample, method = "radix"), , drop = FALSE]
  rownames(selected) <- NULL
  progress$total <- nrow(selected)
  progress$missing_link_n <- sum(!calls$linked)
  progress$incomplete_gt_n <- sum(calls$call_status != "CALLED")
  progress$group_summary <- data.frame(genotype = sort(unique(calls$genotype)))
  for (name in c("dna_n", "linked_n", "selected_n")) progress$group_summary[[name]] <- 0L
  for (i in seq_len(nrow(progress$group_summary))) {
    group <- progress$group_summary$genotype[[i]]
    progress$group_summary$dna_n[[i]] <- sum(calls$genotype == group)
    progress$group_summary$linked_n[[i]] <- sum(calls$genotype == group & calls$linked)
    progress$group_summary$selected_n[[i]] <- sum(selected$genotype == group)
  }
  progress$group_summary$completed <- progress$group_summary$success <- progress$group_summary$failed <- 0L
  publish_progress()
  version <- variant_run("--version")
  executables <- unname(Sys.which(c("bcftools", "samtools", "regtools")))
  bam_states <- lapply(selected$bam, variant_bam_state)
  binding <- list(schema = "regshiny_full_cohort_v1", code = variant_code_identity(),
    variant = as.list(variant[1L, c("record_id", "chrom", "pos1", "ref", "alt", "build", "vcf")]),
    vcf_state = info$state, vcf_reference_header = info$reference_header,
    manifests = resources$provenance[c("vcf_manifest", "sample_manifest", "vcf_manifest_md5", "sample_manifest_md5", "build")],
    config = cfg, selected_samples = selected[, c("vcf_sample", "bam", "genotype", "raw_gt")],
    bam_states = bam_states,
    bcftools_version = version$text, executable_state = variant_file_state(executables[nzchar(executables)]))
  fingerprint <- variant_hash_object(binding)
  if (!is.null(checkpoint_dir)) {
    if (length(checkpoint_dir) != 1L || !nzchar(checkpoint_dir)) stop("Provide one checkpoint directory.")
    if (!dir.exists(checkpoint_dir) && !dir.create(checkpoint_dir, recursive = TRUE, mode = "0700"))
      stop("Could not create checkpoint directory.")
    checkpoint_dir <- normalizePath(checkpoint_dir, mustWork = TRUE)
    run_path <- file.path(checkpoint_dir, "run.rds")
    if (file.exists(run_path)) {
      old <- tryCatch(readRDS(run_path), error = function(e) NULL)
      if (is.null(old) || !identical(old$fingerprint, fingerprint) ||
          !identical(variant_hash_object(old$binding), fingerprint))
        stop("Checkpoint fingerprint mismatch: inputs, configuration, tools or code changed. Use a new checkpoint directory.")
    } else variant_atomic_rds(list(fingerprint = fingerprint, binding = binding,
                                  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE)), run_path)
  }
  selected$status <- "PENDING"; selected$error <- NA_character_
  selected$retained_reads <- selected$junction_count <- selected$mean_depth <- selected$rna_valid_bases <- NA_real_
  selected$native_audit_passed <- NA
  for (name in c("skipped_N", "deleted_D", "low_baseq", "missing_sequence", "missing_quality"))
    selected[[paste0("rna_excluded_", name)]] <- NA_real_
  selected$checkpoint_reused <- FALSE
  alleles <- c(variant$ref, strsplit(variant$alt, ",", fixed = TRUE)[[1L]])
  rna_supported <- all(grepl("^[ACGTN]$", alleles))
  stream <- variant_stream_init(calls, cfg, rna_supported)
  audits <- vector("list", nrow(selected)); warnings <- info$warnings
  if (!rna_supported) warnings <- c(warnings, "RNA base comparison is available only for SNV records. Indels and symbolic alleles retain DNA GT, junction and depth evidence; RNA allele support is unavailable.")
  if (any(!calls$linked)) warnings <- c(warnings, paste(sum(!calls$linked), "VCF samples have no configured RNA BAM and are not analyzed."))
  if (any(calls$call_status != "CALLED")) warnings <- c(warnings, "Missing and partial DNA calls are shown in the audit and are excluded from RNA genotype comparison.")
  compute_one <- function(i) {
    sample_cfg <- cfg; sample_cfg$bam <- selected$bam[[i]]; sample_cfg$demo <- FALSE
    sample_cfg$variant_position <- if (rna_supported) variant$pos1 else NULL
    source_state <- variant_bam_state(sample_cfg$bam)
    if (!identical(source_state, bam_states[[i]])) stop("A BAM or index changed after the source snapshot was frozen.")
    key <- variant_hash_object(list(sample = selected$vcf_sample[[i]], bam = selected$bam[[i]],
                                    genotype = selected$genotype[[i]], raw_gt = selected$raw_gt[[i]]))
    path <- if (is.null(checkpoint_dir)) NULL else file.path(checkpoint_dir, paste0("sample_", key, ".rds"))
    contribution <- NULL; reused <- FALSE; invalidated <- FALSE
    if (isTRUE(resume) && !is.null(path) && file.exists(path)) {
      saved <- tryCatch(readRDS(path), error = function(e) NULL)
      if (!is.null(saved) && identical(saved$fingerprint, fingerprint) && identical(saved$sample_key, key) &&
          identical(saved$checksum, variant_hash_object(saved$contribution)) &&
          identical(saved$contribution$source_state, source_state) && saved$contribution$status %in% c("SUCCESS", "FAILED")) {
        contribution <- saved$contribution
        reused <- TRUE
      } else invalidated <- TRUE
    }
    if (is.null(contribution)) {
      contribution <- tryCatch({
        contigs <- read_bam_contigs(sample_cfg$bam)
        m <- match(variant$chrom, contigs$chrom)
        if (is.na(m)) stop("The BAM has no exact matching chromosome name.")
        if (!is.na(info$contig_length) && contigs$length_bp[[m]] != info$contig_length)
          stop("VCF and BAM contig lengths differ; reference assembly compatibility failed.")
        full_result <- analyze_bam(sample_cfg)
        variant_compact_contribution(full_result, source_state)
      }, error = function(e) list(status = "FAILED", error = conditionMessage(e),
                                  source_state = source_state, audit = NULL, depth = NULL, junctions = NULL))
      if (!identical(source_state, variant_bam_state(sample_cfg$bam)))
        stop("The BAM or its index changed during analysis; no checkpoint was committed.")
      # Per-sample failure is itself a committed, explicitly missing result.
      if (!is.null(path)) variant_atomic_rds(list(fingerprint = fingerprint, sample_key = key,
        checksum = variant_hash_object(contribution), contribution = contribution), path)
    }
    list(contribution = contribution, reused = reused, invalidated = invalidated)
  }
  batches <- split(seq_len(nrow(selected)), ceiling(seq_len(nrow(selected)) / workers))
  for (batch in batches) {
    completed <- if (workers == 1L) lapply(batch, compute_one) else
      parallel::mclapply(batch, compute_one, mc.cores = min(workers, length(batch)),
                        mc.preschedule = FALSE, mc.set.seed = FALSE, mc.cleanup = TRUE)
    if (any(vapply(completed, function(x) is.null(x) || inherits(x, "try-error"), logical(1))))
      stop("An RNA worker exited without a committed contribution; resume the run after resolving the worker failure.")
    for (k in seq_along(batch)) {
    i <- batch[[k]]; entry <- completed[[k]]; contribution <- entry$contribution
    selected$checkpoint_reused[[i]] <- entry$reused
    progress$resumed <- progress$resumed + as.integer(entry$reused)
    progress$invalidated <- progress$invalidated + as.integer(entry$invalidated)
    variant_stream_add(stream, selected$genotype[[i]], contribution)
    selected$status[[i]] <- contribution$status; selected$error[[i]] <- contribution$error
    if (identical(contribution$status, "SUCCESS")) {
      a <- contribution$audit; audits[[i]] <- a
      selected$retained_reads[[i]] <- a$retained_reads; selected$junction_count[[i]] <- a$junction_count
      selected$mean_depth[[i]] <- a$mean_depth; selected$native_audit_passed[[i]] <- a$native_audit_passed
      if (rna_supported) {
        selected$rna_valid_bases[[i]] <- sum(a$rna_bases$count)
        for (name in names(attr(a$rna_bases, "excluded")))
          selected[[paste0("rna_excluded_", name)]][[i]] <- unname(attr(a$rna_bases, "excluded")[[name]])
      }
      warnings <- unique(c(warnings, a$warnings)); progress$success <- progress$success + 1L
    } else progress$failed <- progress$failed + 1L
    progress$completed <- progress$completed + 1L
    progress$current_sample <- selected$vcf_sample[[i]]
    g <- match(selected$genotype[[i]], progress$group_summary$genotype)
    progress$group_summary$completed[[g]] <- progress$group_summary$completed[[g]] + 1L
    field <- if (contribution$status == "SUCCESS") "success" else "failed"
    progress$group_summary[[field]][[g]] <- progress$group_summary[[field]][[g]] + 1L
    }
    publish_progress()
    # Only at most 'workers' compact contributions are held between batches.
    completed <- contribution <- entry <- NULL
  }
  if (!identical(info$state, variant_file_state(info$state$path)))
    stop("The VCF or its index changed during analysis; results were not published.")
  variant_check_resources(resources)
  if (!all(vapply(seq_len(nrow(selected)), function(i)
    identical(bam_states[[i]], variant_bam_state(selected$bam[[i]])), logical(1))))
    stop("A BAM or index changed during the cohort run; results were not published.")
  aggregated <- variant_stream_finish(stream, calls, selected)
  if (any(selected$status != "SUCCESS")) warnings <- c(warnings, "One or more selected RNA analyses failed; their evidence remains NA and they are excluded from successful-sample means. Full cohort processing completed with these explicitly recorded failures.")
  provenance <- c(resources$provenance, list(record_id = variant$record_id, sampling_mode = "all",
    vcf_state = info$state, vcf_reference_header = info$reference_header,
    selection = list(mode = "all", order = "genotype, then exact VCF sample ID (radix sort)",
                     max_per_group = NULL, max_total = NULL, replacement_on_failure = FALSE),
    fingerprint = fingerprint, code_identity = binding$code,
    execution = list(workers = workers, allocated_slots = if (is.na(slots)) NULL else slots,
                     source_consistency = "All selected BAM/index states frozen before analysis and checked at completion."),
    checkpoint = list(directory = checkpoint_dir, resume_requested = isTRUE(resume),
                      reused = progress$resumed, invalidated = progress$invalidated,
                      policy = "Per-sample atomic contributions replayed exactly once; unchanged failed checkpoints are retained; resume=FALSE recomputes all samples."),
    sample_result_format = "Compact per-sample audit; detailed reads/events and depth are not retained in the final result.",
    rna_evidence_status = if (rna_supported) "SNV_BASE_COUNTS" else "UNSUPPORTED_NON_SNV",
    count_unit = "retained primary alignment record; paired ends separate",
    aggregate_definition = "Unweighted mean across all successfully analyzed linked called-genotype BAMs, including observed zeros; failures excluded and counted.",
    bcftools_version = version$text, config = cfg,
    completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE), log = c(info$log, attr(calls, "log"))))
  result <- c(list(variant = variant, genotypes = calls, samples = selected, sample_results = audits,
                   provenance = provenance, warnings = unique(warnings)), aggregated)
  progress$state <- "COMPLETE"; publish_progress(); complete <- TRUE
  result
}

analyze_variant <- function(resources, variant, cfg, max_per_group = 3L,
                            sampling_mode = c("all", "preview"), progress_file = NULL,
                            checkpoint_dir = NULL, resume = TRUE, workers = 1L) {
  # Older, already-running Shiny sessions pass an explicit cap but no mode.
  # Preserve that bounded request while new calls without a cap default to all.
  if (missing(sampling_mode) && !missing(max_per_group) &&
      is.null(checkpoint_dir) && is.null(progress_file)) sampling_mode <- "preview"
  sampling_mode <- match.arg(sampling_mode)
  if (sampling_mode == "preview") {
    if (!is.null(checkpoint_dir)) stop("Persistent checkpoints are available for full-cohort mode only.")
    return(analyze_variant_preview(resources, variant, cfg, max_per_group = max_per_group))
  }
  analyze_variant_all(resources, variant, cfg, progress_file = progress_file,
                      checkpoint_dir = checkpoint_dir, resume = resume, workers = workers)
}

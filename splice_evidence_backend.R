# Independent reference overlay: no BAM/VCF recount; no allele-phased RNA claim.
splice_int <- function(x, name, lower = 1L, upper = 2e9) {
  if (length(x) != 1L || is.na(x) || !is.numeric(x) || !is.finite(x) ||
      x != floor(x) || x < lower || x > upper) stop(name, " must be an integer from ", lower, " to ", upper, ".")
  as.integer(x)
}

splice_file_state <- function(paths) {
  x <- file.info(paths)
  data.frame(path = paths, size_bytes = unname(x$size), mtime_epoch = as.numeric(x$mtime), stringsAsFactors = FALSE)
}

splice_hash_file <- function(path) {
  tool <- unname(Sys.which("sha256sum")); args <- path
  if (!nzchar(tool)) { tool <- unname(Sys.which("shasum")); args <- c("-a", "256", path) }
  if (!nzchar(tool)) return(list(algorithm = "MD5", value = unname(tools::md5sum(path))))
  p <- processx::run(tool, args, timeout = 60, error_on_status = FALSE)
  if (p$status != 0L) stop("Cannot hash reference evidence: ", p$stderr)
  value <- strsplit(trimws(p$stdout), "[[:space:]]+")[[1L]][[1L]]
  if (!grepl("^[0-9a-fA-F]{64}$", value)) stop("Unexpected SHA256 output.")
  list(algorithm = "SHA256", value = tolower(value))
}

load_splice_reference <- function(fasta = Sys.getenv("REGSHINY_REFERENCE_FASTA"),
                                  receipt = Sys.getenv("REGSHINY_REFERENCE_RECEIPT")) {
  if (length(fasta) != 1L || is.na(fasta) || !nzchar(fasta) || !file.exists(fasta) || file.access(fasta, 4) != 0L)
    stop("GRCh38 reference FASTA is not configured or readable.")
  fasta <- normalizePath(fasta, mustWork = TRUE)
  if (grepl("\\.gz$", fasta, ignore.case = TRUE)) stop("This overlay requires a plain indexed FASTA.")
  fai <- paste0(fasta, ".fai")
  if (!file.exists(fai) || file.access(fai, 4) != 0L) stop("An existing readable FASTA .fai index is required; no source index is created.")
  state <- splice_file_state(c(fasta, fai))
  x <- utils::read.delim(fai, header = FALSE, sep = "\t", quote = "", comment.char = "", colClasses = "character")
  if (ncol(x) < 5L || !nrow(x)) stop("Invalid FASTA index schema.")
  lengths <- suppressWarnings(as.numeric(x[[2L]]))
  if (anyDuplicated(x[[1L]]) || any(!nzchar(x[[1L]])) || anyNA(lengths) ||
      any(!is.finite(lengths)) || any(lengths < 1 | lengths != floor(lengths))) stop("Invalid FASTA contig names or lengths.")
  fai_hash <- splice_hash_file(fai); prior <- NULL
  if (length(receipt) == 1L && !is.na(receipt) && nzchar(receipt)) {
    if (!file.exists(receipt)) stop("The configured FASTA provenance receipt does not exist.")
    prior <- jsonlite::fromJSON(receipt, simplifyVector = TRUE)
    needed <- c("fasta_path", "size_bytes", "mtime_epoch", "sha256", "fai_sha256")
    if (!all(needed %in% names(prior)) || !identical(normalizePath(prior$fasta_path, mustWork = TRUE), fasta) ||
        !isTRUE(as.numeric(prior$size_bytes) == state$size_bytes[[1L]]) ||
        !isTRUE(as.numeric(prior$mtime_epoch) == state$mtime_epoch[[1L]]) ||
        !grepl("^[0-9a-fA-F]{64}$", prior$sha256) || fai_hash$algorithm != "SHA256" ||
        !identical(tolower(prior$fai_sha256), fai_hash$value))
      stop("The prior FASTA hash receipt does not match current FASTA/FAI identity and stat.")
  }
  if (!identical(state, splice_file_state(c(fasta, fai)))) stop("Reference resources changed during initialization.")
  structure(list(fasta = fasta, fai = fai, build = "GRCh38", contigs = setNames(lengths, x[[1L]]),
    state = state, fai_hash = fai_hash, prior_receipt = prior), class = "regshiny_splice_reference")
}

splice_run_faidx <- function(fasta, region) {
  tool <- unname(Sys.which("samtools"))
  if (!nzchar(tool)) stop("samtools is required for reference sequence queries.")
  p <- processx::run(tool, c("faidx", fasta, region), timeout = 60, error_on_status = FALSE)
  if (p$status != 0L || isTRUE(p$timeout)) stop("Reference query failed: ", p$stderr)
  p$stdout
}

splice_fetch_reference <- function(reference, chrom, start1, end1) {
  if (!inherits(reference, "regshiny_splice_reference")) stop("Load a valid GRCh38 reference FASTA first.")
  if (!identical(reference$state, splice_file_state(c(reference$fasta, reference$fai))))
    stop("Reference FASTA or index changed; reload the reference resource.")
  if (!chrom %in% names(reference$contigs)) stop("Reference has no exact chromosome name ", chrom, ". No aliases are inferred.")
  start1 <- splice_int(start1, "Reference query start")
  end1 <- splice_int(end1, "Reference query end", start1, reference$contigs[[chrom]])
  if (end1 - start1 + 1L > 1000L) stop("Reference display query exceeds 1,000 bp.")
  # A private alias/index isolates even automatic faidx index maintenance.
  work <- tempfile("regshiny_reference_"); dir.create(work, mode = "0700")
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  alias <- file.path(work, "reference.fa")
  if (!file.symlink(reference$fasta, alias) || !file.copy(reference$fai, paste0(alias, ".fai")))
    stop("Could not stage a private reference query.")
  Sys.chmod(paste0(alias, ".fai"), "0600")
  region <- paste0(chrom, ":", start1, "-", end1)
  lines <- strsplit(splice_run_faidx(alias, region), "\n", fixed = TRUE)[[1L]]
  if (!length(lines) || !startsWith(lines[[1L]], ">")) stop("Unexpected faidx response.")
  sequence <- toupper(paste(lines[!startsWith(lines, ">")], collapse = ""))
  if (nchar(sequence) != end1 - start1 + 1L || !grepl("^[ACGTRYSWKMBDHVN]+$", sequence))
    stop("Reference response length/alphabet does not match the requested interval.")
  digest_path <- file.path(work, "query.txt"); writeLines(c(region, sequence), digest_path, useBytes = TRUE)
  hash <- splice_hash_file(digest_path)
  if (!identical(reference$state, splice_file_state(c(reference$fasta, reference$fai))))
    stop("Reference FASTA or index changed during query; sequence evidence was discarded.")
  list(chrom = chrom, start1 = start1, end1 = end1, sequence = sequence, query_hash = hash)
}

splice_revcomp <- function(sequence) {
  paste(rev(strsplit(chartr("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN", toupper(sequence)), "", fixed = TRUE)[[1L]]), collapse = "")
}

splice_empty_sites <- function() data.frame(candidate_id = character(), chrom = character(), site_type = character(),
  strand = character(), start1 = integer(), end1 = integer(), anchor1 = integer(), canonical_expected = character(),
  source = character(), transcript_ids = character(), context_transcript_ids = character(), junction_ids = character(),
  distance_to_variant = integer(), variant_overlaps = logical(), allele_index = integer(), allele = character(),
  ref_motif = character(), alt_motif = character(), ref_canonical = logical(), alt_canonical = logical(),
  change = character(), sequence_only = logical(), display_default = logical(), stringsAsFactors = FALSE)

build_splice_evidence <- function(result, annotation, reference, flank = 25L) {
  flank <- splice_int(flank, "Splice display flank", 10L, 100L)
  v <- result$variant
  if (!is.data.frame(v) || nrow(v) != 1L || !all(c("chrom", "pos1", "ref", "alt", "build") %in% names(v)))
    stop("Splice evidence needs one exact VCF variant record.")
  chrom <- as.character(v$chrom[[1L]]); position <- splice_int(v$pos1[[1L]], "Variant position")
  if (!grepl("^[A-Za-z0-9_][A-Za-z0-9_.-]*$", chrom)) stop("Invalid variant chromosome.")
  ref_allele <- as.character(v$ref[[1L]]); alt_alleles <- strsplit(as.character(v$alt[[1L]]), ",", fixed = TRUE)[[1L]]
  alleles <- c(ref_allele, alt_alleles)
  if (!length(alt_alleles)) stop("The selected variant has no ALT allele.")
  warnings <- character(); reference_status <- "UNAVAILABLE_REFERENCE"; sequence_evidence <- NULL
  chr_length <- if ("contig_length" %in% names(v) && length(v$contig_length) == 1L && is.finite(v$contig_length)) v$contig_length[[1L]] else 2e9
  ref_ok <- inherits(reference, "regshiny_splice_reference")
  if (chr_length == 2e9 && ref_ok && chrom %in% names(reference$contigs) && reference$contigs[[chrom]] >= position)
    chr_length <- reference$contigs[[chrom]]
  if (position > chr_length) stop("Variant position exceeds the declared contig length.")
  cfg <- result$provenance$config
  if (is.null(cfg) || !all(c("chrom", "start1", "end1") %in% names(cfg))) {
    available <- result$depth$pos1[result$depth$chrom == chrom]
    available <- available[is.finite(available)]
    if (!length(available)) stop("The loaded RNA interval cannot be established from configuration or coverage.")
    cfg <- list(chrom = chrom, start1 = min(available), end1 = max(available))
    warnings <- c(warnings, "Loaded RNA interval inferred from existing depth coordinates because saved configuration is absent.")
  }
  roi_start <- splice_int(cfg$start1, "Loaded RNA interval start")
  roi_end <- splice_int(cfg$end1, "Loaded RNA interval end", roi_start)
  if (!identical(as.character(cfg$chrom), chrom) || position < roi_start || position > roi_end)
    stop("The selected variant is outside the loaded RNA result interval.")
  requested_start <- max(1L, position - flank); requested_end <- min(chr_length, position + flank)
  window <- list(chrom = chrom, start1 = max(roi_start, requested_start), end1 = min(roi_end, requested_end),
    variant_pos1 = position, build = as.character(v$build[[1L]]), requested_flank = flank,
    clipped = requested_start < roi_start || requested_end > roi_end,
    clip_reason = if (requested_start < roi_start || requested_end > roi_end) "Display constrained to the loaded RNA analysis interval" else NULL)
  positions <- seq.int(window$start1, window$end1)
  groups <- result$group_summary
  if (!is.data.frame(groups) || !all(c("genotype", "analyzed_n") %in% names(groups)) || anyDuplicated(groups$genotype))
    stop("Genotype sample accounting is missing or ambiguous.")
  if (!"call_status" %in% names(groups)) groups$call_status <- ifelse(grepl("^[0-9]+([/|][0-9]+)*$", groups$genotype), "CALLED", "UNAVAILABLE")
  genotype_indices <- lapply(seq_len(nrow(groups)), function(i) {
    if (groups$call_status[[i]] != "CALLED" || !grepl("^[0-9]+([/|][0-9]+)*$", groups$genotype[[i]])) return(integer())
    index <- suppressWarnings(as.integer(strsplit(groups$genotype[[i]], "[/|]")[[1L]]))
    if (anyNA(index) || any(index >= length(alleles))) return(integer())
    index
  })
  genotypes <- groups
  genotypes$dna_alleles <- vapply(genotype_indices, function(index) if (length(index)) paste(alleles[index + 1L], collapse = "/") else NA_character_, character(1))
  genotypes$ploidy <- vapply(genotype_indices, length, integer(1)); genotypes$motif_summary <- NA_character_
  depth <- result$depth
  if (!is.data.frame(depth) || !all(c("genotype", "chrom", "pos1", "mean_depth") %in% names(depth)))
    stop("Existing per-base genotype coverage is required; no BAM recount is performed.")
  coverage <- do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
    d <- depth[depth$genotype == groups$genotype[[i]] & depth$chrom == chrom, , drop = FALSE]
    if (anyDuplicated(d$pos1)) stop("Existing coverage contains duplicate genotype-position rows.")
    m <- match(positions, d$pos1); measured <- !is.na(m); values <- d$mean_depth[m]
    if (any(values < 0, na.rm = TRUE)) stop("Existing coverage contains negative values.")
    data.frame(genotype = groups$genotype[[i]], pos1 = positions, mean_depth = values,
      analyzed_n = groups$analyzed_n[[i]], coverage_status = ifelse(!measured, "UNMEASURED",
        ifelse(is.na(values), "UNAVAILABLE", ifelse(values == 0, "OBSERVED_ZERO", "OBSERVED"))), stringsAsFactors = FALSE)
  }))
  if (!ref_ok) {
    warnings <- c(warnings, if (inherits(reference, "error")) conditionMessage(reference) else "Reference FASTA is unavailable; coverage and RNA junction counts remain available.")
  } else if (!identical(as.character(v$build[[1L]]), reference$build)) {
    reference_status <- "REFERENCE_BUILD_MISMATCH"; warnings <- c(warnings, "Variant and reference assemblies differ. Motif inference is blocked; no liftover is performed.")
  } else if (!chrom %in% names(reference$contigs)) {
    reference_status <- "REFERENCE_CONTIG_MISSING"; warnings <- c(warnings, "Reference chromosome name does not match the variant exactly.")
  } else if ("contig_length" %in% names(v) && is.finite(v$contig_length[[1L]]) && v$contig_length[[1L]] != reference$contigs[[chrom]]) {
    reference_status <- "REFERENCE_LENGTH_MISMATCH"; warnings <- c(warnings, "Variant and FASTA contig lengths differ. Motif inference is blocked.")
  } else {
    sequence_evidence <- tryCatch(splice_fetch_reference(reference, chrom, max(1L, window$start1 - 1L),
      min(reference$contigs[[chrom]], window$end1 + 1L)), error = function(e) e)
    if (inherits(sequence_evidence, "error")) {
      warnings <- c(warnings, conditionMessage(sequence_evidence)); sequence_evidence <- NULL
    } else if (!all(grepl("^[ACGT]$", alleles))) {
      reference_status <- "UNAVAILABLE_NON_SNV_OR_AMBIGUOUS_ALLELE"
      warnings <- c(warnings, "Motif reconstruction supports unambiguous SNVs only; indels, symbolic and ambiguous alleles are unavailable.")
    } else {
      actual <- substr(sequence_evidence$sequence, position - sequence_evidence$start1 + 1L, position - sequence_evidence$start1 + 1L)
      if (!identical(actual, ref_allele)) {
        reference_status <- "REF_MISMATCH"
        warnings <- c(warnings, paste0("VCF REF ", ref_allele, " differs from FASTA base ", actual, "; all ALT/motif inference is blocked."))
      } else reference_status <- "REF_MATCH"
    }
  }
  subseq <- function(start, end) {
    if (is.null(sequence_evidence) || start < sequence_evidence$start1 || end > sequence_evidence$end1) return(NA_character_)
    substr(sequence_evidence$sequence, start - sequence_evidence$start1 + 1L, end - sequence_evidence$start1 + 1L)
  }
  bases <- data.frame(pos1 = positions, ref_base = vapply(positions, function(p) subseq(p, p), character(1)),
                      alt_base = NA_character_, alt_bases = NA_character_, is_variant = positions == position)
  sequences <- data.frame(allele_index = seq_along(alleles) - 1L, allele = alleles,
                          sequence = NA_character_, available = FALSE, stringsAsFactors = FALSE)
  display_ref <- subseq(window$start1, window$end1)
  if (reference_status == "REF_MATCH") {
    sequences$sequence[[1L]] <- display_ref; sequences$available[[1L]] <- TRUE
    bases$alt_base <- bases$ref_base
    bases$alt_base[bases$is_variant] <- if (length(alt_alleles) == 1L) alt_alleles else NA_character_
    bases$alt_bases[bases$is_variant] <- paste(paste0(seq_along(alt_alleles), ":", alt_alleles), collapse = ";")
    for (i in seq_along(alt_alleles)) {
      sequence <- display_ref
      substr(sequence, position - window$start1 + 1L, position - window$start1 + 1L) <- alt_alleles[[i]]
      sequences$sequence[[i + 1L]] <- sequence; sequences$available[[i + 1L]] <- TRUE
    }
  }
  j <- result$junctions
  required_j <- c("chrom", "intron_start1", "intron_end1", "strand")
  if (!is.data.frame(j) || !all(required_j %in% names(j))) stop("Existing RNA junction boundaries/strand are required.")
  junction_support <- j
  junction_support$junction_id <- paste(j$chrom, j$intron_start1, j$intron_end1, j$strand, sep = ":")
  junction_support$distance_to_variant <- pmin(abs(j$intron_start1 - position), abs(j$intron_end1 - position))
  junction_support$within_window <- j$chrom == chrom & j$intron_start1 <= window$end1 & j$intron_end1 >= window$start1
  junction_support$strand_resolved <- !is.na(j$strand) & j$strand %in% c("+", "-")
  candidates <- list()
  add_site <- function(type, strand, start, end, anchor, source, transcript = "", context = "", junction = "") {
    if (start < 1 || end > chr_length || end < window$start1 || start > window$end1) return(invisible(NULL))
    id <- paste(chrom, start, end, strand, type, sep = ":")
    candidates[[length(candidates) + 1L]] <<- data.frame(candidate_id = id, chrom = chrom, site_type = type,
      strand = strand, start1 = start, end1 = end, anchor1 = anchor,
      canonical_expected = if (type == "donor") "GT" else if (type == "acceptor") "AG" else NA_character_,
      source = source, transcript_ids = transcript, context_transcript_ids = context, junction_ids = junction,
      distance_to_variant = min(abs(c(start, end) - position)), variant_overlaps = start <= position && position <= end,
      stringsAsFactors = FALSE)
  }
  add_intron <- function(start, end, strand, source, transcript = "", junction = "") {
    if (!is.finite(start) || !is.finite(end) || end - start + 1L < 2L) return(invisible(NULL))
    if (strand == "+") {
      add_site("donor", strand, start, start + 1L, start, source, transcript, junction = junction)
      add_site("acceptor", strand, end - 1L, end, end, source, transcript, junction = junction)
    } else if (strand == "-") {
      add_site("donor", strand, end - 1L, end, end, source, transcript, junction = junction)
      add_site("acceptor", strand, start, start + 1L, start, source, transcript, junction = junction)
    } else {
      add_site("boundary_start", "?", start, start + 1L, start, source, junction = junction)
      add_site("boundary_end", "?", end - 1L, end, end, source, junction = junction)
    }
  }
  observed <- unique(junction_support[, c(required_j, "junction_id"), drop = FALSE])
  for (i in seq_len(nrow(observed))) if (observed$chrom[[i]] == chrom)
    add_intron(observed$intron_start1[[i]], observed$intron_end1[[i]], if (is.na(observed$strand[[i]])) "?" else observed$strand[[i]],
               "RNA_OBSERVED", junction = observed$junction_id[[i]])
  annotation_ok <- is.list(annotation) && !inherits(annotation, "error") &&
    is.data.frame(annotation$introns) && is.data.frame(annotation$transcripts) &&
    identical(annotation$metadata$build, "GRCh38") && identical(annotation$metadata$build, as.character(v$build[[1L]])) &&
    isTRUE(annotation$metadata$query$full_transcript_models) &&
    identical(annotation$metadata$query$truncated, FALSE) && identical(annotation$metadata$query$chrom, chrom) &&
    isTRUE(annotation$metadata$query$start1 <= position) && isTRUE(annotation$metadata$query$end1 >= position)
  if (!annotation_ok) warnings <- c(warnings, "Complete matching GENCODE models are unavailable; only observed RNA boundaries are evaluated and no sequence-only candidate sites are scanned.")
  if (annotation_ok) {
    introns <- annotation$introns
    for (i in seq_len(nrow(introns))) if (introns$chrom[[i]] == chrom && introns$strand[[i]] %in% c("+", "-"))
      add_intron(introns$intron_start1[[i]], introns$intron_end1[[i]], introns$strand[[i]], "GENCODE", introns$transcript_id[[i]])
    tx <- annotation$transcripts
    tx <- tx[tx$chrom == chrom & tx$start1 <= position & tx$end1 >= position & tx$strand %in% c("+", "-"), , drop = FALSE]
    for (strand in unique(tx$strand)) for (start in c(position - 1L, position)) for (type in c("donor", "acceptor")) {
      end <- start + 1L
      anchor <- if ((type == "donor" && strand == "+") || (type == "acceptor" && strand == "-")) start else end
      add_site(type, strand, start, end, anchor, "SEQUENCE_ONLY", context = paste(sort(unique(tx$transcript_id[tx$strand == strand])), collapse = ";"))
    }
  }
  sites <- splice_empty_sites()
  if (length(candidates)) {
    raw <- do.call(rbind, candidates)
    raw <- do.call(rbind, lapply(split(raw, raw$candidate_id), function(rows) {
      one <- rows[1L, , drop = FALSE]
      for (field in c("source", "transcript_ids", "context_transcript_ids", "junction_ids"))
        one[[field]] <- paste(sort(unique(rows[[field]][nzchar(rows[[field]])])), collapse = ";")
      sources <- strsplit(one$source, ";", fixed = TRUE)[[1L]]
      if (length(sources) > 1L) one$source <- paste(sources[sources != "SEQUENCE_ONLY"], collapse = ";")
      one
    }))
    evaluated <- list()
    for (k in seq_len(nrow(raw))) for (a in seq_along(alt_alleles)) {
      site <- raw[k, , drop = FALSE]; site$allele_index <- a; site$allele <- alt_alleles[[a]]
      site$sequence_only <- identical(site$source, "SEQUENCE_ONLY")
      site$ref_motif <- site$alt_motif <- NA_character_; site$ref_canonical <- site$alt_canonical <- NA
      site$change <- "UNAVAILABLE_REFERENCE"
      if (reference_status != "REF_MATCH") site$change <- paste0("UNAVAILABLE_", reference_status) else if (!site$strand %in% c("+", "-")) {
        site$change <- "UNAVAILABLE_UNKNOWN_STRAND"
      } else {
        ref <- subseq(site$start1, site$end1); alt <- ref
        if (!is.na(ref) && site$variant_overlaps) substr(alt, position - site$start1 + 1L, position - site$start1 + 1L) <- alt_alleles[[a]]
        if (site$strand == "-" && !is.na(ref)) { ref <- splice_revcomp(ref); alt <- splice_revcomp(alt) }
        site$ref_motif <- ref; site$alt_motif <- alt
        if (is.na(ref) || !grepl("^[ACGT]{2}$", ref)) site$change <- "UNAVAILABLE_AMBIGUOUS_REFERENCE" else {
          site$ref_canonical <- identical(ref, site$canonical_expected[[1L]])
          site$alt_canonical <- identical(alt, site$canonical_expected[[1L]])
          site$change <- if (!site$variant_overlaps) "NOT_OVERLAPPING" else
            if (site$ref_canonical && !site$alt_canonical) "DISRUPTED_CANONICAL" else
            if (!site$ref_canonical && site$alt_canonical) "CREATED_CANONICAL" else
            if (site$ref_canonical && site$alt_canonical) "PRESERVED_CANONICAL" else
            if (identical(ref, alt)) "UNCHANGED_NONCANONICAL" else "CHANGED_NONCANONICAL"
        }
      }
      site$display_default <- site$variant_overlaps && (!site$sequence_only || site$change %in% c("CREATED_CANONICAL", "DISRUPTED_CANONICAL"))
      evaluated[[length(evaluated) + 1L]] <- site
    }
    if (length(evaluated)) {
      sites <- do.call(rbind, evaluated)
      retained_ids <- unique(sites$candidate_id[!sites$sequence_only |
        sites$change %in% c("CREATED_CANONICAL", "DISRUPTED_CANONICAL")])
      sites <- sites[sites$candidate_id %in% retained_ids, , drop = FALSE]
    }
    rownames(sites) <- NULL
  }
  genotype_sites <- data.frame(genotype = character(), candidate_id = character(), site_type = character(), strand = character(),
    motif_alleles = character(), canonical_alleles_n = integer(), ploidy = integer(), interpretation = character(), stringsAsFactors = FALSE)
  if (nrow(sites)) {
    rows <- list()
    for (id in unique(sites$candidate_id)) {
      site <- sites[sites$candidate_id == id, , drop = FALSE]
      motif <- setNames(c(site$ref_motif[[1L]], rep(NA_character_, length(alt_alleles))), as.character(0:length(alt_alleles)))
      motif[as.character(site$allele_index)] <- site$alt_motif
      expected <- site$canonical_expected[[1L]]
      for (i in seq_len(nrow(genotypes))) {
        index <- genotype_indices[[i]]; labels <- motif[as.character(index)]
        valid <- length(index) > 0L && !anyNA(labels) && all(grepl("^[ACGT]{2}$", labels)) &&
          !is.na(expected) && reference_status == "REF_MATCH"
        rows[[length(rows) + 1L]] <- data.frame(genotype = genotypes$genotype[[i]], candidate_id = id,
          site_type = site$site_type[[1L]], strand = site$strand[[1L]],
          motif_alleles = if (valid) paste(labels, collapse = "/") else NA_character_,
          canonical_alleles_n = if (valid) sum(labels == expected) else NA_integer_, ploidy = length(index),
          interpretation = if (!valid) "Motif unavailable for this genotype/site" else
            "Motifs reconstructed from DNA alleles; pooled RNA coverage is not allele-phased", stringsAsFactors = FALSE)
      }
    }
    genotype_sites <- do.call(rbind, rows)
  }
  primary_id <- NA_character_
  if (nrow(sites)) {
    resolved <- sites[sites$strand %in% c("+", "-") & sites$variant_overlaps &
      !startsWith(sites$change, "UNAVAILABLE") & reference_status == "REF_MATCH", , drop = FALSE]
    if (nrow(resolved)) {
      source_rank <- ifelse(grepl("RNA_OBSERVED", resolved$source, fixed = TRUE), 0L,
        ifelse(grepl("GENCODE", resolved$source, fixed = TRUE), 1L, 2L))
      priority <- order(!resolved$variant_overlaps, source_rank, resolved$distance_to_variant,
                        resolved$candidate_id, resolved$allele_index)
      primary_id <- resolved$candidate_id[priority[[1L]]]
    }
  }
  genotypes$primary_candidate_id <- primary_id
  for (i in seq_len(nrow(genotypes))) {
    g <- genotype_sites[genotype_sites$genotype == genotypes$genotype[[i]] &
      !is.na(genotype_sites$motif_alleles) & !is.na(primary_id) & genotype_sites$candidate_id == primary_id, , drop = FALSE]
    genotypes$motif_summary[[i]] <- if (nrow(g)) g$motif_alleles[[1L]] else NA_character_
  }
  metadata <- list(reference_status = reference_status,
    primary_candidate_id = primary_id,
    focus_candidate_id = primary_id,
    focus_allele_index = if (!is.na(primary_id)) {
      focused <- sites[sites$candidate_id == primary_id, , drop = FALSE]
      changed <- focused$allele_index[focused$change %in% c("CREATED_CANONICAL", "DISRUPTED_CANONICAL")]
      if (length(changed)) min(changed) else min(focused$allele_index)
    } else NA_integer_,
    variant = as.list(v[1L, c("chrom", "pos1", "ref", "alt", "build")]),
    reference = if (ref_ok) list(fasta = reference$fasta, state = reference$state, fai_hash = reference$fai_hash,
      prior_receipt = reference$prior_receipt, validation = "Exact contig length and selected VCF REF are checked; no liftover or genome-wide sequence identity inference.") else NULL,
    sequence_query = sequence_evidence, annotation_available = annotation_ok,
    annotation = if (annotation_ok) annotation$metadata else NULL,
    window = window, loaded_rna_interval = cfg[c("chrom", "start1", "end1")],
    selection = result$provenance$selection,
    definitions = list(depth = "Existing unweighted mean depth per successfully analyzed genotype-group sample; NA and observed zeros preserved.",
      motifs = "GT donor and AG acceptor in transcript orientation; negative-strand sequence is reverse complemented. Other splice motifs are not classified as canonical here.",
      candidates = "Sequence-only gains/losses require overlapping GENCODE strand context. They are sequence hypotheses, not observed junctions or proof of altered splicing.",
      genotype = "Allele letters and motif pairs come from DNA GT. Heterozygous RNA reads are not phased to either allele.",
      nonoverlap = "A selected variant outside a site's two-base motif is not classified as affecting that motif."),
    completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
  list(window = window, coverage = coverage, bases = bases, sequences = sequences, genotypes = genotypes,
    sites = sites, genotype_sites = genotype_sites, junction_support = junction_support, metadata = metadata,
    warnings = unique(warnings))
}

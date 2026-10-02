# Reg_Shiny computational backend.
# Coordinates: GUI = 1-based closed intron interval; internal = 0-based half-open.
# Counting unit = retained primary alignment record (paired ends count separately).
# Requires processx and data.table; external programs: samtools and regtools.

int_scalar <- function(x, name, lower = 0, upper = 2e9) {
  if (length(x) != 1L || is.na(x) || !is.numeric(x) ||
      !is.finite(x) || x != floor(x) || x < lower || x > upper)
    stop(name, " must be an integer from ", lower, " to ", upper, ".", call. = FALSE)
  as.integer(x)
}

validate_config <- function(cfg) {
  cfg$start1 <- int_scalar(cfg$start1, "Analysis start", 1)
  cfg$end1 <- int_scalar(cfg$end1, "Analysis end", cfg$start1)
  if (cfg$end1 - cfg$start1 + 1 > 250000)
    stop("The maximum analysis interval is 250,000 bp. Select a smaller interval.", call. = FALSE)
  if (length(cfg$chrom) != 1L || is.na(cfg$chrom) ||
      !grepl("^[A-Za-z0-9_.-]+$", cfg$chrom))
    stop("Chromosome names may contain only letters, numbers, dots, underscores and hyphens.", call. = FALSE)
  cfg$mapq <- int_scalar(cfg$mapq, "MAPQ", 0, 255)
  cfg$baseq <- int_scalar(cfg$baseq, "BaseQ", 0, 93)
  cfg$anchor <- int_scalar(cfg$anchor, "Anchor", 1, 1000)
  cfg$min_intron <- int_scalar(cfg$min_intron, "Minimum intron length", 1)
  cfg$max_intron <- int_scalar(cfg$max_intron, "Maximum intron length", cfg$min_intron)
  if (!cfg$strand_mode %in% c("XS", "RF", "FR")) stop("Invalid strandedness mode.")
  cfg$exclude_duplicates <- isTRUE(cfg$exclude_duplicates)
  cfg$nh1_only <- isTRUE(cfg$nh1_only)
  cfg$demo <- isTRUE(cfg$demo)
  cfg
}


parse_bam_contigs <- function(header_lines) {
  sq <- strsplit(header_lines[startsWith(header_lines, "@SQ\t")], "\t", fixed = TRUE)
  if (!length(sq)) stop("The BAM header has no reference sequences (@SQ records).")
  value <- function(fields, prefix) {
    hit <- fields[startsWith(fields, prefix)]
    if (length(hit) != 1L) stop("Malformed reference sequence in the BAM header.")
    substring(hit, nchar(prefix) + 1L)
  }
  chrom <- vapply(sq, value, character(1), prefix = "SN:")
  length_bp <- suppressWarnings(as.numeric(vapply(sq, value, character(1), prefix = "LN:")))
  if (any(!nzchar(chrom)) || anyDuplicated(chrom) || anyNA(length_bp) ||
      any(!is.finite(length_bp)) || any(length_bp < 1) || any(length_bp != floor(length_bp)))
    stop("Invalid reference names or lengths in the BAM header.")
  data.frame(chrom = chrom, length_bp = length_bp, stringsAsFactors = FALSE)
}

read_bam_contigs <- function(bam) {
  if (length(bam) != 1L || !file.exists(bam)) stop("Select an existing BAM file.")
  samtools <- unname(Sys.which("samtools"))
  if (!nzchar(samtools)) stop("samtools is not available in PATH.")
  p <- processx::run(samtools, c("view", "-H", bam), timeout = 20, error_on_status = FALSE)
  if (p$status != 0L) stop("Cannot read the BAM header: ", p$stderr)
  parse_bam_contigs(strsplit(p$stdout, "\n", fixed = TRUE)[[1L]])
}

choose_bam_contig <- function(contigs, current = NULL) {
  if (length(current) == 1L && !is.na(current) && current %in% contigs$chrom) return(current)
  preferred <- c("chr1", "1")
  found <- preferred[preferred %in% contigs$chrom]
  if (length(found)) found[[1L]] else contigs$chrom[[1L]]
}

empty_events <- function() data.frame(
  read_id = integer(), chrom = character(), intron_start0 = integer(),
  intron_end0 = integer(), strand = character(), stringsAsFactors = FALSE)

junction_key <- function(chrom, start0, end0, strand) {
  paste(chrom, start0, end0, strand, sep = "|")
}

# Parse CIGAR without treating N or D as covered reference bases.
parse_cigar <- function(cigar, pos1) {
  empty <- matrix(integer(), ncol = 2L, dimnames = list(NULL, c("start0", "end0")))
  if (identical(cigar, "*")) return(list(blocks = empty, introns = empty))
  tokens <- regmatches(cigar, gregexpr("[0-9]+[MIDNSHP=X]", cigar, perl = TRUE))[[1L]]
  if (!length(tokens) || paste0(tokens, collapse = "") != cigar)
    stop("Unsupported or invalid CIGAR: ", cigar, call. = FALSE)
  lens <- as.integer(sub("[MIDNSHP=X]$", "", tokens))
  if (anyNA(lens) || any(lens <= 0)) stop("Invalid CIGAR operation length.")
  ops <- substring(tokens, nchar(tokens))
  pos0 <- as.integer(pos1) - 1L
  blocks <- list(); introns <- list()
  for (k in seq_along(ops)) {
    op <- ops[[k]]; n <- lens[[k]]
    if (op %in% c("M", "=", "X")) blocks[[length(blocks) + 1L]] <- c(pos0, pos0 + n)
    if (op == "N") introns[[length(introns) + 1L]] <- c(pos0, pos0 + n)
    if (op %in% c("M", "=", "X", "D", "N")) pos0 <- pos0 + n
    # I, S consume query only; H, P consume neither query nor reference.
  }
  bind <- function(x) if (length(x)) {
    z <- do.call(rbind, x); colnames(z) <- c("start0", "end0"); z
  } else empty
  list(blocks = bind(blocks), introns = bind(introns))
}

# Match RegTools flag-based RF/FR convention, including '?' for inconsistent flags.
# No FASTA/intron-motif override is offered in this prototype.
infer_regtools_strand <- function(flag, xs, mode) {
  if (mode == "XS") return(ifelse(xs %in% c("+", "-"), xs, "?"))
  bit <- function(mask) bitwAnd(as.integer(flag), as.integer(mask)) != 0L
  flip <- identical(mode, "RF")
  first <- xor(xor(flip, bit(64)), bit(16))
  second <- xor(xor(flip, bit(128)), bit(32))
  ifelse(first != second, "?", ifelse(first, "+", "-"))
}

parse_sam_records <- function(lines, cfg) {
  z <- strsplit(lines, "\t", fixed = TRUE)
  if (length(z) && any(lengths(z) < 11L)) stop("A SAM record contains fewer than 11 fields.")
  get_col <- function(k) vapply(z, function(x) x[[k]], character(1))
  get_tag <- function(prefix) vapply(z, function(x) {
    optional <- if (length(x) > 11L) x[12:length(x)] else character()
    hit <- optional[startsWith(optional, prefix)]
    if (length(hit)) substring(hit[[1L]], nchar(prefix) + 1L) else NA_character_
  }, character(1))
  reads <- data.frame(
    read_id = seq_along(lines), qname = get_col(1),
    flag = as.integer(get_col(2)), chrom = get_col(3),
    pos1 = as.integer(get_col(4)), mapq = as.integer(get_col(5)),
    cigar = get_col(6), nh = suppressWarnings(as.integer(get_tag("NH:i:"))),
    xs = get_tag("XS:A:"), stringsAsFactors = FALSE)
  reads$strand <- infer_regtools_strand(reads$flag, reads$xs, cfg$strand_mode)
  parsed <- lapply(seq_len(nrow(reads)), function(i)
    parse_cigar(reads$cigar[[i]], reads$pos1[[i]]))
  overlap <- vapply(parsed, function(x) {
    b <- x$blocks
    nrow(b) > 0L && any(b[, 1] < cfg$end1 & b[, 2] > cfg$start1 - 1L)
  }, logical(1))
  keep <- overlap
  if (isTRUE(cfg$nh1_only)) keep <- keep & !is.na(reads$nh) & reads$nh == 1L
  kept <- which(keep)
  ev <- lapply(kept, function(i) {
    x <- parsed[[i]]$introns
    if (!nrow(x)) return(NULL)
    data.frame(read_id = reads$read_id[[i]], chrom = reads$chrom[[i]],
               intron_start0 = x[, 1], intron_end0 = x[, 2],
               strand = reads$strand[[i]], stringsAsFactors = FALSE)
  })
  ev <- Filter(Negate(is.null), ev)
  list(reads = reads[keep, , drop = FALSE],
       events = if (length(ev)) do.call(rbind, ev) else empty_events(),
       keep = keep, span_only_excluded = sum(!overlap),
       missing_nh = sum(is.na(reads$nh)))
}

read_regtools_bed <- function(path) {
  empty <- data.frame(chrom = character(), chromStart = integer(), chromEnd = integer(),
    name = character(), score = integer(), strand = character(), thickStart = integer(),
    thickEnd = integer(), itemRgb = character(), blockCount = integer(),
    blockSizes = character(), blockStarts = character(),
    intron_start0 = integer(), intron_end0 = integer(),
    intron_start1 = integer(), intron_end1 = integer(), key = character())
  if (!file.exists(path) || file.info(path)$size == 0) return(empty)
  b <- data.table::fread(path, header = FALSE, sep = "\t", data.table = FALSE,
                        colClasses = "character", quote = "", na.strings = NULL)
  if (ncol(b) != 12L) stop("RegTools output is not BED12.")
  names(b) <- names(empty)[seq_len(12L)]
  for (n in c("chromStart", "chromEnd", "score", "thickStart", "thickEnd", "blockCount"))
    b[[n]] <- as.integer(b[[n]])
  if (anyNA(b$score) || any(b$score < 0) || any(b$blockCount != 2L))
    stop("Invalid score or blockCount in RegTools BED12 output.")
  sizes <- strsplit(b$blockSizes, ",", fixed = TRUE)
  left <- vapply(sizes, function(x) as.integer(x[[1L]]), integer(1))
  right <- vapply(sizes, function(x) as.integer(x[[2L]]), integer(1))
  b$intron_start0 <- b$chromStart + left
  b$intron_end0 <- b$chromEnd - right
  b$intron_start1 <- b$intron_start0 + 1L
  b$intron_end1 <- b$intron_end0
  b$key <- junction_key(b$chrom, b$intron_start0, b$intron_end0, b$strand)
  b
}

make_demo_sam <- function(path) {
  rec <- function(name, flag, pos, mq, cigar, len = 50L) {
    paste(c(name, flag, "chrDemo", pos, mq, cigar, "*", 0, 0,
            paste(rep("A", len), collapse = ""), paste(rep("I", len), collapse = ""),
            "NH:i:1", "XS:A:+"), collapse = "\t")
  }
  writeLines(c("@HD\tVN:1.6\tSO:unsorted", "@SQ\tSN:chrDemo\tLN:1000",
    rec("exact01", 0, 181, 60, "20M100N30M"),
    rec("exact02", 0, 171, 60, "30M100N20M"),
    rec("near01", 0, 183, 60, "20M100N30M"),
    rec("plain01", 0, 351, 60, "50M"),
    rec("lowmapq01", 0, 181, 5, "20M100N30M"),
    rec("secondary01", 256, 181, 60, "20M100N30M"),
    rec("supplementary01", 2048, 181, 60, "20M100N30M"),
    rec("span_only01", 0, 1, 60, "50M500N50M", 100L)), path)
  invisible(path)
}

default_demo_config <- function() list(
  demo = TRUE, bam = "", chrom = "chrDemo", start1 = 101L, end1 = 450L,
  mapq = 20L, baseq = 0L, anchor = 8L, min_intron = 70L, max_intron = 500000L,
  strand_mode = "XS", exclude_duplicates = FALSE, nh1_only = FALSE)

# Count observed RNA bases at one reference position in the exact retained SAM
# snapshot. SAM SEQ is already in reference orientation, including reverse reads.
# These are alignment observations, not DNA genotypes or a variant caller.
count_rna_bases <- function(records, position1, min_baseq = 0L) {
  position1 <- int_scalar(position1, "Variant position", 1)
  min_baseq <- int_scalar(min_baseq, "RNA base quality", 0, 93)
  counts <- stats::setNames(integer(5), c("A", "C", "G", "T", "N"))
  excluded <- c(skipped_N = 0L, deleted_D = 0L, low_baseq = 0L,
                missing_sequence = 0L, missing_quality = 0L)
  for (record in records) {
    fields <- strsplit(record, "\t", fixed = TRUE)[[1L]]
    if (length(fields) < 11L) stop("Malformed SAM record in RNA base audit.")
    ref <- as.integer(fields[[4L]]); query <- 1L
    tokens <- regmatches(fields[[6L]], gregexpr("[0-9]+[MIDNSHP=X]", fields[[6L]]))[[1L]]
    for (token in tokens) {
      op <- substring(token, nchar(token)); n <- as.integer(substring(token, 1L, nchar(token)-1L))
      if (op %in% c("M", "=", "X", "D", "N") && position1 >= ref && position1 < ref+n) {
        if (op %in% c("D", "N")) {
          key <- if (op == "D") "deleted_D" else "skipped_N"
          excluded[[key]] <- excluded[[key]] + 1L
        } else {
          q <- query + position1-ref
          if (fields[[10L]] == "*" || nchar(fields[[10L]]) < q) {
            excluded[["missing_sequence"]] <- excluded[["missing_sequence"]] + 1L
          } else if (fields[[11L]] == "*" || nchar(fields[[11L]]) < q) {
            excluded[["missing_quality"]] <- excluded[["missing_quality"]] + 1L
          } else if (utf8ToInt(substr(fields[[11L]], q, q))-33L < min_baseq) {
            excluded[["low_baseq"]] <- excluded[["low_baseq"]] + 1L
          } else {
            base <- toupper(substr(fields[[10L]], q, q))
            if (!base %in% names(counts)) base <- "N"
            counts[[base]] <- counts[[base]] + 1L
          }
        }
        break
      }
      if (op %in% c("M", "=", "X", "D", "N")) ref <- ref+n
      if (op %in% c("M", "=", "X", "I", "S")) query <- query+n
    }
  }
  answer <- data.frame(base = names(counts), count = unname(counts))
  attr(answer, "excluded") <- excluded
  answer
}

analyze_bam <- function(cfg) {
  cfg <- validate_config(cfg)
  tools <- Sys.which(c("samtools", "regtools"))
  if (any(!nzchar(tools))) stop("Missing executables in PATH: ", paste(names(tools)[!nzchar(tools)], collapse = ", "))
  work <- tempfile("regtools_shiny_"); dir.create(work, mode = "0700")
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  logs <- character(); warnings <- character()
  run <- function(tool, args, stdout = "|") {
    # Argument-vector invocation; never interpolate inputs into a shell command.
    p <- processx::run(unname(tools[[tool]]), as.character(args), stdout = stdout,
                      stderr = "|", timeout = 600, error_on_status = FALSE,
                      cleanup_tree = TRUE)
    logs <<- c(logs, paste(c(tool, vapply(as.character(args), shQuote, character(1))), collapse = " "),
               p$stderr)
    if (isTRUE(p$timeout) || is.na(p$status) || p$status != 0L)
      stop(tool, " failed: ", p$stderr, call. = FALSE)
    p$stdout
  }
  help <- processx::run(unname(tools[["regtools"]]), c("junctions", "extract", "-h"),
                       error_on_status = FALSE, timeout = 20)
  help_text <- paste(help$stdout, help$stderr)
  strand_arg <- if (grepl("XS,|XS tags|XS.*RF", help_text)) cfg$strand_mode else
    c(XS = "0", RF = "1", FR = "2")[[cfg$strand_mode]]
  ver_sam <- run("samtools", "--version")
  ver_reg <- processx::run(unname(tools[["regtools"]]), "--version",
                          error_on_status = FALSE, timeout = 20)
  reg_version <- paste(ver_reg$stdout, ver_reg$stderr)
  if (cfg$demo) {
    sam <- file.path(work, "demo.sam"); make_demo_sam(sam)
    bam <- file.path(work, "demo.bam")
    run("samtools", c("sort", "-o", bam, sam)); run("samtools", c("index", bam))
  } else {
    if (!file.exists(cfg$bam) || !grepl("\\.bam$", cfg$bam, ignore.case = TRUE))
      stop("Select an existing BAM file.")
    bam <- normalizePath(cfg$bam, mustWork = TRUE)
    idx <- c(paste0(bam, ".bai"), sub("\\.bam$", ".bai", bam, ignore.case = TRUE),
             paste0(bam, ".csi"), sub("\\.bam$", ".csi", bam, ignore.case = TRUE))
    if (!any(file.exists(idx))) stop("No BAI or CSI index was found next to the BAM. A matching index is required.")
    if (any(file.info(idx[file.exists(idx)])$mtime < file.info(bam)$mtime))
      warnings <- c(warnings, "The index is older than the BAM. Verify that the index matches the BAM.")
  }
  bam_stat_before <- file.info(bam)[, c("size", "mtime"), drop = FALSE]
  run("samtools", c("quickcheck", "-v", bam))
  hdr <- strsplit(run("samtools", c("view", "-H", bam)), "\n", fixed = TRUE)[[1L]]
  sq <- strsplit(hdr[startsWith(hdr, "@SQ\t")], "\t", fixed = TRUE)
  names_seq <- vapply(sq, function(x) sub("^SN:", "", x[startsWith(x, "SN:")][1L]), character(1))
  sizes_seq <- vapply(sq, function(x) as.numeric(sub("^LN:", "", x[startsWith(x, "LN:")][1L])), numeric(1))
  hit <- match(cfg$chrom, names_seq)
  if (is.na(hit)) stop("The selected BAM has no chromosome named ", cfg$chrom, ". Select a chromosome from this BAM and check the reference assembly.")
  if (cfg$end1 > sizes_seq[[hit]]) stop("The analysis end exceeds the chromosome length.")
  region <- paste0(cfg$chrom, ":", cfg$start1, "-", cfg$end1)
  flag_mask <- 4L + 256L + 512L + 2048L + if (cfg$exclude_duplicates) 1024L else 0L
  opts <- c("-q", cfg$mapq, "-F", flag_mask)
  fetched_count <- as.numeric(trimws(run("samtools", c("view", "-c", opts, bam, region))))
  if (!is.finite(fetched_count) || fetched_count > 100000)
    stop("The interval exceeds 100,000 candidate alignments. Select a smaller interval; reads will not be subsampled.")
  selected_sam <- file.path(work, "selected.sam")
  run("samtools", c("view", "-h", opts, "-o", selected_sam, bam, region))
  if (file.info(selected_sam)$size > 256 * 1024^2)
    stop("The regional SAM exceeds 256 MiB. Select a smaller interval.")
  text <- readLines(selected_sam, warn = FALSE)
  headers <- text[startsWith(text, "@")]; records <- text[!startsWith(text, "@")]
  parsed <- parse_sam_records(records, cfg)
  rna_bases <- if (!is.null(cfg$variant_position)) {
    if (cfg$variant_position < cfg$start1 || cfg$variant_position > cfg$end1)
      stop("Variant position must be inside the analysis interval.")
    count_rna_bases(records[parsed$keep], cfg$variant_position, cfg$baseq)
  } else NULL
  retained_sam <- file.path(work, "retained.sam")
  writeLines(c(headers, records[parsed$keep]), retained_sam)
  local_bam <- file.path(work, "retained.bam")
  run("samtools", c("sort", "-o", local_bam, retained_sam))
  run("samtools", c("index", local_bam))
  bed_path <- file.path(work, "junctions.bed")
  run("regtools", c("junctions", "extract", "-s", strand_arg,
                     "-a", cfg$anchor, "-m", cfg$min_intron, "-M", cfg$max_intron,
                     "-r", region, "-o", bed_path, local_bam))
  j <- read_regtools_bed(bed_path)
  # Full alignment records can have junctions outside the requested ROI.
  j <- j[j$chrom == cfg$chrom & j$intron_start0 >= cfg$start1 - 1L &
           j$intron_end0 <= cfg$end1, , drop = FALSE]
  e <- parsed$events
  e$key <- junction_key(e$chrom, e$intron_start0, e$intron_end0, e$strand)
  e <- e[e$key %in% j$key, , drop = FALSE]
  # Do not impose an extra per-read anchor filter: that changes native RegTools semantics.
  observed <- table(e$key)
  observed_counts <- as.integer(observed[match(j$key, names(observed))])
  observed_counts[is.na(observed_counts)] <- 0L
  if (any(observed_counts != j$score))
    stop("Audit failed: read-level CIGAR counts differ from native RegTools scores. Check tool versions, strandedness and unusual CIGAR strings. Unverified metrics are not reported.")
  e$intron_start1 <- e$intron_start0 + 1L; e$intron_end1 <- e$intron_end0
  depth_path <- file.path(work, "depth.tsv")
  # Input is already filtered. Re-include DUP to respect the user's duplicate choice.
  # No -s: overlapping paired ends are two reads; no -J: D is not covered sequence.
  run("samtools", c("depth", "-a", "-g", "1024", "-q", cfg$baseq,
                     "-Q", "0", "-r", region, local_bam), stdout = depth_path)
  depth <- data.frame(chrom = cfg$chrom, pos1 = seq.int(cfg$start1, cfg$end1), depth = 0L)
  if (file.info(depth_path)$size > 0) {
    d <- data.table::fread(depth_path, header = FALSE, data.table = FALSE)
    names(d) <- c("chrom", "pos1", "depth")
    m <- match(d$pos1, depth$pos1); ok <- !is.na(m) & d$chrom == cfg$chrom
    depth$depth[m[ok]] <- d$depth[ok]
  }
  if (!nrow(j)) warnings <- c(warnings, "No junctions passed the current filters. This does not establish that the region has no splicing.")
  if (any(parsed$reads$mapq == 255L)) warnings <- c(warnings,
    "Some retained records have MAPQ=255, meaning mapping quality is unavailable. Do not automatically interpret these as uniquely mapped reads.")
  if (any(parsed$reads$strand == "?")) warnings <- c(warnings,
    "Some reads have unknown strand (?). Metrics combine all strands; a missing XS tag is not interpreted as positive strand.")
  if (cfg$strand_mode != "XS") warnings <- c(warnings,
    "RF/FR follows RegTools flag rules. Verify library direction using a known control, especially for single-end data.")
  if (cfg$nh1_only && parsed$missing_nh > 0) warnings <- c(warnings,
    paste0("Requiring NH==1 also excluded ", parsed$missing_nh, " candidate records with missing NH."))
  if (!identical(bam_stat_before, file.info(bam)[, c("size", "mtime"), drop = FALSE]))
    stop("The BAM size or modification time changed during analysis. Retry after the file is stable.")
  cfg$bam <- if (cfg$demo) "synthetic_demo" else bam
  list(config = cfg, reads = parsed$reads, events = e, junctions = j, depth = depth, rna_bases = rna_bases,
       native_audit_passed = TRUE, candidate_alignments = fetched_count,
       span_only_excluded = parsed$span_only_excluded,
       warnings = unique(warnings), log = logs,
       versions = list(samtools = ver_sam, regtools = reg_version,
                       regtools_extract_help = help_text,
                       R = R.version.string, processx = as.character(utils::packageVersion("processx")),
                       data_table = as.character(utils::packageVersion("data.table"))),
       source_metadata = list(size_bytes = bam_stat_before$size,
                              mtime = as.character(bam_stat_before$mtime)),
       completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
}

summarize_target <- function(result, start1, end1, delta = 5L, flank = 50L) {
  s <- int_scalar(start1, "Target intron start", 1)
  t <- int_scalar(end1, "Target intron end", s)
  delta <- int_scalar(delta, "Near tolerance", 0, 1000)
  c <- result$config
  if (s - delta <= c$start1 || t + delta >= c$end1)
    stop("The target and near tolerance must leave at least 1 bp on each side inside the loaded interval. Expand the analysis interval and rerun.")
  e <- result$events
  exact <- e$intron_start1 == s & e$intron_end1 == t
  near <- abs(e$intron_start1 - s) <= delta & abs(e$intron_end1 - t) <= delta & !exact
  n <- nrow(result$reads)
  exact_n <- length(unique(e$read_id[exact]))
  near_n <- length(unique(e$read_id[near]))
  depth_mean <- function(a, b) {
    v <- result$depth$depth[result$depth$pos1 >= a & result$depth$pos1 <= b]
    if (length(v)) mean(v) else NA_real_
  }
  data.frame(
    chrom = c$chrom, target_intron_start1 = s, target_intron_end1 = t,
    near_delta_bp = delta, denominator_start1 = c$start1, denominator_end1 = c$end1,
    count_unit = "primary_alignment_record; paired_ends_separate",
    strand_scope = "all_strands_including_unknown",
    denominator_reads = n,
    junction_reads = length(unique(e$read_id)),
    mean_read_depth = mean(result$depth$depth),
    left_flank_mean_depth = depth_mean(max(c$start1, s - flank), s - 1L),
    right_flank_mean_depth = depth_mean(t + 1L, min(c$end1, t + flank)),
    junction_exact_count = exact_n,
    junction_near_count = near_n,
    exact_or_near_reads = length(unique(e$read_id[exact | near])),
    junction_per_100_reads = if (n) 100 * exact_n / n else NA_real_,
    stringsAsFactors = FALSE)
}

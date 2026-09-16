# RegTools Shiny MVP — computational backend.
# Coordinates: GUI = 1-based closed intron interval; internal = 0-based half-open.
# Counting unit = retained primary alignment record (paired ends count separately).
# Requires processx and data.table; external programs: samtools and regtools.

int_scalar <- function(x, name, lower = 0, upper = 2e9) {
  if (length(x) != 1L || is.na(x) || !is.numeric(x) ||
      !is.finite(x) || x != floor(x) || x < lower || x > upper)
    stop(name, " 必须是 ", lower, " 到 ", upper, " 之间的整数。", call. = FALSE)
  as.integer(x)
}

validate_config <- function(cfg) {
  cfg$start1 <- int_scalar(cfg$start1, "区间起点", 1)
  cfg$end1 <- int_scalar(cfg$end1, "区间终点", cfg$start1)
  if (cfg$end1 - cfg$start1 + 1 > 250000)
    stop("原型限单次区间长度 250,000 bp；请缩小区间。", call. = FALSE)
  if (length(cfg$chrom) != 1L || is.na(cfg$chrom) ||
      !grepl("^[A-Za-z0-9_.-]+$", cfg$chrom))
    stop("染色体名称只支持字母、数字、点、下划线、短横线。", call. = FALSE)
  cfg$mapq <- int_scalar(cfg$mapq, "MAPQ", 0, 255)
  cfg$baseq <- int_scalar(cfg$baseq, "BaseQ", 0, 93)
  cfg$anchor <- int_scalar(cfg$anchor, "Anchor", 1, 1000)
  cfg$min_intron <- int_scalar(cfg$min_intron, "最短内含子", 1)
  cfg$max_intron <- int_scalar(cfg$max_intron, "最长内含子", cfg$min_intron)
  if (!cfg$strand_mode %in% c("XS", "RF", "FR")) stop("无效的链模式。")
  cfg$exclude_duplicates <- isTRUE(cfg$exclude_duplicates)
  cfg$nh1_only <- isTRUE(cfg$nh1_only)
  cfg$demo <- isTRUE(cfg$demo)
  cfg
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
    stop("不支持或无效的 CIGAR: ", cigar, call. = FALSE)
  lens <- as.integer(sub("[MIDNSHP=X]$", "", tokens))
  if (anyNA(lens) || any(lens <= 0)) stop("CIGAR 长度无效。")
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
  if (length(z) && any(lengths(z) < 11L)) stop("SAM 行少于 11 列。")
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
  if (ncol(b) != 12L) stop("RegTools 输出不是 BED12。")
  names(b) <- names(empty)[seq_len(12L)]
  for (n in c("chromStart", "chromEnd", "score", "thickStart", "thickEnd", "blockCount"))
    b[[n]] <- as.integer(b[[n]])
  if (anyNA(b$score) || any(b$score < 0) || any(b$blockCount != 2L))
    stop("RegTools BED12 的 score 或 blockCount 无效。")
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

analyze_bam <- function(cfg) {
  cfg <- validate_config(cfg)
  tools <- Sys.which(c("samtools", "regtools"))
  if (any(!nzchar(tools))) stop("PATH 中缺少：", paste(names(tools)[!nzchar(tools)], collapse = ", "))
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
      stop(tool, " 执行失败：", p$stderr, call. = FALSE)
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
      stop("请选择存在的 BAM 文件。")
    bam <- normalizePath(cfg$bam, mustWork = TRUE)
    idx <- c(paste0(bam, ".bai"), sub("\\.bam$", ".bai", bam, ignore.case = TRUE),
             paste0(bam, ".csi"), sub("\\.bam$", ".csi", bam, ignore.case = TRUE))
    if (!any(file.exists(idx))) stop("BAM 旁边没有 BAI/CSI；请先排序并建立索引。")
    if (any(file.info(idx[file.exists(idx)])$mtime < file.info(bam)$mtime))
      warnings <- c(warnings, "索引修改时间早于 BAM：请核实索引与 BAM 是否匹配。")
  }
  bam_stat_before <- file.info(bam)[, c("size", "mtime"), drop = FALSE]
  run("samtools", c("quickcheck", "-v", bam))
  hdr <- strsplit(run("samtools", c("view", "-H", bam)), "\n", fixed = TRUE)[[1L]]
  sq <- strsplit(hdr[startsWith(hdr, "@SQ\t")], "\t", fixed = TRUE)
  names_seq <- vapply(sq, function(x) sub("^SN:", "", x[startsWith(x, "SN:")][1L]), character(1))
  sizes_seq <- vapply(sq, function(x) as.numeric(sub("^LN:", "", x[startsWith(x, "LN:")][1L])), numeric(1))
  hit <- match(cfg$chrom, names_seq)
  if (is.na(hit)) stop("BAM 没有染色体 ", cfg$chrom, "；检查 chr 前缀和参考基因组。")
  if (cfg$end1 > sizes_seq[[hit]]) stop("分析终点超过染色体长度。")
  region <- paste0(cfg$chrom, ":", cfg$start1, "-", cfg$end1)
  flag_mask <- 4L + 256L + 512L + 2048L + if (cfg$exclude_duplicates) 1024L else 0L
  opts <- c("-q", cfg$mapq, "-F", flag_mask)
  fetched_count <- as.numeric(trimws(run("samtools", c("view", "-c", opts, bam, region))))
  if (!is.finite(fetched_count) || fetched_count > 100000)
    stop("区间候选 alignment 超过原型上限 100,000；请缩小区间，不会偷偷抽样。")
  selected_sam <- file.path(work, "selected.sam")
  run("samtools", c("view", "-h", opts, "-o", selected_sam, bam, region))
  if (file.info(selected_sam)$size > 256 * 1024^2)
    stop("局部 SAM 超过 256 MiB；请缩小区间。")
  text <- readLines(selected_sam, warn = FALSE)
  headers <- text[startsWith(text, "@")]; records <- text[!startsWith(text, "@")]
  parsed <- parse_sam_records(records, cfg)
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
    stop("核验失败：逐 read CIGAR 计数与原生 RegTools score 不一致。请检查工具版本、链模式和异常 CIGAR；不输出未核实指标。")
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
  if (!nrow(j)) warnings <- c(warnings, "没有检出通过当前过滤的 junction；这不等于证明该区域不存在剪接。")
  if (any(parsed$reads$mapq == 255L)) warnings <- c(warnings,
    "保留了 MAPQ=255 的记录；255 表示 mapping quality 不可用，不应自动等同于唯一比对。")
  if (any(parsed$reads$strand == "?")) warnings <- c(warnings,
    "部分 reads 链方向未知（?）。原型汇总所有链；不把缺失 XS 解释为正链。")
  if (cfg$strand_mode != "XS") warnings <- c(warnings,
    "RF/FR 按 RegTools flags 规则推断；请用已知方向的对照核实你的文库，尤其是单端数据。")
  if (cfg$nh1_only && parsed$missing_nh > 0) warnings <- c(warnings,
    paste0("启用 NH==1 后，缺失 NH 的 ", parsed$missing_nh, " 条候选记录也被排除。"))
  if (!identical(bam_stat_before, file.info(bam)[, c("size", "mtime"), drop = FALSE]))
    stop("分析期间 BAM 大小或修改时间发生变化；请在文件稳定后重新运行。")
  cfg$bam <- if (cfg$demo) "synthetic_demo" else bam
  list(config = cfg, reads = parsed$reads, events = e, junctions = j, depth = depth,
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
  s <- int_scalar(start1, "目标内含子起点", 1)
  t <- int_scalar(end1, "目标内含子终点", s)
  delta <- int_scalar(delta, "near 容差", 0, 1000)
  c <- result$config
  if (s - delta <= c$start1 || t + delta >= c$end1)
    stop("目标及 near 容差两侧必须各留至少 1 bp 在已读取区间内；请扩大区间后重新读取。")
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

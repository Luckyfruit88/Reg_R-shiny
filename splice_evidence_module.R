# Read-only splice evidence view of a saved WGS result. No BAM is opened here.

splice_called_groups <- function(model) {
  g <- model$genotypes
  if (is.null(g) || !nrow(g)) return(data.frame())
  if ("call_status" %in% names(g)) g <- g[g$call_status == "CALLED", , drop = FALSE]
  g
}

splice_primary_site <- function(model) {
  s <- model$sites
  if (is.null(s) || !nrow(s) || !"variant_overlaps" %in% names(s)) return(NULL)
  s <- s[!is.na(s$variant_overlaps) & s$variant_overlaps, , drop = FALSE]
  if (!nrow(s)) return(NULL)
  focus <- model$metadata$focus_candidate_id
  if ("focus_candidate_id" %in% names(model$metadata)) {
    if (length(focus) != 1L || is.na(focus) || !focus %in% s$candidate_id) return(NULL)
    rows <- s[s$candidate_id == focus, , drop = FALSE]
    allele <- model$metadata$focus_allele_index
    if (length(allele) == 1L && !is.na(allele) && allele %in% rows$allele_index)
      return(rows[match(allele, rows$allele_index), , drop = FALSE])
    return(rows[1L, , drop = FALSE])
  }
  source <- if ("source" %in% names(s)) as.character(s$source) else rep("", nrow(s))
  priority <- ifelse(grepl("RNA|observed", source, ignore.case = TRUE), 0L,
    ifelse(grepl("GENCODE|annotation", source, ignore.case = TRUE), 1L, 2L))
  distance <- if ("distance_to_variant" %in% names(s)) s$distance_to_variant else abs(s$anchor1 - model$window$variant_pos1)
  known <- !is.na(s$strand) & s$strand %in% c("+", "-")
  evaluable <- known & !is.na(s$ref_motif) & !is.na(s$alt_motif) & !grepl("^UNAVAILABLE", s$change)
  s[order(!evaluable, !known, priority, distance, s$start1, s$site_type, s$strand)[[1L]], , drop = FALSE]
}

splice_alt_label <- function(model, site) {
  if (is.null(site) || !"allele_index" %in% names(site)) return("ALT")
  index <- site$allele_index[[1L]]
  alt <- model$metadata$variant$alt
  letters <- if (length(alt) == 1L && !is.na(alt)) strsplit(as.character(alt), ",", fixed = TRUE)[[1L]] else character()
  paste0("ALT ", index, if (is.finite(index) && index >= 1L && index <= length(letters)) paste0("=", letters[[index]]) else "")
}

splice_near_junctions <- function(model) {
  j <- model$junction_support
  if (is.null(j) || !nrow(j)) return(data.frame())
  g <- splice_called_groups(model)
  j <- j[j$genotype %in% g$genotype, , drop = FALSE]
  w <- model$window
  inside <- (j$intron_start1 >= w$start1 & j$intron_start1 <= w$end1) |
    (j$intron_end1 >= w$start1 & j$intron_end1 <= w$end1)
  j <- j[inside, , drop = FALSE]
  if (!"junction_id" %in% names(j)) j$junction_id <- paste(j$chrom, j$intron_start1, j$intron_end1, j$strand, sep = ":")
  j$dna_alleles <- g$dna_alleles[match(j$genotype, g$genotype)]
  j
}

splice_rna_boundaries <- function(model) {
  j <- splice_near_junctions(model)
  empty <- data.frame(pos1 = numeric(), role = character(), strand = character(), delta = numeric())
  if (!nrow(j)) return(empty)
  j <- j[is.finite(j$count_sum) & j$count_sum > 0, , drop = FALSE]
  if (!nrow(j)) return(empty)
  w <- model$window
  left <- data.frame(pos1 = j$intron_start1,
    role = ifelse(j$strand == "+", "donor", ifelse(j$strand == "-", "acceptor", "boundary")), strand = j$strand)
  right <- data.frame(pos1 = j$intron_end1,
    role = ifelse(j$strand == "+", "acceptor", ifelse(j$strand == "-", "donor", "boundary")), strand = j$strand)
  b <- unique(rbind(left, right))
  b <- b[b$pos1 >= w$start1 & b$pos1 <= w$end1, , drop = FALSE]
  b$delta <- b$pos1 - w$variant_pos1
  b[order(abs(b$delta), b$pos1, b$role, b$strand), , drop = FALSE]
}

splice_exon_regions <- function(annotation, window) {
  empty <- data.frame(start1 = numeric(), end1 = numeric())
  if (is.null(annotation) || inherits(annotation, "error") || is.null(annotation$exons) || !nrow(annotation$exons)) return(empty)
  e <- annotation$exons
  e <- e[e$chrom == window$chrom & e$start1 <= window$end1 & e$end1 >= window$start1, , drop = FALSE]
  if (!nrow(e)) return(empty)
  x <- data.frame(start1 = pmax(e$start1, window$start1), end1 = pmin(e$end1, window$end1))
  x <- unique(x[order(x$start1, x$end1), , drop = FALSE])
  merged <- list(); start <- x$start1[[1L]]; end <- x$end1[[1L]]
  for (i in seq_len(nrow(x))[-1L]) {
    if (x$start1[[i]] <= end + 1) end <- max(end, x$end1[[i]]) else {
      merged[[length(merged) + 1L]] <- c(start, end)
      start <- x$start1[[i]]; end <- x$end1[[i]]
    }
  }
  merged[[length(merged) + 1L]] <- c(start, end)
  z <- as.data.frame(do.call(rbind, merged)); names(z) <- c("start1", "end1")
  z
}

splice_motif_label <- function(model, genotype) {
  if (!identical(model$metadata$reference_status, "REF_MATCH")) return("unavailable")
  primary <- splice_primary_site(model)
  evidence <- model$genotype_sites
  if (!is.null(primary) && !is.null(evidence) && nrow(evidence) && "motif_alleles" %in% names(evidence)) {
    hit <- evidence[evidence$genotype == genotype & evidence$candidate_id == primary$candidate_id[[1L]], , drop = FALSE]
    value <- unique(as.character(hit$motif_alleles))
    value <- value[!is.na(value) & nzchar(value)]
    if (length(value)) return(paste(value, collapse = "; "))
  }
  g <- model$genotypes
  i <- match(genotype, g$genotype)
  if (is.na(i)) return("unavailable")
  if ("motif_summary" %in% names(g) && length(g$motif_summary[[i]]) == 1L &&
      !is.na(g$motif_summary[[i]]) && nzchar(g$motif_summary[[i]])) return(as.character(g$motif_summary[[i]]))
  "unavailable"
}

# Public plotting function reused by the UI, downloads and independent checks.
plot_splice_coverage <- function(model, annotation = NULL) {
  w <- model$window; g <- splice_called_groups(model)
  old <- graphics::par(mar = c(4.8, 9.5, 4.2, 3.2))
  on.exit(graphics::par(old))
  if (!nrow(g)) { graphics::plot.new(); graphics::text(.5, .5, "No called DNA genotype groups in this saved result"); return(invisible(NULL)) }
  positions <- seq.int(w$start1, w$end1)
  d <- model$coverage
  observed <- d$mean_depth[is.finite(d$mean_depth) & d$genotype %in% g$genotype]
  ymax <- if (length(observed) && max(observed) > 0) max(observed) else 1
  n <- nrow(g); top <- 2.0 + (n - 1L) * 1.3 + 1.0
  graphics::plot(NA, xlim = c(w$start1 - .5, w$end1 + .5), ylim = c(0, top), xaxs = "i", yaxs = "i",
    xlab = paste0(w$chrom, " · genomic position (1-based)"), ylab = "", yaxt = "n",
    main = "RNA coverage at base resolution · saved genotype groups")
  graphics::mtext(paste0("One shared height scale: 0-", format(round(ymax, 3), trim = TRUE),
    " mean RNA depth · each tile is one genomic base"), side = 3, line = .5, cex = .8)
  colors <- grDevices::hcl.colors(max(3L, n), "Dark 3")[seq_len(n)]
  exon <- splice_exon_regions(annotation, w)
  for (i in seq_len(n)) {
    baseline <- 2 + (n - i) * 1.3
    if (nrow(exon)) graphics::rect(exon$start1 - .5, baseline, exon$end1 + .5, baseline + .84,
      col = "#EEEEEE", border = NA)
    row <- d[d$genotype == g$genotype[[i]], , drop = FALSE]
    values <- row$mean_depth[match(positions, row$pos1)]
    usable <- is.finite(values)
    dna <- if (!is.na(g$dna_alleles[[i]]) && nzchar(g$dna_alleles[[i]])) g$dna_alleles[[i]] else "DNA alleles unavailable"
    motif <- splice_motif_label(model, g$genotype[[i]])
    if (nchar(motif) > 46L) motif <- paste0(substr(motif, 1L, 43L), "...")
    label <- paste0(dna, " (", g$genotype[[i]], ")\nDNA motif: ", motif,
      "\nsuccessful RNA N = ", g$analyzed_n[[i]])
    graphics::text(w$start1 - .5 - length(positions) * .018, baseline + .43,
      label, adj = c(1, .5), cex = .76, xpd = NA, col = colors[[i]])
    if (any(usable)) {
      graphics::segments(w$start1 - .5, baseline, w$end1 + .5, baseline, col = "#C3C3C3")
      graphics::rect(positions[usable] - .44, baseline, positions[usable] + .44,
        baseline + .82 * values[usable] / ymax, col = colors[[i]], border = NA)
      if (all(values[usable] == 0) && all(usable)) graphics::text(mean(c(w$start1, w$end1)), baseline + .4,
        "Measured zero depth across this window", col = colors[[i]], cex = .77)
      if (any(!usable)) graphics::points(positions[!usable], rep(baseline + .4, sum(!usable)),
        pch = 4, cex = .45, col = "#929292")
    } else graphics::text(mean(c(w$start1, w$end1)), baseline + .4,
      "RNA depth unavailable · missing is not zero", col = "#727272", cex = .8)
  }
  base <- model$bases
  letters <- if (!is.null(base) && nrow(base)) as.character(base$ref_base[match(positions, base$pos1)]) else rep(NA_character_, length(positions))
  palette <- c(A = "#207B52", C = "#326EA8", G = "#A87617", T = "#AF3D47", N = "#666666")
  valid <- !is.na(letters) & nzchar(letters)
  graphics::text(w$start1 - .5 - length(positions) * .018, .42, "Reference (+)", adj = 1, cex = .78, xpd = NA)
  if (any(valid)) {
    graphics::rect(positions[valid] - .46, .16, positions[valid] + .46, .67, col = "#F7F7F7", border = "#DFDFDF")
    base_colors <- unname(palette[letters[valid]]); base_colors[is.na(base_colors)] <- "#666666"
    cex <- min(.95, max(.28, 48 / length(positions)))
    graphics::text(positions[valid], .42, letters[valid], col = base_colors, family = "mono", cex = cex)
  } else graphics::text(mean(c(w$start1, w$end1)), .42, "Reference sequence unavailable", col = "#777777", cex = .85)
  primary <- splice_primary_site(model)
  if (!is.null(primary) && isTRUE(primary$variant_overlaps[[1L]])) {
    graphics::rect(max(w$start1 - .5, primary$start1[[1L]] - .5), .12,
      min(w$end1 + .5, primary$end1[[1L]] + .5), .71, border = "#A66100", lwd = 2)
    motif_known <- identical(model$metadata$reference_status, "REF_MATCH") &&
      !is.na(primary$ref_motif[[1L]]) && !is.na(primary$alt_motif[[1L]])
    sequence_only <- "sequence_only" %in% names(primary) && isTRUE(primary$sequence_only[[1L]])
    motif_text <- paste0(if (sequence_only) "Sequence-only motif candidate: " else "", primary$site_type[[1L]], " (", primary$strand[[1L]], "), ", splice_alt_label(model, primary), ": ",
      if (motif_known) paste0(primary$ref_motif[[1L]], " -> ", primary$alt_motif[[1L]]) else "DNA motif unavailable")
    graphics::text(mean(c(primary$start1[[1L]], primary$end1[[1L]])), 1.52,
      motif_text, cex = .77, col = "#8D5600", xpd = FALSE)
  }
  boundaries <- head(splice_rna_boundaries(model), 6L)
  if (nrow(boundaries)) for (i in seq_len(nrow(boundaries))) {
    b <- boundaries[i, , drop = FALSE]
    graphics::abline(v = b$pos1, lty = 3, col = "#596A7A", lwd = .8)
    delta <- if (b$delta > 0) paste0("+", b$delta) else as.character(b$delta)
    graphics::text(b$pos1, .93 + (i %% 2L) * .22,
      paste0(if (b$role == "acceptor") "A" else if (b$role == "donor") "D" else "B", " ", delta),
      cex = .66, col = "#465B70")
  }
  graphics::abline(v = w$variant_pos1, lty = 2, col = "#B13F32", lwd = 1.5)
  invisible(list(ymax = ymax, genotypes = g$genotype, positions = positions,
    boundary_positions = boundaries$pos1, primary_site = primary))
}

splice_junction_display <- function(model, max_junctions = 12L) {
  j <- splice_near_junctions(model)
  if (!nrow(j)) return(j)
  ids <- unique(j$junction_id)
  distance <- vapply(ids, function(id) {
    r <- j[j$junction_id == id, , drop = FALSE]
    min(abs(c(r$intron_start1, r$intron_end1) - model$window$variant_pos1))
  }, numeric(1))
  support <- vapply(ids, function(id) {
    x <- j$count_sum[j$junction_id == id]; if (any(is.finite(x))) max(x[is.finite(x)]) else -Inf
  }, numeric(1))
  kept <- head(ids[order(distance, -support, ids)], max_junctions)
  j[j$junction_id %in% kept, , drop = FALSE][order(match(j$junction_id[j$junction_id %in% kept], kept)), , drop = FALSE]
}

# A junction occupies one shared matrix row in every genotype; absent cells stay NA.
plot_splice_junction_support <- function(model, max_junctions = 12L) {
  j <- splice_junction_display(model, max_junctions)
  g <- splice_called_groups(model)
  old <- graphics::par(mar = c(3.1, 14.5, 4.8, 1.3))
  on.exit(graphics::par(old))
  if (!nrow(j) || !nrow(g)) {
    graphics::plot.new(); graphics::text(.5, .5, "No saved RNA junction endpoint in this display window")
    return(invisible(list(junction_ids = character())))
  }
  ids <- unique(j$junction_id)
  maximum <- if (any(is.finite(j$count_mean))) max(0, j$count_mean[is.finite(j$count_mean)]) else 0
  colors <- grDevices::colorRampPalette(c("#FFFFFF", "#D9E7F2", "#367CA9", "#123D62"))(100L)
  graphics::plot(NA, xlim = c(.5, nrow(g) + .5), ylim = c(.5, length(ids) + .5),
    xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "",
    main = "Observed RNA junction support · one shared scale")
  graphics::axis(3, at = seq_len(nrow(g)), labels = g$dna_alleles, tick = FALSE, line = -.2, cex.axis = .9)
  for (i in seq_along(ids)) {
    rows <- j[j$junction_id == ids[[i]], , drop = FALSE]
    y <- length(ids) - i + 1L
    start <- rows$intron_start1[[1L]]; end <- rows$intron_end1[[1L]]; strand <- rows$strand[[1L]]
    near <- if (abs(start - model$window$variant_pos1) <= abs(end - model$window$variant_pos1)) start else end
    delta <- near - model$window$variant_pos1
    label <- paste0(start, "-", end, " (", strand, ")\nnearest endpoint ", if (delta > 0) "+" else "", delta, " bp")
    graphics::text(.43, y, label, adj = c(1, .5), cex = .7, xpd = NA)
    for (k in seq_len(nrow(g))) {
      r <- rows[rows$genotype == g$genotype[[k]], , drop = FALSE]
      known <- nrow(r) == 1L && is.finite(r$count_mean[[1L]]) && is.finite(r$count_sum[[1L]]) &&
        is.finite(r$analyzed_n[[1L]]) && r$analyzed_n[[1L]] > 0
      fill <- if (known) colors[[if (maximum > 0) max(1L, min(100L, 1L + round(99 * r$count_mean[[1L]] / maximum))) else 1L]] else "#ECECEC"
      graphics::rect(k - .47, y - .44, k + .47, y + .44, col = fill, border = "#B9C2CB")
      if (known) {
        label <- paste0(format(round(r$count_mean[[1L]], 3), trim = TRUE), "\n", r$count_sum[[1L]], " / ", r$analyzed_n[[1L]])
        text_col <- if (maximum > 0 && r$count_mean[[1L]] / maximum > .62) "white" else "#233544"
      } else { label <- "NA\nunavailable"; text_col <- "#676767" }
      graphics::text(k, y, label, cex = .79, col = text_col)
    }
  }
  graphics::mtext(paste0("Each cell: mean support, then total / successful N · shared color range 0-", format(round(maximum, 3), trim = TRUE)),
    side = 1, line = 1.1, cex = .76)
  invisible(list(junction_ids = ids, max_mean = maximum, genotypes = g$genotype))
}

splice_evidence_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::selectInput(ns("flank"), "Base-resolution window around the loaded variant",
      c("±10 bp" = "10", "±25 bp" = "25", "±50 bp" = "50", "±100 bp" = "100"), selected = "25"),
    shiny::uiOutput(ns("status")),
    shiny::uiOutput(ns("coverage_plot_ui")),
    shiny::p("Tiles show mean RNA depth across successfully analyzed samples in each DNA genotype group, with one shared height scale. They are not individual reads or reconstructed RNA haplotypes. Gray background marks an exon in at least one loaded GENCODE model; × marks missing depth."),
    shiny::uiOutput(ns("motif_note")),
    shiny::p("Reference (+) letters and coordinates follow the genomic positive strand. GT/AG motifs follow the stated transcript strand; negative-strand motifs are reverse complemented. DNA genotype labels retain the VCF genomic alleles. Motif substitution applies only the selected variant, not neighboring variants or a phased haplotype; it is separate from measured RNA support and does not establish causality."),
    shiny::uiOutput(ns("junction_plot_ui")),
    shiny::uiOutput(ns("junction_note")),
    shiny::tags$details(shiny::tags$summary("Detailed nearby junction counts"), DT::DTOutput(ns("junction_table"))),
    shiny::tags$details(shiny::tags$summary("All nearby splice-site predictions and provenance"),
      DT::DTOutput(ns("sites_table")), shiny::verbatimTextOutput(ns("metadata"))),
    shiny::downloadButton(ns("coverage_png"), "Coverage PNG"),
    shiny::downloadButton(ns("coverage_pdf"), "Coverage PDF"),
    shiny::downloadButton(ns("junction_csv"), "RNA junction evidence CSV"),
    shiny::downloadButton(ns("evidence_json"), "Sequence / motif evidence JSON"))
}

splice_evidence_server <- function(id, result, reference, model_source,
                                   annotation_resources = NULL, annotation_source = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    compact_result <- shiny::reactive({
      r <- result()
      cfg <- r$provenance$config
      if (is.null(cfg)) cfg <- r$config
      list(variant = r$variant, group_summary = r$group_summary, depth = r$depth,
        junctions = r$junctions, provenance = list(config = cfg, selection = r$provenance$selection))
    })
    context <- shiny::reactive({
      r <- compact_result(); cfg <- r$provenance$config
      list(chrom = cfg$chrom, start1 = cfg$start1, end1 = cfg$end1, build = r$variant$build[[1L]])
    })
    context_key <- shiny::reactive({
      x <- context()
      paste(x$build, x$chrom, x$start1, x$end1, sep = ":")
    })
    annotation_cache <- shiny::reactiveVal(NULL)
    annotation_task <- shiny::ExtendedTask$new(function(resources, interval, source, key) {
      promises::future_promise({
        value <- tryCatch({
          e <- new.env(parent = globalenv()); sys.source(source, envir = e)
          e$query_annotation(resources, interval$chrom, interval$start1, interval$end1, interval$build)
        }, error = function(e) e)
        list(key = key, value = value)
      }, seed = TRUE)
    })
    shiny::observeEvent(context_key(), {
      key <- context_key()
      old <- annotation_cache()
      if (!is.null(old) && identical(old$key, key)) return()
      if (is.null(annotation_resources) || inherits(annotation_resources, "error") || is.null(annotation_source)) {
        annotation_cache(list(key = key, value = if (inherits(annotation_resources, "error")) annotation_resources else NULL))
      } else annotation_task$invoke(annotation_resources, context(), annotation_source, key)
    })
    shiny::observeEvent(annotation_task$status(), {
      if (annotation_task$status() == "success") annotation_cache(annotation_task$result())
    }, ignoreInit = TRUE)
    flank <- shiny::reactive({
      x <- if (is.null(input$flank)) 25L else suppressWarnings(as.integer(input$flank))
      shiny::validate(shiny::need(length(x) == 1L && !is.na(x) && x %in% c(10L, 25L, 50L, 100L), "Choose a supported base-resolution window."))
      x
    })
    # Bind a displayed model to the actual saved evidence, not just the latest
    # worker invocation. This also hides old results while a new ROI is waiting
    # for its annotation query and no replacement model has been invoked yet.
    request_stamp <- shiny::reactive({
      r <- compact_result(); cfg <- r$provenance$config; v <- r$variant
      start <- max(1, cfg$start1, v$pos1[[1L]] - flank())
      end <- min(cfg$end1, v$pos1[[1L]] + flank())
      list(context_key = context_key(), flank = flank(), variant = v,
        group_summary = r$group_summary, config = cfg,
        depth = r$depth[r$depth$pos1 >= start & r$depth$pos1 <= end & r$depth$chrom == v$chrom[[1L]], , drop = FALSE],
        junctions = r$junctions)
    })
    serial <- shiny::reactiveVal(0L)
    model_task <- shiny::ExtendedTask$new(function(snapshot, annotation, fasta, half_window, source, request_serial, stamp) {
      promises::future_promise({
        e <- new.env(parent = globalenv()); sys.source(source, envir = e)
        list(model = e$build_splice_evidence(snapshot, annotation, fasta, flank = half_window),
          annotation = annotation, request_serial = request_serial, request_stamp = stamp)
      }, seed = TRUE)
    })
    shiny::observeEvent(list(compact_result(), flank(), annotation_cache()), {
      cache <- annotation_cache()
      shiny::req(!is.null(cache), identical(cache$key, context_key()))
      serial(serial() + 1L)
      model_task$invoke(compact_result(), cache$value, reference, flank(), model_source, serial(), request_stamp())
    })
    answer <- shiny::reactive({
      x <- model_task$result()
      cache <- annotation_cache()
      shiny::req(!is.null(cache), identical(cache$key, context_key()),
        identical(x$request_serial, serial()), identical(x$request_stamp, request_stamp()))
      x
    })
    model <- shiny::reactive(answer()$model)
    output$status <- shiny::renderUI({
      s <- model_task$status()
      if (s == "initial" || s == "running") return(shiny::p(if (s == "running" || annotation_task$status() == "running")
        "Loading reference and splice context for the saved RNA result..." else "Load a genotype comparison to view base-resolution splice evidence."))
      if (s == "error") {
        message <- tryCatch(model_task$result(), error = function(e) conditionMessage(e))
        return(shiny::div(class = "alert alert-warning", "Splice evidence unavailable: ", as.character(message), ". The original RNA allele table remains available."))
      }
      m <- model(); w <- m$window
      shiny::tagList(shiny::p("Loaded result window: ", w$chrom, ":", w$start1, "–", w$end1,
        " · selected variant ", w$variant_pos1, " · ", w$build,
        " · reference status: ", m$metadata$reference_status),
        if (length(m$warnings)) shiny::div(class = "alert alert-warning", paste(unique(m$warnings), collapse = "; ")))
    })
    output$coverage_plot_ui <- shiny::renderUI(shiny::plotOutput(session$ns("coverage_plot"),
      height = max(340L, 160L + 105L * nrow(splice_called_groups(model())))))
    output$coverage_plot <- shiny::renderPlot(plot_splice_coverage(model(), answer()$annotation))
    output$motif_note <- shiny::renderUI({
      m <- model(); s <- splice_primary_site(m); b <- splice_rna_boundaries(m)
      shiny::tagList(
        if (!is.null(s) && identical(m$metadata$reference_status, "REF_MATCH")) shiny::p(
          if ("sequence_only" %in% names(s) && isTRUE(s$sequence_only[[1L]])) "Sequence-only motif candidate: " else "Highlighted boundary: ",
          s$site_type[[1L]], " (", s$strand[[1L]], ") at ",
          s$start1[[1L]], "–", s$end1[[1L]], " · transcript-oriented motif ", s$ref_motif[[1L]], " → ", s$alt_motif[[1L]],
          " · ", splice_alt_label(m, s), " · ", s$change[[1L]], ". This is the indicated target-variant DNA substitution; all ALT rows remain in the table."),
        shiny::p("Red dashed line: selected variant. Dotted RNA boundary markers: A = acceptor, D = donor, B = unresolved strand; offsets are relative to the variant. ",
          min(6L, nrow(b)), " of ", nrow(b), " supported nearby boundaries are marked; complete evidence remains in the table."))
    })
    output$junction_plot_ui <- shiny::renderUI({
      n <- length(unique(splice_junction_display(model())$junction_id))
      shiny::plotOutput(session$ns("junction_plot"), height = max(230L, 120L + 62L * n))
    })
    output$junction_plot <- shiny::renderPlot(plot_splice_junction_support(model()))
    output$junction_note <- shiny::renderUI({
      j <- splice_near_junctions(model()); all_n <- length(unique(j$junction_id))
      shiny::p("Showing ", min(12L, all_n), " of ", all_n,
        " junctions with an endpoint in the displayed window, nearest to the selected variant first. Each junction has the same row across genotypes. Measured zero is retained; NA is unavailable. The drawing cap does not filter this full table. These are saved RNA junction counts, independent of the DNA motif prediction.")
    })
    output$junction_table <- DT::renderDT(splice_near_junctions(model()), rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE), selection = "none")
    output$sites_table <- DT::renderDT(model()$sites, rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE), selection = "none")
    output$metadata <- shiny::renderPrint(print(model()$metadata))
    output$coverage_png <- shiny::downloadHandler("splice_base_resolution.png", function(file) {
      grDevices::png(file, width = 1800, height = max(650, 280 + 180 * nrow(splice_called_groups(model()))), res = 160)
      on.exit(grDevices::dev.off())
      plot_splice_coverage(model(), answer()$annotation)
    })
    output$coverage_pdf <- shiny::downloadHandler("splice_base_resolution.pdf", function(file) {
      if (capabilities("cairo")) grDevices::cairo_pdf(file, width = 12, height = max(4.5, 2 + 1.2 * nrow(splice_called_groups(model())))) else
        grDevices::pdf(file, width = 12, height = max(4.5, 2 + 1.2 * nrow(splice_called_groups(model()))), useDingbats = FALSE)
      on.exit(grDevices::dev.off())
      plot_splice_coverage(model(), answer()$annotation)
    })
    output$junction_csv <- shiny::downloadHandler("nearby_rna_junction_evidence.csv", function(file)
      utils::write.csv(splice_near_junctions(model()), file, row.names = FALSE, na = "NA"))
    output$evidence_json <- shiny::downloadHandler("splice_evidence.json", function(file) {
      x <- model()
      a <- answer()$annotation
      x$display <- list(depth_scale = "shared across called genotypes", boundary_marker_cap = 6L,
        junction_matrix_cap = 12L, dna_motif = "target-variant substitution only; not phased haplotype or causal inference",
        annotation_metadata = if (is.null(a) || inherits(a, "error")) NULL else a$metadata)
      jsonlite::write_json(x, file, auto_unbox = TRUE, pretty = TRUE, na = "null")
    })
    invisible(list(model = model, answer = answer, model_task = model_task,
      annotation_task = annotation_task, compact_result = compact_result, flank = flank))
  })
}

# GENCODE reference models are an independent, read-only display layer.
# Annotation queries reuse the app's future pool; read/junction metrics are untouched.
annotation_track_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::card(bslib::card_header("GENCODE v48 · reference transcript models"),
    shiny::uiOutput(ns("status")), shiny::uiOutput(ns("track_ui")),
    shiny::uiOutput(ns("display_note")),
    shiny::p("Exons are gray; annotated CDS is thick blue; annotated UTR is thin teal. Arrows show transcript strand. Models are reference annotations, not measured isoform expression or PSI."))
}

annotation_detail_ui <- function(id, allow_confirmation = FALSE) {
  ns <- shiny::NS(id)
  shiny::tagList(
    if (allow_confirmation) shiny::tagList(
      shiny::checkboxInput(ns("confirm_grch38"), "Confirm this BAM uses GRCh38 for GENCODE v48", FALSE),
      shiny::helpText("Only needed for a BAM outside the verified RNA manifest. The selected chromosome's header length must also match GRCh38. This controls annotation only; no liftOver or read reanalysis is performed.")),
    shiny::uiOutput(ns("reference_status")),
    shiny::selectizeInput(ns("priority"), "Transcripts to draw first (optional)", choices = character(), multiple = TRUE),
    shiny::helpText("The track draws at most 20 transcripts: selected IDs first, then MANE_Select, Ensembl_canonical, basic, and finally gene symbol / versioned ID. These GENCODE tags set display priority only; they do not measure isoform expression. Tables retain every overlapping transcript model and its full feature coordinates."),
    bslib::navset_card_tab(
      bslib::nav_panel("Genes", DT::DTOutput(ns("genes"))),
      bslib::nav_panel("Transcripts", DT::DTOutput(ns("transcripts"))),
      bslib::nav_panel("Exons / features", DT::DTOutput(ns("exons")), DT::DTOutput(ns("features"))),
      bslib::nav_panel("Junction annotation",
        DT::DTOutput(ns("junctions")), DT::DTOutput(ns("matches")),
        shiny::p("Junction annotation compares exact 1-based intron boundaries and strand. UNKNOWN_STRAND is not a confirmed match. NOT_IN_ANNOTATION does not establish a novel splice event. Native RNA counts are unchanged.")),
      bslib::nav_panel("Source / export",
        shiny::downloadButton(ns("download_transcripts"), "Transcripts CSV"),
        shiny::downloadButton(ns("download_exons"), "Exons CSV"),
        shiny::downloadButton(ns("download_matches"), "Junction matches CSV"),
        shiny::downloadButton(ns("download_annotation"), "Full annotation JSON"),
        shiny::verbatimTextOutput(ns("metadata"))))
  )
}

annotation_pick_transcripts <- function(transcripts, selected = character(), cap = 20L, window = NULL) {
  if (!is.null(window)) transcripts <- transcripts[transcripts$end1 >= window[[1L]] & transcripts$start1 <= window[[2L]], , drop = FALSE]
  if (!nrow(transcripts)) return(transcripts)
  selected <- unique(as.character(selected))
  rank <- match(transcripts$transcript_id, selected)
  rank[is.na(rank)] <- length(selected) + 1L
  tags <- if ("tags" %in% names(transcripts)) as.character(transcripts$tags) else rep("", nrow(transcripts))
  tags[is.na(tags)] <- ""
  has_tag <- function(tag) grepl(paste0("(^|;)", tag, "(;|$)"), tags)
  tag_rank <- ifelse(has_tag("MANE_Select"), 0L,
    ifelse(has_tag("Ensembl_canonical"), 1L, ifelse(has_tag("basic"), 2L, 3L)))
  ordered <- order(rank, tag_rank, transcripts$gene_name, transcripts$transcript_id, method = "radix", na.last = TRUE)
  transcripts[head(ordered, cap), , drop = FALSE]
}

annotation_plot_models <- function(annotation, window, selected = character(), marker = NULL, cap = 20L) {
  tx <- annotation_pick_transcripts(annotation$transcripts, selected, cap, window)
  old <- graphics::par(mar = c(4.6, 4.1, 2.2, 2.1))
  on.exit(graphics::par(old))
  graphics::plot(NA, xlim = window, ylim = c(.35, max(1L, nrow(tx)) + .85),
    xlab = "Genomic position (1-based)", ylab = "", yaxt = "n", main = "GENCODE v48 · transcript structure")
  if (!nrow(tx)) {
    graphics::text(mean(window), 1, "No annotated transcript spans overlap this display window")
    return(invisible(tx))
  }
  if (length(marker) == 1L && is.finite(marker) && marker >= window[[1L]] && marker <= window[[2L]])
    graphics::abline(v = marker, lty = 2, col = "#A64242", lwd = 1.3)
  exons <- annotation$exons; features <- annotation$features
  span <- diff(window)
  for (i in seq_len(nrow(tx))) {
    row <- nrow(tx) - i + 1
    start <- max(window[[1L]], tx$start1[[i]])
    end <- min(window[[2L]], tx$end1[[i]])
    label <- paste0(substr(tx$gene_name[[i]], 1L, 38L), " · ", tx$transcript_id[[i]], " (", tx$strand[[i]], ")")
    graphics::text(window[[1L]], row + .32, label, adj = c(0, .5), cex = .71)
    if (start > end) next
    graphics::segments(start, row, end, row, col = "#737373", lwd = 1)
    if (end > start && tx$strand[[i]] %in% c("+", "-")) {
      centers <- seq(start, end, length.out = max(3L, min(12L, ceiling((end - start) / max(1, span) * 12L))))
      centers <- centers[centers > start & centers < end]
      width <- min((end - start) / 8, span / 100)
      direction <- if (tx$strand[[i]] == "+") 1 else -1
      graphics::arrows(centers - direction * width / 2, row,
        centers + direction * width / 2, row, length = .055, col = "#555555", angle = 25)
    }
    e <- exons[exons$transcript_id == tx$transcript_id[[i]] & exons$end1 >= window[[1L]] & exons$start1 <= window[[2L]], , drop = FALSE]
    if (nrow(e)) graphics::rect(pmax(e$start1, window[[1L]]), row - .10,
      pmin(e$end1, window[[2L]]), row + .10, col = "#ADB5BD", border = "#697079")
    if (!is.null(features) && nrow(features)) {
      f <- features[features$transcript_id == tx$transcript_id[[i]] & features$end1 >= window[[1L]] & features$start1 <= window[[2L]], , drop = FALSE]
      for (kind in c("UTR", "CDS")) {
        z <- f[f$feature == kind, , drop = FALSE]
        if (!nrow(z)) next
        halfheight <- if (kind == "CDS") .16 else .065
        graphics::rect(pmax(z$start1, window[[1L]]), row - halfheight,
          pmin(z$end1, window[[2L]]), row + halfheight,
          col = if (kind == "CDS") "#245A91" else "#278C8C", border = NA)
      }
    }
  }
  invisible(tx)
}

annotation_server <- function(id, resources, annotation_file, region, display_range,
                              junctions, marker = shiny::reactive(NULL),
                              allow_confirmation = FALSE) {
  shiny::moduleServer(id, function(input, output, session) {
    resource_ok <- !is.null(resources) && !inherits(resources, "error")
    resource_message <- if (inherits(resources, "error")) conditionMessage(resources) else
      if (is.null(resources)) "GENCODE annotation is not configured." else NULL
    unavailable <- shiny::reactiveVal(NULL)
    submitted <- shiny::reactiveVal(NULL)
    request_serial <- shiny::reactiveVal(0L)
    confirmed_source <- shiny::reactiveVal(NULL)
    if (allow_confirmation) {
      shiny::observeEvent(region()$source_id, {
        confirmed_source(NULL)
        shiny::updateCheckboxInput(session, "confirm_grch38", value = FALSE)
      }, ignoreNULL = TRUE)
      shiny::observeEvent(input$confirm_grch38, {
        confirmed_source(if (isTRUE(input$confirm_grch38)) region()$source_id else NULL)
      })
    }
    annotation_task <- shiny::ExtendedTask$new(function(res, interval, native_junctions, source_file, serial) {
      promises::future_promise({
        e <- new.env(parent = globalenv())
        sys.source(source_file, envir = e)
        a <- e$query_annotation(res, interval$chrom, interval$start1, interval$end1, build = interval$build)
        matches <- e$annotate_junctions(native_junctions, a)
        list(annotation = a, junctions = matches$junctions, matches = matches$matches,
          junction_annotation_metadata = matches$metadata, reference_context = interval, request_serial = serial)
      }, seed = TRUE)
    })
    effective_region <- shiny::reactive({
      x <- region()
      if (isTRUE(x$demo) || identical(x$chrom, "chrDemo")) {
        x$unavailable <- "GENCODE is unavailable for the synthetic chrDemo reference."
        return(x)
      }
      unknown_build <- is.null(x$build) || length(x$build) != 1L || is.na(x$build) || x$build %in% c("unknown", "UNKNOWN", "")
      if (!unknown_build && !identical(x$build, "GRCh38")) {
        x$unavailable <- paste("The configured result assembly", x$build, "differs from GENCODE v48 GRCh38; no coordinate conversion is performed.")
        return(x)
      }
      if (!identical(x$build, "GRCh38") && allow_confirmation && isTRUE(input$confirm_grch38) &&
          identical(confirmed_source(), x$source_id)) {
        if (!isTRUE(x$reference_compatible)) {
          x$unavailable <- if (!is.null(x$reference_note)) x$reference_note else
            "The selected BAM chromosome length has not been verified against GRCh38."
          return(x)
        }
        x$build <- "GRCh38"
        x$build_basis <- "User-confirmed GRCh38; selected chromosome header length matches the reference."
      }
      if (!identical(x$build, "GRCh38"))
        x$unavailable <- if (allow_confirmation) paste("This BAM is outside the verified GRCh38 manifest. Confirm its reference in the GENCODE v48 tab.", x$reference_note) else
          "GENCODE v48 requires a GRCh38 result; no coordinate conversion is performed."
      x
    })
    shiny::observeEvent(list(effective_region(), junctions()), {
      if (!resource_ok) return()
      x <- effective_region()
      unavailable(x$unavailable)
      if (!is.null(x$unavailable)) return()
      submitted(x)
      request_serial(request_serial() + 1L)
      annotation_task$invoke(resources, x, junctions(), annotation_file, request_serial())
    }, ignoreNULL = FALSE)
    answer <- shiny::reactive({
      shiny::req(resource_ok, is.null(effective_region()$unavailable))
      value <- annotation_task$result()
      shiny::req(identical(value$request_serial, request_serial()))
      value
    })
    annotation <- shiny::reactive(answer()$annotation)
    shiny::observeEvent(annotation_task$status(), {
      if (annotation_task$status() != "success") return()
      a <- annotation()
      t <- a$transcripts
      labels <- paste(t$gene_name, t$transcript_id, t$strand, sep = " · ")
      keep <- shiny::isolate(input$priority)
      keep <- keep[keep %in% t$transcript_id]
      shiny::updateSelectizeInput(session, "priority", choices = stats::setNames(t$transcript_id, labels), selected = keep, server = TRUE)
    }, ignoreInit = TRUE)
    status_content <- shiny::reactive({
      if (!resource_ok) return(shiny::div(class = "alert alert-secondary", "GENCODE unavailable: ", resource_message, " RNA metrics remain available."))
      x <- tryCatch(effective_region(), error = function(e) NULL)
      if (is.null(x)) return(shiny::p("Load an RNA result to view its GENCODE annotation."))
      if (!is.null(x$unavailable)) return(shiny::div(class = "alert alert-secondary", x$unavailable, " RNA metrics are unchanged."))
      state <- annotation_task$status()
      if (state == "initial") return(shiny::p("GENCODE annotation is ready for the loaded interval."))
      if (state == "running") return(shiny::p("Reading the indexed GENCODE interval..."))
      if (state == "error") {
        message <- tryCatch(annotation_task$result(), error = function(e) conditionMessage(e))
        return(shiny::div(class = "alert alert-warning", "GENCODE query unavailable: ", as.character(message), ". RNA metrics are unchanged."))
      }
      a <- annotation()
      shiny::p("GENCODE v48 · GRCh38 · ", x$chrom, ":", x$start1, "–", x$end1,
        " · ", nrow(a$genes), " genes · ", nrow(a$transcripts), " transcript models. ",
        if (!is.null(x$build_basis)) x$build_basis)
    })
    output$status <- shiny::renderUI(status_content())
    output$reference_status <- shiny::renderUI(status_content())
    selected_transcripts <- shiny::reactive(annotation_pick_transcripts(annotation()$transcripts, input$priority, 20L, display_range()))
    output$track_ui <- shiny::renderUI({
      n <- nrow(selected_transcripts())
      shiny::plotOutput(session$ns("track"), height = max(180L, 46L * n + 85L))
    })
    output$track <- shiny::renderPlot({
      annotation_plot_models(annotation(), display_range(), input$priority, marker(), 20L)
    })
    output$display_note <- shiny::renderUI({
      a <- annotation(); n <- nrow(selected_transcripts()); w <- display_range()
      visible <- sum(a$transcripts$end1 >= w[[1L]] & a$transcripts$start1 <= w[[2L]])
      marker_note <- if (length(marker()) == 1L && is.finite(marker())) {
        if (marker() >= w[[1L]] && marker() <= w[[2L]]) paste0(" Dashed red line: selected variant at ", marker(), ".") else
          paste0(" Selected variant at ", marker(), " lies outside this display window.")
      } else NULL
      shiny::p("Drawing ", n, " of ", visible, " transcript models overlapping this display window; ", nrow(a$transcripts), " models remain in the full tables. Default display priority: MANE_Select, Ensembl_canonical, basic, then gene / ID. Selection changes only the display. Full features may extend beyond the window.",
        marker_note)
    })
    for (name in c("genes", "transcripts", "exons", "features")) local({
      field <- name
      output[[field]] <- DT::renderDT(annotation()[[field]], rownames = FALSE, selection = "none",
        options = list(pageLength = 10, scrollX = TRUE))
    })
    output$junctions <- DT::renderDT(answer()$junctions, rownames = FALSE, selection = "none",
      options = list(pageLength = 10, scrollX = TRUE))
    output$matches <- DT::renderDT(answer()$matches, rownames = FALSE, selection = "none",
      options = list(pageLength = 10, scrollX = TRUE))
    output$metadata <- shiny::renderPrint({
      print(annotation()$metadata)
      cat("\nCoordinates: 1-based inclusive. Transcript models are not expression measurements.\n")
    })
    output$download_transcripts <- shiny::downloadHandler("gencode_v48_transcripts.csv", function(file)
      utils::write.csv(annotation()$transcripts, file, row.names = FALSE, na = "NA"))
    output$download_exons <- shiny::downloadHandler("gencode_v48_exons.csv", function(file)
      utils::write.csv(annotation()$exons, file, row.names = FALSE, na = "NA"))
    output$download_matches <- shiny::downloadHandler("gencode_v48_junction_matches.csv", function(file)
      utils::write.csv(answer()$matches, file, row.names = FALSE, na = "NA"))
    output$download_annotation <- shiny::downloadHandler("gencode_v48_interval_annotation.json", function(file) {
      a <- answer()
      a$display <- list(window = display_range(), marker = marker(), cap = 20L,
        priority = input$priority, plotted_transcript_ids = selected_transcripts()$transcript_id,
        rule = "selected IDs first, then MANE_Select, Ensembl_canonical, basic, gene symbol and versioned transcript ID; display-only cap")
      jsonlite::write_json(a, file, auto_unbox = TRUE, pretty = TRUE, na = "null")
    })
    invisible(list(answer = answer, annotation = annotation, task = annotation_task,
      effective_region = effective_region, selected_transcripts = selected_transcripts))
  })
}

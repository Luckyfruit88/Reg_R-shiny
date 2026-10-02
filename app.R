# Start from this directory: Rscript -e 'shiny::runApp(".", host="127.0.0.1", port=3838)'
# Load deployment settings before package initialization and background workers.
# Shiny sources app.R with the application directory as the working directory.
# An explicit R_ENVIRON_USER remains supported for existing/custom launchers.
local({
  override <- Sys.getenv("R_ENVIRON_USER", unset = "")
  config <- if (nzchar(override)) path.expand(override) else
    file.path(getwd(), "config", "runtime.Renviron")
  if (!file.exists(config)) {
    if (nzchar(override)) stop("R_ENVIRON_USER configuration does not exist: ", config)
    # Preserve unconfigured/demo use; WGS resource validation reports missing inputs.
  } else {
    config <- normalizePath(config, mustWork = TRUE)
    if (!isTRUE(file_test("-f", config)) || file.access(config, 4) != 0L)
      stop("Reg_Shiny runtime configuration must be a readable regular file: ", config)
    loaded <- withCallingHandlers(readRenviron(config), warning = function(w)
      stop("Invalid Reg_Shiny runtime configuration syntax: ", config, call. = FALSE))
    if (!isTRUE(loaded))
      stop("Cannot read Reg_Shiny runtime configuration: ", config)
  }
})

needed <- c("shiny", "bslib", "DT", "data.table", "processx", "future", "promises", "jsonlite")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "), ". Install the dependencies listed in the README.")
if (utils::packageVersion("shiny") < "1.8.1") stop("shiny >= 1.8.1 is required.")
source("runtime_setup.R", local = TRUE)
regshiny_runtime_setup()
library(shiny)
library(bslib)
backend_file <- normalizePath("backend.R", mustWork = TRUE)
backend <- new.env(parent = globalenv())
sys.source(backend_file, envir = backend)
variant_file <- normalizePath("variant_backend.R", mustWork = TRUE)
sys.source(variant_file, envir = backend)
job_source <- normalizePath("variant_jobs.R", mustWork = TRUE)
job_api <- new.env(parent = backend)
sys.source(job_source, envir = job_api)
annotation_file <- normalizePath("annotation_backend.R", mustWork = TRUE)
annotation_backend <- new.env(parent = globalenv())
sys.source(annotation_file, envir = annotation_backend)
source("annotation_module.R", local = TRUE)
splice_file <- normalizePath("splice_evidence_backend.R", mustWork = TRUE)
splice_backend <- new.env(parent = globalenv())
sys.source(splice_file, envir = splice_backend)
source("splice_evidence_module.R", local = TRUE)
source("variant_module.R", local = TRUE)
source("data_sources_module.R", local = TRUE)
app_dir <- normalizePath(getwd(), mustWork = TRUE)
source_file <- normalizePath("data_sources.R", mustWork = TRUE)
source_api <- new.env(parent = globalenv())
sys.source(source_file, envir = source_api)
profile_config <- Sys.getenv("REGSHINY_PROFILE_CONFIG", "")
startup_spec <- tryCatch(source_api$startup_profile_spec(app_dir,
  config_file = if (nzchar(profile_config)) profile_config else NULL), error = function(e) e)
ui_workers <- regshiny_ui_worker_count()
# One dedicated compute worker for NSLOTS 1/2; two for NSLOTS >= 3.
# I(1) keeps the one-worker case asynchronous instead of switching to sequential.
future::plan(future::multisession, workers = I(ui_workers))
onStop(function() future::plan(future::sequential))

single_bam_ui <- function(id, bam_choices = character()) {
  ns <- NS(id)
  initial_source <- if (length(bam_choices)) "server" else "demo"
  layout_sidebar(
  fillable = FALSE,
  sidebar = sidebar(width = 345,
    h5("1 · Data and analysis interval"),
    selectInput(ns("source"), "Data source", c("Server BAM files" = "server", "Synthetic demo (chrDemo)" = "demo"), selected = initial_source),
    conditionalPanel("input.source == 'server'", ns = ns,
      selectInput(ns("bam"), "BAM file (read only)", choices = bam_choices),
      helpText("A matching BAI or CSI index must be next to the BAM. Choose another SCC file collection in Data sources.")),
    selectInput(ns("chrom"), "Chromosome", choices = character()),
    uiOutput(ns("source_help")),
    numericInput(ns("roi_start"), "Analysis start (1-based)", 101, min = 1),
    numericInput(ns("roi_end"), "Analysis end (inclusive)", 450, min = 1),
    tags$details(tags$summary("Filters and RegTools settings"),
      numericInput(ns("mapq"), "Minimum MAPQ", 20, min = 0, max = 255),
      numericInput(ns("baseq"), "Minimum base quality (depth only)", 0, min = 0, max = 93),
      numericInput(ns("anchor"), "RegTools minimum anchor (bp)", 8, min = 1),
      numericInput(ns("min_intron"), "Minimum intron length (bp)", 70, min = 1),
      numericInput(ns("max_intron"), "Maximum intron length (bp)", 500000, min = 1),
      selectInput(ns("strand_mode"), "Library strandedness", c("XS tags" = "XS", "RF / first-strand" = "RF", "FR / second-strand" = "FR")),
      checkboxInput(ns("exclude_dup"), "Exclude marked duplicates", FALSE),
      checkboxInput(ns("nh1"), "Require NH == 1 (exclude missing NH)", FALSE),
      helpText("Always exclude unmapped, secondary, supplementary and QC-failed alignments. A MAPQ threshold does not universally identify uniquely mapped reads.")),
    input_task_button(ns("run"), "Run analysis / apply filters"),
    hr(),
    h5("2 · Explore the loaded result"),
    numericInput(ns("target_start"), "Target intron start (1-based)", 201, min = 1),
    numericInput(ns("target_end"), "Target intron end (inclusive)", 300, min = 1),
    sliderInput(ns("delta"), "Near tolerance at both boundaries (±bp)", min = 0, max = 50, value = 5, step = 1),
    helpText("Near matches must satisfy both boundaries and exclude the exact junction. Select a row in the Junctions table to set the target."),
    uiOutput(ns("display_control"))
  ),
  uiOutput(ns("status")),
  uiOutput(ns("metrics")),
  p(textOutput(ns("denominator_note"))),
  navset_card_tab(
    nav_panel("Plots",
      card(card_header("Read depth · all positions, including zero coverage"), plotOutput(ns("depth_plot"), height = 250)),
      card(card_header("Junction arcs · top 30 junctions in the display window"), plotOutput(ns("junction_plot"), height = 260)),
      annotation_track_ui(ns("single_annotation")),
      p("Zoom changes the display only. The analysis interval and denominator stay fixed. Arcs are not transcript annotations.")),
    nav_panel("Junctions", DT::DTOutput(ns("junction_table"))),
    nav_panel("GENCODE v48", annotation_detail_ui(ns("single_annotation"), allow_confirmation = TRUE)),
    nav_panel("Read evidence",
      checkboxInput(ns("evidence_all"), "Show all accepted junction evidence (otherwise show exact / near matches)", FALSE),
      DT::DTOutput(ns("evidence_table")),
      p("Each row represents one CIGAR N operation in an alignment. A read_id may appear more than once. Exports omit SEQ and QUAL.")),
    nav_panel("Export and audit",
      downloadButton(ns("download_summary"), "Metrics CSV"),
      downloadButton(ns("download_junctions"), "Junction CSV"),
      downloadButton(ns("download_depth"), "Depth TSV"),
      downloadButton(ns("download_evidence"), "All read evidence CSV"),
      downloadButton(ns("download_metadata"), "Parameters / logs JSON"),
      hr(), verbatimTextOutput(ns("audit")),
      h5("Metric definitions"),
      p("Reads are retained primary alignment records; paired ends count separately. These are not fragment or UMI counts. Summary metrics combine all strands, including unknown strand; the junction table retains strand information."),
      p("The denominator counts reads with M, = or X blocks overlapping the fixed analysis interval. An interval spanned only by N or D does not count as covered."),
      p("Junction reads are the union of read IDs supporting accepted RegTools junctions. Exact and near counts each use a union of matching read IDs; near excludes exact junction events."),
      p("Junctions per 100 reads = 100 × exact count / denominator reads. A zero denominator gives NA. This is not PSI."),
      p("Mean depth averages all positions in the fixed interval, including introns and zero coverage. The export also includes mean depth in up to 50 bp on each side of the target.")),
    nav_panel("Scope and limitations",
      p("This local-analysis prototype is for research and is not clinically validated. BAM files are read on the server and are not modified."),
      p("Use only on a trusted server or behind an authenticated gateway. The app does not implement authentication, cross-user permissions, cancellation or a shared job queue. Do not expose it directly to the public internet."),
      p("Limits: 250 kb and 100,000 candidate alignments per analysis. Exceeding a limit produces an error, not subsampling. Target and tolerance changes reuse loaded results; filter changes require a new analysis."),
      p("This tab analyzes one BAM. Use WGS genotype comparison for matched multi-sample exploration. CRAM, fragment/UMI counts, strand-specific denominators and formal PSI are not supported. GENCODE v48 reference models are available for compatible GRCh38 inputs."))
  )
)
}

single_bam_server <- function(id, bundle, is_active = reactive(TRUE)) {
  force(bundle)
  moduleServer(id, function(input, output, session) {
  bams <- unname(bundle$bam_choices)
  variant_resources <- bundle$variant_resources
  annotation_resources <- bundle$annotation_resources
  contigs <- reactiveVal(NULL)
  contig_error <- reactiveVal(NULL)
  previous_source <- reactiveVal(NULL)
  observeEvent(list(input$source, input$bam), {
    req(is_active())
    demo <- identical(input$source, "demo")
    changed_source <- !identical(previous_source(), input$source)
    previous_source(input$source)
    info <- tryCatch({
      if (demo) data.frame(chrom = "chrDemo", length_bp = 1000)
      else {
        if (is.null(input$bam) || !input$bam %in% unname(bams))
          stop("Select an available BAM file.")
        backend$read_bam_contigs(input$bam)
      }
    }, error = function(e) e)
    if (inherits(info, "error")) {
      contigs(NULL)
      contig_error(conditionMessage(info))
      updateSelectInput(session, "chrom", choices = character())
      return()
    }
    contig_error(NULL)
    contigs(info)
    chosen <- backend$choose_bam_contig(info, if (changed_source) NULL else isolate(input$chrom))
    freezeReactiveValue(input, "chrom")
    updateSelectInput(session, "chrom", choices = info$chrom, selected = chosen)
    if (changed_source) {
      coords <- if (demo) c(101, 450, 201, 300) else
        if (chosen %in% c("chr1", "1") && info$length_bp[match(chosen, info$chrom)] >= 17000)
          c(14000, 17000, 14830, 14969) else rep(NA_real_, 4)
      for (i in seq_along(coords))
        updateNumericInput(session, c("roi_start", "roi_end", "target_start", "target_end")[[i]],
                           value = coords[[i]])
    }
  }, ignoreNULL = FALSE)
  output$source_help <- renderUI({
    if (!is.null(contig_error()))
      return(div(class = "alert alert-danger", contig_error()))
    if (identical(input$source, "demo"))
      return(helpText("chrDemo is an artificial reference used only by the synthetic demo. Select Server BAM files for real chromosomes."))
    info <- contigs()
    if (is.null(info)) return(helpText("Reading chromosome names from the selected BAM..."))
    helpText("Chromosome choices come from this BAM's header. The initial chr1:14000-17000 interval is an example, not automatic gene or junction detection. Enter coordinates appropriate for your BAM and reference assembly.")
  })

  current_config <- reactive({
    list(demo = identical(input$source, "demo"),
      bam = if (identical(input$source, "demo")) "" else input$bam,
      chrom = trimws(input$chrom), start1 = input$roi_start, end1 = input$roi_end,
      mapq = input$mapq, baseq = input$baseq, anchor = input$anchor,
      min_intron = input$min_intron, max_intron = input$max_intron,
      strand_mode = input$strand_mode, exclude_duplicates = input$exclude_dup, nh1_only = input$nh1)
  })
  launched_config <- reactiveVal(NULL)
  launched_contigs <- reactiveVal(NULL)
  task <- ExtendedTask$new(function(cfg, source_file) {
    promises::future_promise({
      e <- new.env(parent = globalenv())
      sys.source(source_file, envir = e)
      e$analyze_bam(cfg)
    }, seed = TRUE)
  }) |> bind_task_button("run")
  observeEvent(input$run, {
    req(is_active())
    cfg <- current_config()
    if (!cfg$demo && (is.null(cfg$bam) || !cfg$bam %in% unname(bams))) {
      showNotification("No available BAM is selected. Configure REGTOOLS_BAM_DIR and restart the app.", type = "error")
      bslib::update_task_button("run", state = "ready", session = session)
      return()
    }
    info <- contigs()
    problem <- tryCatch({
      if (!is.null(contig_error())) stop(contig_error())
      if (is.null(info) || length(cfg$chrom) != 1L || !cfg$chrom %in% info$chrom)
        stop("Select a chromosome from the selected BAM.")
      backend$validate_config(cfg)
      if (cfg$end1 > info$length_bp[match(cfg$chrom, info$chrom)])
        stop("The analysis end exceeds the selected chromosome length.")
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(problem)) {
      showNotification(problem, type = "error", duration = 10)
      bslib::update_task_button("run", state = "ready", session = session)
      return()
    }
    cfg$data_source_id <- bundle$id
    cfg$data_source_label <- bundle$label
    launched_config(cfg)
    launched_contigs(info)
    task$invoke(cfg, backend_file)
  })
  result <- reactive({req(is_active()); task$result()})
  metrics <- reactive({
    r <- result()
    out <- tryCatch(backend$summarize_target(r, input$target_start, input$target_end, input$delta),
                    error = function(e) e)
    validate(need(!inherits(out, "error"), if (inherits(out, "error")) conditionMessage(out) else ""))
    out
  })
  output$status <- renderUI({
    s <- task$status()
    if (s == "initial") return(p(if (identical(input$source, "demo")) "Synthetic demo ready: 4 reads, exact=2, near=1, and 50 junctions per 100 reads are expected." else "Choose a BAM, chromosome and analysis interval, then click Run analysis. The initial interval is an example."))
    if (s == "running") return(p("Analyzing the submitted parameter snapshot. Unverified counts are hidden until the analysis completes."))
    if (s == "error") {
      err <- tryCatch(task$result(), error = function(e) conditionMessage(e))
      return(div(class = "alert alert-danger", "Analysis failed: ", as.character(err)))
    }
    r <- task$result()
    previous <- launched_config(); previous$data_source_id <- previous$data_source_label <- NULL
    stale <- !identical(current_config(), previous)
    tagList(
      div(class = if (stale) "alert alert-warning" else "alert alert-success",
          if (stale) "Inputs have changed. The results below are from the previous analysis. Run analysis to apply the new inputs." else "Analysis complete. Native RegTools scores match the read-level CIGAR counts."),
      if (length(r$warnings)) div(class = "alert alert-warning", paste(r$warnings, collapse = "; ")))
  })
  output$metrics <- renderUI({
    m <- metrics()
    num <- function(x, digits = 0L) if (is.na(x)) "NA" else format(round(x, digits), big.mark = ",", trim = TRUE)
    layout_column_wrap(width = 1/3,
      value_box("Junction reads (union within interval)", num(m$junction_reads)),
      value_box("Read depth (interval mean)", num(m$mean_read_depth, 3)),
      value_box("Junction exact count", num(m$junction_exact_count)),
      value_box("Junction near count (excluding exact)", num(m$junction_near_count)),
      value_box("Junction per 100 reads", num(m$junction_per_100_reads, 2)))
  })
  output$denominator_note <- renderText({
    m <- metrics()
    paste0("Fixed denominator: ", m$denominator_reads, " primary alignments; interval ",
      m$chrom, ":", m$denominator_start1, "-", m$denominator_end1,
      ". The per-100 metric uses the exact count, not the near count or mean depth.")
  })
  output$display_control <- renderUI({
    r <- result()
    sliderInput("display", "Display window (analysis interval stays fixed)", min = r$config$start1,
                max = r$config$end1, value = c(r$config$start1, r$config$end1), step = 1, sep = "")
  })
  view_range <- reactive({
    r <- result(); x <- input$display
    if (is.null(x) || length(x) != 2L || x[[1L]] < r$config$start1 || x[[2L]] > r$config$end1)
      c(r$config$start1, r$config$end1) else x
  })
  annotation_server("single_annotation", resources = annotation_resources, annotation_file = annotation_file,
    region = reactive({
      r <- result(); cfg <- r$config
      x <- list(chrom = cfg$chrom, start1 = cfg$start1, end1 = cfg$end1,
        demo = isTRUE(cfg$demo), source_id = cfg$bam, build = "unknown",
        reference_compatible = FALSE)
      if (x$demo) return(x)
      verified <- !inherits(variant_resources, "error") && cfg$bam %in% variant_resources$samples$bam
      if (verified) {
        x$build <- variant_resources$build
        x$build_basis <- "BAM is linked in the configured DNA/RNA reference manifest."
      }
      expected_length <- if (!inherits(annotation_resources, "error"))
        tryCatch(annotation_backend$annotation_reference_length(annotation_resources, cfg$chrom), error = function(e) NA_real_) else NA_real_
      if ((!length(expected_length) || !is.finite(expected_length)) && !inherits(variant_resources, "error") &&
          identical(variant_resources$build, "GRCh38") &&
          "length_bp" %in% names(variant_resources$vcf_registry)) {
        expected_length <- suppressWarnings(as.numeric(variant_resources$vcf_registry$length_bp[match(cfg$chrom, variant_resources$vcf_registry$chrom)]))
      }
      info <- launched_contigs()
      observed_length <- info$length_bp[match(cfg$chrom, info$chrom)]
      if (length(expected_length) == 1L && is.finite(expected_length) && length(observed_length) == 1L && is.finite(observed_length)) {
        x$reference_compatible <- identical(as.numeric(observed_length), as.numeric(expected_length))
        x$reference_note <- if (x$reference_compatible) paste("The selected chromosome header length matches GRCh38:", expected_length, "bp.") else
          paste("GENCODE unavailable: BAM chromosome length", observed_length, "differs from the GRCh38 length", expected_length, ".")
        if (!x$reference_compatible) x$unavailable <- x$reference_note
      } else x$reference_note <- "A reference chromosome length is unavailable; annotation requires compatible reference metadata."
      x
    }), display_range = view_range, junctions = reactive(result()$junctions), allow_confirmation = TRUE)
  output$depth_plot <- renderPlot({
    r <- result(); w <- view_range()
    d <- r$depth[r$depth$pos1 >= w[[1L]] & r$depth$pos1 <= w[[2L]], , drop = FALSE]
    plot(d$pos1, d$depth, type = "l", xlab = paste0(r$config$chrom, " · 1-based position"),
         ylab = "Read depth", main = "Aligned-base coverage", ylim = c(0, max(1, d$depth)))
    abline(v = c(input$target_start, input$target_end), lty = 2)
  })
  junction_rows <- reactive({
    r <- result(); j <- r$junctions
    j <- j[order(-j$score, j$intron_start1), , drop = FALSE]
    exact <- j$intron_start1 == input$target_start & j$intron_end1 == input$target_end
    near <- abs(j$intron_start1 - input$target_start) <= input$delta &
      abs(j$intron_end1 - input$target_end) <= input$delta & !exact
    j$target_class <- ifelse(exact, "exact", ifelse(near, "near", "other"))
    j
  })
  output$junction_plot <- renderPlot({
    w <- view_range(); j <- junction_rows()
    j <- head(j[j$intron_start1 >= w[[1L]] & j$intron_end1 <= w[[2L]], , drop = FALSE], 30L)
    plot(NA, xlim = w, ylim = c(0, 1.25), xlab = "Intron positions (1-based)", ylab = "",
         yaxt = "n", main = "Native RegTools junction support")
    if (!nrow(j)) { text(mean(w), .5, "No junctions fully inside this display window"); return(invisible(NULL)) }
    for (i in seq_len(nrow(j))) {
      t <- seq(0, 1, length.out = 80)
      height <- .2 + .75 * j$score[[i]] / max(j$score)
      x <- j$intron_start1[[i]] + t * (j$intron_end1[[i]] - j$intron_start1[[i]])
      y <- 4 * height * t * (1 - t)
      lines(x, y, lwd = if (j$target_class[[i]] == "exact") 3 else 1,
            lty = if (j$target_class[[i]] == "near") 2 else 1)
      text(mean(range(x)), height + .035, labels = paste0(j$score[[i]], " ", j$strand[[i]]), cex = .8)
    }
  })
  output$junction_table <- DT::renderDT({
    j <- junction_rows()
    j[, c("chrom", "intron_start1", "intron_end1", "strand", "score", "target_class",
           "intron_start0", "intron_end0", "blockSizes"), drop = FALSE]
  }, rownames = FALSE, selection = "single", options = list(pageLength = 10, scrollX = TRUE))
  observeEvent(input$junction_table_rows_selected, {
    req(is_active())
    i <- input$junction_table_rows_selected
    j <- junction_rows()
    if (length(i) == 1L && i <= nrow(j)) {
      updateNumericInput(session, "target_start", value = j$intron_start1[[i]])
      updateNumericInput(session, "target_end", value = j$intron_end1[[i]])
    }
  })
  all_evidence <- reactive({
    r <- result(); e <- r$events
    extra <- setdiff(names(r$reads), c("read_id", "chrom", "strand"))
    ans <- cbind(e, r$reads[match(e$read_id, r$reads$read_id), extra, drop = FALSE])
    rownames(ans) <- NULL
    ans
  })
  output$evidence_table <- DT::renderDT({
    e <- all_evidence()
    if (!input$evidence_all) e <- e[
      abs(e$intron_start1 - input$target_start) <= input$delta &
      abs(e$intron_end1 - input$target_end) <= input$delta, , drop = FALSE]
    e[, setdiff(names(e), "key"), drop = FALSE]
  }, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  output$audit <- renderPrint({
    r <- result()
    cat("RegTools score / CIGAR evidence:", if (r$native_audit_passed) "PASS" else "FAIL", "\n")
    cat("Candidate alignments:", r$candidate_alignments, "\n")
    cat("Excluded span-only overlaps (no M/= /X coverage):", r$span_only_excluded, "\n")
    cat("Retained alignments:", nrow(r$reads), "\n")
    cat("Completed at:", r$completed_at, "\n\n")
    print(metrics())
  })
  output$download_summary <- downloadHandler("junction_metrics.csv", function(file)
    utils::write.csv(metrics(), file, row.names = FALSE, na = "NA", fileEncoding = "UTF-8"))
  output$download_junctions <- downloadHandler("junctions.csv", function(file)
    utils::write.csv(junction_rows(), file, row.names = FALSE, fileEncoding = "UTF-8"))
  output$download_depth <- downloadHandler("depth.tsv", function(file)
    utils::write.table(result()$depth, file, sep = "\t", quote = FALSE, row.names = FALSE))
  output$download_evidence <- downloadHandler("junction_read_evidence.csv", function(file)
    utils::write.csv(all_evidence(), file, row.names = FALSE, na = "NA", fileEncoding = "UTF-8"))
  output$download_metadata <- downloadHandler("parameters_and_audit.json", function(file) {
    r <- result()
    x <- list(config = r$config, target_metrics = metrics(), versions = r$versions,
      source_metadata = r$source_metadata, completed_at = r$completed_at,
      native_audit_passed = r$native_audit_passed, warnings = r$warnings, commands = r$log,
      definitions = list(count_unit = "primary alignment records, paired ends separate",
        denominator = "retained reads with M/= /X blocks overlapping fixed analysis ROI",
        near = "both boundaries within delta; exclude exact events; union read IDs",
        per100 = "100 * exact / denominator; NA when denominator is zero",
        depth = "samtools depth; zero-inclusive mean; BQ affects depth only; pair overlaps counted twice"))
    jsonlite::write_json(x, file, auto_unbox = TRUE, pretty = TRUE, na = "null")
  })
  invisible(list(result = result, task = task, current_config = current_config))
  })
}

ui <- page_navbar(
  title = "Reg_Shiny", id = "workspace", selected = "Data sources", fillable = FALSE,
  header = uiOutput("active_data_banner"),
  nav_panel("Data sources", data_sources_ui("sources")),
  nav_panel("WGS genotype comparison", uiOutput("variant_workspace")),
  nav_panel("Single BAM / demo", uiOutput("single_workspace"))
)

server <- function(input, output, session) {
  unavailable <- simpleError("Connect a dataset in the Data sources tab.")
  demo_bundle <- list(id = "unconfigured", label = "Unconfigured / synthetic demo", bam_choices = character(),
    variant_resources = unavailable, annotation_resources = unavailable, splice_reference = unavailable, job_root = NULL)
  current_view <- reactiveVal(list(generation = 0L, variant_id = NULL, single_id = "single_0", bundle = demo_bundle))
  single_bam_server("single_0", demo_bundle, is_active = reactive(identical(current_view()$generation, 0L)))
  sources <- data_sources_server("sources", source_api, source_file, app_dir, startup_spec = startup_spec)
  generation <- 0L
  observeEvent(sources$bundle(), {
    b <- sources$bundle(); req(!is.null(b))
    generation <<- generation + 1L
    local({
      frozen_bundle <- b
      this_generation <- generation
      variant_id <- paste0("variant_", this_generation)
      single_id <- paste0("single_", this_generation)
      current_view(list(generation = this_generation, variant_id = variant_id, single_id = single_id, bundle = frozen_bundle))
      enabled <- reactive(identical(current_view()$generation, this_generation))
      variant_server(variant_id, backend, frozen_bundle$variant_resources, backend_file, variant_file,
        job_api = job_api, job_source = job_source, job_root = frozen_bundle$job_root,
        annotation_resources = frozen_bundle$annotation_resources, annotation_file = annotation_file,
        splice_reference = frozen_bundle$splice_reference, splice_file = splice_file,
        source_identity = list(id = frozen_bundle$id, label = frozen_bundle$label), is_active = enabled)
      single_bam_server(single_id, frozen_bundle, is_active = enabled)
    })
    bslib::nav_select("workspace", if (inherits(b$variant_resources, "error")) "Single BAM / demo" else "WGS genotype comparison", session = session)
  }, ignoreNULL = TRUE)
  output$active_data_banner <- renderUI({
    b <- current_view()$bundle
    if (identical(b$id, "unconfigured")) return(helpText("Prepare a dataset in Data sources. The synthetic demo remains available while preparation runs."))
    div(class = "alert alert-secondary", strong(paste("Connected dataset:", b$label)), " · source identity ", code(b$id))
  })
  output$variant_workspace <- renderUI({
    v <- current_view()
    if (is.null(v$variant_id)) return(p("Connect FHS or custom SCC files in Data sources to start a WGS lookup."))
    variant_ui(v$variant_id)
  })
  output$single_workspace <- renderUI({
    v <- current_view(); single_bam_ui(v$single_id, v$bundle$bam_choices)
  })
}
shinyApp(ui, server)

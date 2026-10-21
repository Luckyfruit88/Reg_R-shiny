# Session-scoped connections to files already on SCC. No data upload and no
# process-wide environment mutation is used to change an active dataset.

regshiny_ui_worker_count <- function(slots = Sys.getenv("NSLOTS", "1"), override = Sys.getenv("REGSHINY_UI_WORKERS", "")) {
  n <- suppressWarnings(as.integer(slots))
  if (length(n) != 1L || is.na(n) || n < 1L || length(slots) != 1L || !grepl("^[1-9][0-9]*$", slots)) n <- 1L
  budget <- if (n >= 3L) 2L else 1L
  if (!nzchar(override)) return(budget)
  requested <- suppressWarnings(as.integer(override))
  if (length(requested) != 1L || is.na(requested) || requested < 1L || requested > budget ||
      !grepl("^[12]$", override)) stop("REGSHINY_UI_WORKERS must be 1 or 2 and may not exceed the current SCC core budget.")
  requested
}

data_sources_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_columns(col_widths = c(6, 6),
    bslib::card(bslib::card_header("Connect files already on SCC"),
      shiny::radioButtons(ns("kind"), "Dataset", c("FHS on SCC (default)" = "fhs", "Custom SCC files" = "custom"), selected = "fhs"),
      shiny::conditionalPanel(sprintf("input['%s'] === 'fhs'", ns("kind")),
        shiny::p("FHS is prepared automatically when the app opens. Your account must have permission to read its DNA and RNA data. Each user gets a private prepared profile and job history."),
        shiny::tags$details(shiny::tags$summary("FHS source locations"), shiny::verbatimTextOutput(ns("fhs_paths")))),
      shiny::conditionalPanel(sprintf("input['%s'] === 'custom'", ns("kind")),
        shiny::textInput(ns("label"), "Dataset name", "Custom"),
        shiny::textInput(ns("build"), "Declared reference assembly", "GRCh38"),
        shiny::textInput(ns("vcf"), "Indexed VCF / BCF file (optional for BAM-only use)", placeholder = "/restricted/projectnb/.../cohort.vcf.gz"),
        shiny::textInput(ns("vcf_registry"), "Or a VCF registry TSV", placeholder = "Columns: chrom, vcf, build"),
        shiny::selectInput(ns("mapping_mode"), "DNA-to-RNA sample matching",
          c("Explicit mapping TSV" = "explicit", "Verify exact BAM @RG SM against VCF IDs" = "read_group")),
        shiny::conditionalPanel(sprintf("input['%s'] === 'explicit'", ns("mapping_mode")),
          shiny::textInput(ns("mapping"), "Mapping TSV (vcf_sample, bam; optional build)", placeholder = "BAM paths must be absolute SCC paths")),
        shiny::textInput(ns("bam_dir"), "BAM directory (optional if mapping lists the BAMs)", placeholder = "Direct .bam children; no recursive discovery"),
        shiny::textInput(ns("bam"), "Or one BAM file", placeholder = "/restricted/projectnb/.../sample.bam"),
        shiny::helpText("Genotype comparison requires an exact verified sample mapping. File names are never used to guess DNA sample identity. For Single BAM, VCF and mapping may be left blank."),
        shiny::tags$details(shiny::tags$summary("Reference sequence and GENCODE (optional)"),
          shiny::textInput(ns("reference_fasta"), "Indexed reference FASTA", placeholder = "Matching .fai must exist"),
          shiny::textInput(ns("reference_receipt"), "Reference provenance receipt JSON (optional)"),
          shiny::textInput(ns("annotation_db"), "Existing Reg_Shiny GENCODE SQLite index"),
          shiny::textInput(ns("annotation_gtf"), "Or GENCODE v48 GTF / GTF.gz to index once"),
          shiny::helpText("GENCODE v48 requires GRCh38. Missing optional reference resources are reported without disabling saved RNA counts."))),
      bslib::input_task_button(ns("prepare"), "Prepare and connect"),
      shiny::p("Enter server paths; files are read on SCC. Nothing is uploaded through this page. A failed preparation leaves the current dataset connected.")),
    bslib::card(bslib::card_header("Connection and preparation"),
      shiny::uiOutput(ns("status")), shiny::uiOutput(ns("active")),
      shiny::hr(), shiny::h5("Your prepared datasets"),
      shiny::selectInput(ns("saved"), "Saved profile", choices = character()),
      shiny::actionButton(ns("refresh"), "Refresh list"),
      shiny::actionButton(ns("load"), "Connect selected profile"),
      shiny::helpText("Reconnecting validates the frozen sources again. A successful switch opens a fresh viewer; lookups and plots from another dataset are never relabelled as the new source. Persistent SCC jobs continue under their original dataset.")))
}

data_sources_server <- function(id, source_api, source_file, app_dir, startup_spec = NULL, state_root = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    active_bundle <- shiny::reactiveVal(NULL)
    submitted <- shiny::reactiveVal(NULL)
    revision <- shiny::reactiveVal(0L)
    applied_revision <- shiny::reactiveVal(0L)
    activation_error <- shiny::reactiveVal(NULL)
    progress_path <- shiny::reactiveVal(NULL)
    progress_files <- character()
    saved_profiles <- shiny::reactiveVal(data.frame(profile_dir = character(), id = character(), label = character()))
    initial <- if (is.null(startup_spec)) tryCatch(source_api$default_fhs_spec(), error = function(e) e) else startup_spec
    source_task <- shiny::ExtendedTask$new(function(operation, payload, code_file, application_dir, root, progress, serial) {
      promises::future_promise({
        e <- new.env(parent = globalenv()); sys.source(code_file, envir = e)
        bundle <- if (identical(operation, "prepare"))
          e$prepare_data_profile(payload, application_dir, state_root = root, progress_file = progress, build_annotation = TRUE) else
          e$load_data_profile(payload, application_dir, state_root = root)
        list(bundle = bundle, revision = serial)
      }, seed = TRUE)
    }) |> bslib::bind_task_button("prepare")
    refresh_profiles <- function() {
      found <- tryCatch(source_api$list_data_profiles(app_dir, state_root = state_root), error = function(e) e)
      if (inherits(found, "error")) return(invisible(NULL))
      saved_profiles(found)
      choices <- if (nrow(found)) stats::setNames(found$profile_dir, paste(found$label, found$id, sep = " · ")) else character()
      current <- shiny::isolate(input$saved)
      shiny::updateSelectInput(session, "saved", choices = choices,
        selected = if (length(current) == 1L && current %in% found$profile_dir) current else NULL)
      invisible(found)
    }
    launch <- function(operation, payload, label) {
      activation_error(NULL)
      revision(shiny::isolate(revision()) + 1L)
      path <- tempfile("regshiny-profile-progress-", fileext = ".json")
      progress_files <<- c(progress_files, path)
      progress_path(path)
      submitted(list(operation = operation, label = label, revision = revision()))
      source_task$invoke(operation, payload, source_file, app_dir, state_root, path, revision())
    }
    build_spec <- function() {
      if (identical(input$kind, "fhs")) return(source_api$default_fhs_spec())
      if (!identical(input$kind, "custom")) stop("Choose FHS or Custom SCC files.")
      value <- function(id) { x <- input[[id]]; if (is.null(x)) "" else trimws(x) }
      source_api$custom_source_spec(label = value("label"), build = value("build"),
        vcf = value("vcf"), vcf_registry = value("vcf_registry"), mapping = value("mapping"),
        bam_dir = value("bam_dir"), bam = value("bam"), annotation_db = value("annotation_db"),
        annotation_gtf = value("annotation_gtf"), reference_fasta = value("reference_fasta"),
        reference_receipt = value("reference_receipt"), mapping_mode = value("mapping_mode"))
    }
    shiny::observeEvent(input$prepare, {
      tryCatch({
        spec <- build_spec()
        launch("prepare", spec, spec$label)
      }, error = function(e) {
        shiny::showNotification(conditionMessage(e), type = "error", duration = 12)
        bslib::update_task_button("prepare", state = "ready", session = session)
      })
    })
    shiny::observeEvent(input$refresh, refresh_profiles())
    shiny::observeEvent(input$load, {
      tryCatch({
        path <- input$saved
        rows <- saved_profiles()
        if (length(path) != 1L || !path %in% rows$profile_dir) stop("Select one of your prepared profiles.")
        launch("load", path, rows$label[match(path, rows$profile_dir)])
      }, error = function(e) shiny::showNotification(conditionMessage(e), type = "error", duration = 12))
    })
    shiny::observeEvent(TRUE, {
      refresh_profiles()
      if (inherits(initial, "error")) return()
      shiny::updateRadioButtons(session, "kind", selected = initial$kind)
      if (identical(initial$kind, "custom")) {
        for (field in c("label", "build", "vcf", "vcf_registry", "mapping", "bam_dir", "bam", "annotation_db", "annotation_gtf", "reference_fasta", "reference_receipt"))
          if (length(initial[[field]]) == 1L) shiny::updateTextInput(session, field, value = initial[[field]])
        shiny::updateSelectInput(session, "mapping_mode", selected = initial$mapping_mode)
      }
      launch("prepare", initial, initial$label)
    }, once = TRUE)
    shiny::observeEvent(source_task$status(), {
      if (source_task$status() != "success") return()
      value <- source_task$result()
      if (!identical(value$revision, revision()) || identical(value$revision, applied_revision())) return()
      bundle <- value$bundle
      required <- c("id", "label", "variant_resources", "bam_choices", "annotation_resources", "splice_reference", "job_root")
      if (!is.list(bundle) || !all(required %in% names(bundle)) || length(bundle$id) != 1L ||
          is.na(bundle$id) || !nzchar(bundle$id)) {
        activation_error("The prepared dataset has an invalid resource bundle; the existing connection was retained.")
        shiny::showNotification(activation_error(), type = "error", duration = 12)
        return()
      }
      active_bundle(bundle)
      applied_revision(value$revision)
      refresh_profiles()
    }, ignoreInit = TRUE)
    progress <- shiny::reactivePoll(1000, session,
      checkFunc = function() {
        path <- progress_path()
        if (is.null(path) || !file.exists(path)) return(NULL)
        info <- file.info(path)
        paste(path, info$size, as.numeric(info$mtime))
      }, valueFunc = function() {
        path <- progress_path()
        if (is.null(path) || !file.exists(path)) return(NULL)
        tryCatch(jsonlite::read_json(path, simplifyVector = TRUE), error = function(e) NULL)
      })
    output$fhs_paths <- shiny::renderPrint({
      spec <- source_api$default_fhs_spec()
      print(spec[intersect(c("rna_map", "bam_map", "vcf_dir", "build", "annotation_gtf", "annotation_db", "reference_fasta"), names(spec))])
    })
    output$status <- shiny::renderUI({
      if (!is.null(activation_error())) return(shiny::div(class = "alert alert-warning", activation_error()))
      if (inherits(initial, "error") && source_task$status() == "initial")
        return(shiny::div(class = "alert alert-warning", "Startup dataset is unavailable: ", conditionMessage(initial), ". Choose Custom SCC files to connect other data."))
      state <- source_task$status()
      request <- submitted()
      if (state == "initial") return(shiny::p("Ready to prepare a dataset."))
      if (state == "running") {
        p <- progress()
        return(shiny::div(class = "alert alert-info", "Preparing ", request$label, "…",
          if (!is.null(p)) shiny::tagList(shiny::tags$br(), p$stage, " · ", p$message),
          shiny::tags$br(), "The current connection changes only after validation succeeds."))
      }
      if (state == "error") {
        msg <- tryCatch(source_task$result(), error = function(e) conditionMessage(e))
        return(shiny::div(class = "alert alert-warning", "Preparation failed: ", as.character(msg),
          if (is.null(active_bundle())) " No dataset is connected; Single BAM demo remains available." else " The previous dataset remains connected."))
      }
      shiny::div(class = "alert alert-success", "Preparation and source validation completed. The connected dataset is shown below.")
    })
    output$active <- shiny::renderUI({
      b <- active_bundle()
      if (is.null(b)) return(shiny::p("No data connection yet."))
      vr <- b$variant_resources
      shiny::tagList(shiny::h5(paste("Connected:", b$label)), shiny::p("Source identity: ", shiny::code(b$id)),
        shiny::p("Available BAM files: ", length(b$bam_choices)),
        if (inherits(vr, "error")) shiny::p("WGS genotype lookup unavailable: ", conditionMessage(vr)) else
          shiny::p("Reference: ", vr$build, " · registered chromosomes: ", nrow(vr$vcf_registry), " · exact DNA/RNA mappings: ", nrow(vr$samples)),
        if (length(b$issues)) shiny::div(class = "alert alert-warning", paste(unlist(b$issues, use.names = FALSE), collapse = "; ")),
        shiny::helpText("Changing the form does not modify this connection. Use Prepare and connect to apply another dataset."))
    })
    session$onSessionEnded(function() unlink(progress_files[file.exists(progress_files)]))
    invisible(list(bundle = active_bundle, task = source_task, submitted = submitted,
      refresh_profiles = refresh_profiles, saved_profiles = saved_profiles))
  })
}

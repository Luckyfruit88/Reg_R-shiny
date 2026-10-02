# Variant-centered interface inside Reg_Shiny. Preview workers share the app's
# future pool. Full comparisons use persistent SCC jobs and survive disconnects.
variant_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    fillable = FALSE,
    sidebar = bslib::sidebar(width = 345,
      h5("1 · Find a WGS variant"),
      uiOutput(ns("resource_status")),
      textInput(ns("query"), "Variant coordinate (1-based)", placeholder = "Example: chr21:14288395:G:A"),
      helpText("Use the reference assembly shown above. Lookup matches the VCF record's exact POS; REF and ALT identify the record."),
      bslib::input_task_button(ns("lookup"), "Find variant"),
      selectInput(ns("record"), "Select the REF / ALT record", choices = character()),
      hr(), h5("2 · Compare matched RNA samples"),
      numericInput(ns("flank"), "Bases on each side of the variant", 1000, min = 50, max = 124999, step = 50),
      helpText("Widen the interval to inspect longer introns. Large / high-depth regions can reach the 100,000 candidate-alignment limit per sample; these failures are reported without subsampling."),
      radioButtons(ns("sampling_mode"), "RNA sample scope",
        c("All matched BAMs" = "all", "Optional exploratory preview" = "preview"), selected = "all"),
      conditionalPanel(sprintf("input['%s'] === 'preview'", ns("sampling_mode")),
        numericInput(ns("sample_cap"), "Preview: maximum RNA samples per genotype", 3, min = 1, max = 10),
        helpText("Preview only: first sample IDs in sorted order within each called genotype, up to this cap and 30 total. Selection does not depend on RNA outcomes.")),
      conditionalPanel(sprintf("input['%s'] === 'all'", ns("sampling_mode")),
        helpText("All matched BAMs with complete DNA genotype calls are processed. The SCC job and sample checkpoints persist after the browser is closed; return to a saved job below to see progress or results.")),
      tags$details(tags$summary("Read and junction filters"),
        numericInput(ns("mapq"), "Minimum MAPQ", 20, min = 0, max = 255),
        numericInput(ns("baseq"), "Minimum base quality (depth / RNA bases)", 0, min = 0, max = 93),
        numericInput(ns("anchor"), "RegTools minimum anchor (bp)", 8, min = 1, max = 1000),
        numericInput(ns("min_intron"), "Minimum intron length (bp)", 70, min = 1),
        numericInput(ns("max_intron"), "Maximum intron length (bp)", 500000, min = 1),
        selectInput(ns("strand_mode"), "Library strandedness", c("XS tags" = "XS", "RF / first-strand" = "RF", "FR / second-strand" = "FR")),
        checkboxInput(ns("exclude_dup"), "Exclude marked duplicates", FALSE),
        checkboxInput(ns("nh1"), "Require NH == 1 (exclude missing NH)", FALSE),
        helpText("Always exclude unmapped, secondary, supplementary and QC-failed alignments. Paired ends count separately.")),
      bslib::input_task_button(ns("run"), "Compare genotypes"),
      hr(),
      h5("Saved full-comparison jobs"),
      selectInput(ns("saved_job"), "Your recent jobs", choices = character()),
      actionButton(ns("refresh_jobs"), "Refresh jobs"),
      actionButton(ns("open_job"), "Open saved job"),
      actionButton(ns("resume_job"), "Resume interrupted job"),
      helpText("Opening displays the job's submitted variant and parameters. Resume reuses verified sample checkpoints; completed results remain available. Search again to submit a new comparison."),
      hr(),
      helpText("DNA GT defines each group. Missing and partial calls remain in the audit and are excluded from RNA comparison. RNA coverage and junctions are descriptive evidence, not a formal splicing association test.")
    ),
    uiOutput(ns("lookup_status")),
    uiOutput(ns("job_progress")),
    uiOutput(ns("status")),
    uiOutput(ns("selected_variant")),
    bslib::navset_card_tab(
      bslib::nav_panel("Genotype comparison",
        h5("Observed genotype groups and sample accounting"),
        DT::DTOutput(ns("groups")),
        p("dna_n counts calls in the VCF; available_n counts calls linked to the RNA manifest. selected_n is the actual number requested, analyzed_n succeeded, and failed_n failed. Missing RNA links and incomplete DNA calls are separate, potentially overlapping exclusions. An unavailable BAM or index remains a recorded failure. A genotype absent from the VCF is not invented, and a missing call is not reference homozygous."),
        bslib::card(bslib::card_header("Mean RNA read depth per successful sample"), plotOutput(ns("depth_plot"), height = 320)),
        bslib::card(bslib::card_header("Junction support by genotype · top 20 per group"), uiOutput(ns("junction_chart"))),
        annotation_track_ui(ns("annotation")),
        p("Depth includes zero-coverage positions. Junction means include successful samples with zero support; failed samples are excluded and shown above. These raw counts are not normalized for library size. Arcs are not transcript annotations.")),
      bslib::nav_panel("Junctions", DT::DTOutput(ns("junctions")),
        p("Intron boundaries are 1-based and inclusive. count_sum is total retained alignment support; count_mean is support per successful RNA sample under the submitted sample scope.")),
      bslib::nav_panel("GENCODE v48", annotation_detail_ui(ns("annotation"))),
      bslib::nav_panel("RNA allele evidence",
        uiOutput(ns("rna_status")),
        splice_evidence_ui(ns("splice")),
        h5("RNA base observations at the selected variant"),
        DT::DTOutput(ns("rna_bases")),
        p("Inspect the exclusion columns: no valid base observations can reflect lack of overlap, a skipped intron / deletion, low base quality, or missing sequence / quality. Missing RNA base data is distinguished from an observed zero."),
        p("Base observations use retained RNA alignments at the selected SNV position; this is not a variant caller. They describe RNA evidence and do not replace the DNA genotype call. Zero observations can reflect no RNA coverage and do not establish absence of a DNA mutation. Indel / symbolic-allele base comparison is unavailable.")),
      bslib::nav_panel("DNA calls and sample audit",
        h5("Selected RNA samples · successes and failures"), DT::DTOutput(ns("samples")),
        h5("All DNA calls at the selected record"), DT::DTOutput(ns("genotypes")),
        p("raw_gt preserves the VCF call, including phasing and missing alleles. Comparison groups combine allele order / phase while retaining haploid, multiallelic and incomplete calls as distinct labels.")),
      bslib::nav_panel("Export and provenance",
        downloadButton(ns("download_groups"), "Genotype counts CSV"),
        downloadButton(ns("download_junctions"), "Junctions CSV"),
        downloadButton(ns("download_depth"), "Depth TSV"),
        downloadButton(ns("download_rna"), "RNA bases CSV"),
        downloadButton(ns("download_samples"), "Sample audit CSV"),
        downloadButton(ns("download_genotypes"), "DNA calls CSV"),
        downloadButton(ns("download_audit"), "Parameters / provenance JSON"),
        hr(), verbatimTextOutput(ns("audit")),
        p("Exports include controlled sample metadata. Use them within the project's authorized data environment.")),
      bslib::nav_panel("Scope",
        p("This extends Reg_Shiny from a fixed candidate list to indexed WGS VCF coordinate lookup with explicitly matched RNA BAMs. Only records present in the configured VCFs can be found. No record found does not establish a reference genotype or absence of variation."),
        p("The declared assembly must match across VCFs, BAMs and sample manifests. The app checks available contig lengths; it does not perform liftOver or prove sequence-level reference identity."),
        p("The default processes all matched RNA BAMs with complete DNA genotype calls in a persistent SCC job. Optional preview uses a deterministic bounded subset. Missing and partial DNA calls remain in the audit with unavailable RNA evidence. An interval is limited to 249,999 bp in this interface and 100,000 candidate alignments per sample. Limit failures remain failures; they never trigger silent sample reduction."),
        p("Use the shared project launcher through SCC's authenticated environment. The app does not supply its own authentication and must not be exposed directly to the public internet. Input VCF and BAM files are read only."),
        p("The Single BAM / demo tab preserves direct interval exploration and exact / near target-junction metrics. Neither mode provides formal PSI or causal inference."))
    )
  )
}

variant_server <- function(id, backend, resources, backend_file, variant_file,
                           job_api = NULL, job_source = NULL, job_root = NULL,
                           annotation_resources = NULL, annotation_file = NULL,
                           splice_reference = NULL, splice_file = NULL,
                           source_identity = NULL, is_active = reactive(TRUE)) {
  force(resources); force(source_identity)
  shiny::moduleServer(id, function(input, output, session) {
    resource_ok <- !inherits(resources, "error")
    output$resource_status <- renderUI({
      if (!resource_ok) return(div(class = "alert alert-warning",
        "WGS resources are unavailable: ", conditionMessage(resources),
        " The Single BAM / demo tab remains available."))
      tagList(if (!is.null(source_identity)) helpText(paste("Dataset:", source_identity$label, "·", source_identity$id)),
        strong(paste("Declared assembly:", resources$build)),
        helpText(paste("Configured chromosomes:", paste(resources$vcf_registry$chrom, collapse = ", "))),
        helpText(paste("Explicit RNA sample mappings:", nrow(resources$samples))),
        helpText("Reference identity is declared in the resource manifests; available contig lengths are checked during analysis."))
    })
    lookup_query <- reactiveVal(NULL)
    lookup_revision <- reactiveVal(0L)
    lookup_task <- shiny::ExtendedTask$new(function(res, query, source_file, variant_source) {
      promises::future_promise({
        e <- new.env(parent = globalenv())
        sys.source(source_file, envir = e)
        sys.source(variant_source, envir = e)
        e$lookup_variants(res, query)
      }, seed = TRUE)
    }) |> bslib::bind_task_button("lookup")
    observeEvent(input$lookup, {
      req(is_active())
      if (!resource_ok) {
        showNotification(conditionMessage(resources), type = "error", duration = 10)
        bslib::update_task_button("lookup", state = "ready", session = session)
        return()
      }
      query <- trimws(input$query)
      if (!nzchar(query)) {
        showNotification("Enter a variant coordinate such as chr1:123456.", type = "error")
        bslib::update_task_button("lookup", state = "ready", session = session)
        return()
      }
      lookup_query(query)
      lookup_revision(lookup_revision() + 1L)
      updateSelectInput(session, "record", choices = character())
      lookup_task$invoke(resources, query, backend_file, variant_file)
    })
    records <- reactive({req(is_active()); lookup_task$result()})
    observeEvent(lookup_task$status(), {
      req(is_active())
      if (lookup_task$status() != "success") return()
      rows <- records()
      labels <- if (nrow(rows)) paste0(rows$chrom, ":", rows$pos1, " ", rows$ref, " > ", rows$alt,
        " | ", rows$id, " | FILTER=", rows$filter) else character()
      # Multiple exact-position records require an explicit allele choice.
      values <- as.character(rows$record_id)
      choices <- if (length(values) > 1L) c("Select a REF / ALT record" = "", stats::setNames(values, labels)) else stats::setNames(values, labels)
      updateSelectInput(session, "record", choices = choices, selected = if (length(values) == 1L) values else "")
    }, ignoreInit = TRUE)
    output$lookup_status <- renderUI({
      s <- lookup_task$status()
      if (s == "initial") return(p("Search for a WGS coordinate to load the DNA variant record, then compare the matched RNA samples."))
      if (s == "running") return(div(class = "alert alert-info", "Looking up the submitted coordinate in indexed VCF data..."))
      if (s == "error") {
        msg <- tryCatch(lookup_task$result(), error = function(e) conditionMessage(e))
        return(div(class = "alert alert-danger", "Variant lookup failed: ", as.character(msg)))
      }
      rows <- records()
      if (!identical(trimws(input$query), lookup_query()))
        return(div(class = "alert alert-warning", "The coordinate input has changed. Run Find variant before comparing the new coordinate."))
      if (!nrow(rows)) return(div(class = "alert alert-warning", "No VCF record has this exact position / allele combination. This is not evidence of a reference genotype or absence of variation."))
      div(class = "alert alert-info", nrow(rows), " matching record(s). ",
        if (nrow(rows) > 1L) "Choose the intended REF / ALT record explicitly." else "The exact REF / ALT record is selected.")
    })
    selected_record <- reactive({
      rows <- records()
      req(length(input$record) == 1L, nzchar(input$record))
      out <- rows[as.character(rows$record_id) == input$record, , drop = FALSE]
      validate(need(nrow(out) == 1L, "Select exactly one REF / ALT record."))
      out
    })
    controls <- reactive(list(query = trimws(input$query), record_id = input$record, lookup_revision = lookup_revision(),
      data_source_id = if (!is.null(source_identity)) source_identity$id else NULL,
      sampling_mode = if (is.null(input$sampling_mode)) "all" else input$sampling_mode,
      flank = input$flank, sample_cap = if (identical(input$sampling_mode, "preview")) input$sample_cap else NULL, mapq = input$mapq,
      baseq = input$baseq, anchor = input$anchor, min_intron = input$min_intron,
      max_intron = input$max_intron, strand_mode = input$strand_mode,
      exclude_duplicates = input$exclude_dup, nh1_only = input$nh1))
    launched_controls <- reactiveVal(NULL)
    launched_config <- reactiveVal(NULL)
    active_mode <- reactiveVal(NULL)
    active_job <- reactiveVal(NULL)
    active_job_variant <- reactiveVal(NULL)
    restored_job <- reactiveVal(FALSE)
    full_result <- reactiveVal(NULL)
    full_result_error <- reactiveVal(NULL)
    loaded_job <- reactiveVal(NULL)
    observeEvent(is_active(), {
      if (!isTRUE(is_active())) {
        # Release loaded full results and stop polling; the scheduler job itself
        # remains under its original immutable dataset and is never cancelled.
        full_result(NULL); full_result_error(NULL); active_job(NULL); loaded_job(NULL)
      }
    })
    results_stale <- reactive(!identical(controls(), launched_controls()))
    analysis_task <- shiny::ExtendedTask$new(function(res, variant, cfg, cap, source_file, variant_source) {
      promises::future_promise({
        e <- new.env(parent = globalenv())
        sys.source(source_file, envir = e)
        sys.source(variant_source, envir = e)
        e$analyze_variant(res, variant, cfg, max_per_group = cap, sampling_mode = "preview")
      }, seed = TRUE)
    }) |> bslib::bind_task_button("run")
    observeEvent(input$run, {
      req(is_active())
      problem <- tryCatch({
        if (!resource_ok) stop(conditionMessage(resources))
        if (lookup_task$status() != "success" || !identical(trimws(input$query), lookup_query()))
          stop("Run Find variant for the current coordinate first.")
        rows <- records()
        if (length(input$record) != 1L || !nzchar(input$record)) stop("Select a REF / ALT record first.")
        variant <- rows[as.character(rows$record_id) == input$record, , drop = FALSE]
        if (nrow(variant) != 1L) stop("Select exactly one returned REF / ALT record.")
        flank <- backend$int_scalar(input$flank, "Flank", 50, 124999)
        mode <- if (is.null(input$sampling_mode)) "all" else input$sampling_mode
        if (!mode %in% c("all", "preview")) stop("Select All matched BAMs or the optional preview.")
        cap <- if (identical(mode, "preview")) backend$int_scalar(input$sample_cap, "Samples per genotype", 1, 10) else NULL
        end1 <- variant$pos1[[1L]] + flank
        if ("contig_length" %in% names(variant) && is.finite(variant$contig_length[[1L]]))
          end1 <- min(end1, variant$contig_length[[1L]])
        cfg <- list(demo = FALSE, bam = "", chrom = variant$chrom[[1L]],
          start1 = max(1, variant$pos1[[1L]] - flank), end1 = end1,
          mapq = input$mapq, baseq = input$baseq, anchor = input$anchor,
          min_intron = input$min_intron, max_intron = input$max_intron,
          strand_mode = input$strand_mode, exclude_duplicates = input$exclude_dup, nh1_only = input$nh1)
        cfg <- backend$validate_config(cfg)
        if (!is.null(source_identity)) {
          cfg$data_source_id <- source_identity$id
          cfg$data_source_label <- source_identity$label
        }
        if (identical(mode, "all")) {
          if (is.null(job_api)) stop("Persistent SCC full-comparison jobs are not configured. No preview or sampling has been substituted.")
          receipt <- job_api$submit_variant_job(resources, variant, cfg,
            backend_file = backend_file, variant_file = variant_file,
            job_source = job_source, job_root = job_root, ui_controls = controls())
          if (is.null(receipt$job_dir) || !nzchar(receipt$job_dir)) stop("The SCC submission did not return a saved job directory.")
          launched_controls(controls()); launched_config(cfg); restored_job(FALSE)
          active_mode("all")
          active_job(receipt$job_dir)
          active_job_variant(variant)
          full_result(NULL); full_result_error(NULL); loaded_job(NULL)
          refresh_job_list()
          bslib::update_task_button("run", state = "ready", session = session)
        } else {
          launched_controls(controls()); launched_config(cfg); restored_job(FALSE)
          active_mode("preview")
          active_job(NULL)
          analysis_task$invoke(resources, variant, cfg, cap, backend_file, variant_file)
        }
        NULL
      }, error = function(e) conditionMessage(e))
      if (!is.null(problem)) {
        showNotification(problem, type = "error", duration = 10)
        bslib::update_task_button("run", state = "ready", session = session)
      }
    })
    job_list <- reactiveVal(data.frame(job_dir = character(), label = character()))
    refresh_job_list <- function() {
      if (is.null(job_api) || !isTRUE(isolate(is_active()))) return(invisible(NULL))
      found <- tryCatch(job_api$list_variant_jobs(job_root), error = function(e) e)
      if (inherits(found, "error")) {
        showNotification(paste("Cannot list saved jobs:", conditionMessage(found)), type = "error")
        return(invisible(NULL))
      }
      job_list(found)
      choices <- if (nrow(found)) stats::setNames(found$job_dir, found$label) else character()
      updateSelectInput(session, "saved_job", choices = choices,
        selected = if (!is.null(active_job()) && active_job() %in% found$job_dir) active_job() else NULL)
      invisible(found)
    }
    observeEvent(input$refresh_jobs, refresh_job_list())
    observeEvent(TRUE, refresh_job_list(), once = TRUE)
    poll_cache <- new.env(parent = emptyenv())
    poll_cache$value <- NULL
    job_state <- reactivePoll(2000, session,
      checkFunc = function() {
        path <- active_job()
        if (!isTRUE(is_active()) || is.null(path) || is.null(job_api)) { poll_cache$value <- NULL; return(NULL) }
        state <- tryCatch(job_api$read_variant_job(path, job_root = job_root), error = function(e)
          list(status = "STATUS_UNAVAILABLE", error = conditionMessage(e), result_ready = FALSE))
        poll_cache$value <- state
        jsonlite::toJSON(state, auto_unbox = TRUE, null = "null", na = "null")
      },
      valueFunc = function() poll_cache$value)
    open_saved_job <- function(path) {
      req(is_active())
      if (is.null(job_api)) stop("Persistent SCC jobs are not configured.")
      if (length(path) != 1L || !nzchar(path) || !path %in% job_list()$job_dir)
        stop("Select one of your listed jobs first.")
      state <- job_api$read_variant_job(path, job_root = job_root)
      if (!is.null(source_identity) && !identical(state$request$cfg$data_source_id, source_identity$id))
        stop("This saved job belongs to a different or unidentified dataset. Reconnect its original prepared profile before opening it.")
      active_mode("all"); active_job(path); restored_job(TRUE)
      full_result(NULL); full_result_error(NULL); loaded_job(NULL)
      launched_controls(state$request$ui_controls)
      if (!is.null(state$request$cfg)) launched_config(state$request$cfg)
      if (!is.null(state$request$variant)) active_job_variant(state$request$variant)
      # Reopening the same completed directory need not change the poll token.
      # Load it explicitly through the manager's validation gate in that case.
      if (isTRUE(state$result_ready)) {
        answer <- job_api$load_variant_job_result(path, job_root = job_root)
        full_result(answer); loaded_job(path)
        active_job_variant(answer$variant)
        if (!is.null(answer$provenance$config)) launched_config(answer$provenance$config)
      }
    }
    observeEvent(input$open_job, {
      req(is_active())
      tryCatch(open_saved_job(input$saved_job), error = function(e)
        showNotification(conditionMessage(e), type = "error", duration = 10))
    })
    observeEvent(input$resume_job, {
      req(is_active())
      tryCatch({
        path <- input$saved_job
        if (is.null(job_api) || length(path) != 1L || !path %in% job_list()$job_dir)
          stop("Select one of your listed jobs first.")
        if (!is.null(source_identity)) {
          state <- job_api$read_variant_job(path, job_root = job_root)
          if (!identical(state$request$cfg$data_source_id, source_identity$id))
            stop("This saved job belongs to a different or unidentified dataset. Reconnect its original prepared profile before resuming it.")
        }
        job_api$resume_variant_job(path, job_root = job_root)
        open_saved_job(path)
        refresh_job_list()
      }, error = function(e) showNotification(conditionMessage(e), type = "error", duration = 10))
    })
    observeEvent(job_state(), {
      req(is_active())
      state <- job_state()
      if (is.null(state)) return()
      if (!isTRUE(state$result_ready)) {
        full_result(NULL); loaded_job(NULL)
        return()
      }
      if (identical(loaded_job(), active_job())) return()
      answer <- tryCatch(job_api$load_variant_job_result(active_job(), job_root = job_root), error = function(e) e)
      if (inherits(answer, "error")) {
        full_result_error(conditionMessage(answer))
        return()
      }
      full_result(answer); full_result_error(NULL); loaded_job(active_job())
      active_job_variant(answer$variant)
      if (!is.null(answer$provenance$config)) launched_config(answer$provenance$config)
    }, ignoreNULL = TRUE)
    output$job_progress <- renderUI({
      if (!identical(active_mode(), "all")) return(NULL)
      state <- job_state()
      if (is.null(state)) return(div(class = "alert alert-info", "Reading the submitted SCC job..."))
      p <- state$progress
      num <- function(x) if (is.null(x) || length(x) != 1L || is.na(x)) "pending" else format(x, big.mark = ",", trim = TRUE)
      v <- active_job_variant()
      tagList(div(class = "alert alert-info",
        strong(paste0("All matched BAMs · ", state$status)),
        if (!is.null(state$job_id)) paste0(" · SCC job ", state$job_id),
        tags$br(), if (!is.null(v)) paste(v$chrom[[1L]], v$pos1[[1L]], v$ref[[1L]], v$alt[[1L]], sep = ":"),
        tags$br(), "Selected: ", num(p$total), " · completed: ", num(p$completed),
        " · succeeded: ", num(p$success), " · failed: ", num(p$failed),
        " · restored checkpoints: ", num(p$resumed),
        tags$br(), "Missing RNA link: ", num(p$missing_link_n),
        " · incomplete DNA call: ", num(p$incomplete_gt_n), " (exclusions may overlap)",
        tags$br(), "The saved SCC job continues independently of this browser. Results appear only after validation.",
        if (restored_job()) tags$p("Showing a saved job's submitted snapshot. The current search controls do not alter this job.") else
          if (results_stale()) tags$p("Current inputs differ from the submitted job. The running job retains its original variant and parameters.")),
        if (!is.null(p$group_summary) && is.null(full_result()))
          tagList(h5("Actual sample counts by genotype · job progress"), DT::DTOutput(session$ns("progress_groups"))),
        if (!is.null(state$error)) div(class = "alert alert-danger", state$error),
        if (!is.null(full_result_error())) div(class = "alert alert-danger", "Result validation failed: ", full_result_error()))
    })
    output$progress_groups <- DT::renderDT({
      req(identical(active_mode(), "all"))
      g <- job_state()$progress$group_summary
      req(!is.null(g))
      g
    }, rownames = FALSE, selection = "none", options = list(dom = "t", paging = FALSE, scrollX = TRUE))
    result <- reactive({
      req(is_active())
      if (identical(active_mode(), "all")) {
        req(!is.null(full_result()))
        full_result()
      } else analysis_task$result()
    })
    result_mode <- reactive({
      mode <- result()$provenance$selection$mode
      if (is.null(mode)) active_mode() else mode
    })
    annotation_server("annotation", resources = annotation_resources, annotation_file = annotation_file,
      region = reactive({
        r <- result(); cfg <- launched_config()
        list(chrom = cfg$chrom, start1 = cfg$start1, end1 = cfg$end1,
          build = r$variant$build[[1L]], demo = FALSE,
          build_basis = "Reference assembly declared by the matched DNA/RNA manifests.")
      }),
      display_range = reactive(c(launched_config()$start1, launched_config()$end1)),
      junctions = reactive(result()$junctions), marker = reactive(result()$variant$pos1[[1L]]))
    output$status <- renderUI({
      if (identical(active_mode(), "all")) {
        if (is.null(full_result())) return(NULL)
      } else {
        s <- analysis_task$status()
        if (s == "initial") return(NULL)
        if (s == "running") return(div(class = "alert alert-info", "Analyzing the submitted exploratory preview. Samples run sequentially within this task; results appear after the audit completes."))
        if (s == "error") {
          msg <- tryCatch(analysis_task$result(), error = function(e) conditionMessage(e))
          return(div(class = "alert alert-danger", "Genotype comparison failed: ", as.character(msg)))
        }
      }
      r <- result()
      stale <- !restored_job() && results_stale()
      any_failed <- any(r$group_summary$failed_n > 0, na.rm = TRUE)
      any_success <- any(r$group_summary$analyzed_n > 0, na.rm = TRUE)
      full <- identical(result_mode(), "all")
      tagList(div(class = if (stale || any_failed || !any_success) "alert alert-warning" else "alert alert-success",
        if (stale) "Inputs have changed. Displayed results and downloads belong to the previous submitted variant / parameter snapshot. " else "",
        if (restored_job()) "Displaying the saved job's submitted variant and parameters. ",
        if (!any_success) "No RNA sample completed successfully; inspect the sample audit. " else
          if (full && any_failed) "All matched BAMs were attempted; RNA evidence is partial because some samples failed. " else
            if (full) "All matched BAMs completed successfully. " else "Exploratory subset analysis completed. ",
        if (any_failed) "Some selected samples failed and are excluded from group means; failures remain in the sample audit."),
        if (length(r$warnings)) div(class = "alert alert-warning", paste(unique(r$warnings), collapse = "; ")))
    })
    output$selected_variant <- renderUI({
      r <- result(); v <- r$variant
      div(class = "alert alert-secondary", strong(paste0(v$chrom[[1L]], ":", v$pos1[[1L]], " · ", v$ref[[1L]], " > ", v$alt[[1L]])),
        if (!is.null(source_identity)) paste0(" · dataset ", source_identity$label, " [", source_identity$id, "]"),
        " · ", v$build[[1L]], " · FILTER=", v$filter[[1L]],
        " · RNA interval ", launched_config()$start1, "–", launched_config()$end1,
        if (identical(result_mode(), "all")) " · ALL MATCHED BAMs with complete DNA calls" else
          paste0(" · cap ", launched_controls()$sample_cap, " per observed genotype · exploratory subset"))
    })
    display_groups <- reactive({
      g <- result()$group_summary
      g$missing_rna_link_n <- g$dna_n - g$available_n
      g$incomplete_dna_call_n <- ifelse(g$call_status == "CALLED", 0L, g$dna_n)
      g
    })
    output$groups <- DT::renderDT(display_groups(), rownames = FALSE, selection = "none",
      options = list(dom = "t", scrollX = TRUE, paging = FALSE))
    for (name in c("junctions", "samples", "genotypes", "rna_bases")) local({
      field <- name
      output[[field]] <- DT::renderDT(result()[[field]], rownames = FALSE, selection = "none",
        options = list(pageLength = 10, scrollX = TRUE))
    })
    group_colors <- reactive({
      groups <- as.character(result()$group_summary$genotype)
      stats::setNames(grDevices::hcl.colors(max(1L, length(groups)), palette = "Dark 3"), groups)
    })
    output$depth_plot <- renderPlot({
      r <- result(); d <- r$depth; cfg <- launched_config(); cols <- group_colors()
      finite <- if (nrow(d)) is.finite(d$mean_depth) else logical()
      plot(NA, xlim = c(cfg$start1, cfg$end1), ylim = c(0, if (any(finite)) max(1, d$mean_depth[finite]) else 1),
        xlab = paste0(r$variant$chrom[[1L]], " · 1-based position"), ylab = "Mean read depth / successful sample")
      if (!any(finite)) { text(mean(c(cfg$start1, cfg$end1)), .5, "No successfully analyzed RNA depth"); return(invisible(NULL)) }
      for (g in unique(d$genotype)) {
        x <- d[d$genotype == g, , drop = FALSE]; x <- x[order(x$pos1), , drop = FALSE]
        lines(x$pos1, x$mean_depth, col = cols[[g]], lwd = 1.8)
      }
      abline(v = r$variant$pos1[[1L]], lty = 2, col = "gray40")
      shown <- unique(as.character(d$genotype[finite]))
      legend("topright", legend = shown, col = cols[shown], lty = 1, lwd = 2, bty = "n", cex = .85)
    })
    output$junction_chart <- renderUI({
      n <- max(1L, nrow(result()$group_summary))
      plotOutput(session$ns("junction_plot"), height = max(260, 170 * n))
    })
    output$junction_plot <- renderPlot({
      r <- result(); gs <- r$group_summary; cfg <- launched_config(); cols <- group_colors()
      if (!nrow(gs)) { plot.new(); text(.5, .5, "No linked RNA genotype groups"); return(invisible(NULL)) }
      old <- par(mfrow = c(nrow(gs), 1L), mar = c(2.8, 4.1, 2.2, 2.1))
      on.exit(par(old))
      for (i in seq_len(nrow(gs))) {
        g <- as.character(gs$genotype[[i]])
        plot(NA, xlim = c(cfg$start1, cfg$end1), ylim = c(0, 1.2), yaxt = "n", ylab = g,
          xlab = "1-based intron positions", main = paste0("Successful RNA samples: ", gs$analyzed_n[[i]], "; failed: ", gs$failed_n[[i]]))
        abline(v = r$variant$pos1[[1L]], lty = 2, col = "gray70")
        j <- r$junctions
        if (nrow(j)) j <- j[j$genotype == g & j$intron_start1 >= cfg$start1 & j$intron_end1 <= cfg$end1 & is.finite(j$count_mean) & j$count_sum > 0, , drop = FALSE]
        if (!nrow(j)) {
          text(mean(c(cfg$start1, cfg$end1)), .5,
            if (gs$analyzed_n[[i]] == 0) "Unavailable: no successful RNA samples" else "No supported junctions fully inside this interval")
          next
        }
        j <- head(j[order(-j$count_sum, j$intron_start1), , drop = FALSE], 20L)
        for (k in seq_len(nrow(j))) {
          t <- seq(0, 1, length.out = 60L)
          h <- .15 + .65 * j$count_mean[[k]] / max(1, r$junctions$count_mean, na.rm = TRUE)
          x <- j$intron_start1[[k]] + t * (j$intron_end1[[k]] - j$intron_start1[[k]])
          lines(x, 4 * h * t * (1 - t), col = cols[[g]], lwd = 1.6)
          text(mean(range(x)), h + .08, paste0(format(round(j$count_mean[[k]], 2), trim = TRUE), " / sample ", j$strand[[k]]), cex = .7)
        }
      }
    })
    output$rna_status <- renderUI({
      r <- result(); v <- r$variant
      is_snv <- identical(r$provenance$rna_evidence_status, "SNV_BASE_COUNTS")
      div(class = if (is_snv) "alert alert-info" else "alert alert-warning",
        if (is_snv) paste0("DNA REF: ", v$ref[[1L]], "; ALT: ", v$alt[[1L]], ". Tracks show mean RNA depth per successful sample; the base-observation table pools counts within each DNA genotype group.") else
          "RNA base comparison is unavailable for this non-SNV record. DNA calls, junctions and read depth remain available.")
    })
    splice_evidence_server("splice", result = result, reference = splice_reference,
      model_source = splice_file, annotation_resources = annotation_resources,
      annotation_source = annotation_file)
    output$audit <- renderPrint({
      r <- result()
      cat(if (identical(result_mode(), "all")) "All matched BAMs with complete DNA calls\n" else "Exploratory deterministic sample subset\n")
      cat("Submitted variant:", paste(r$variant$chrom, r$variant$pos1, r$variant$ref, r$variant$alt, sep = ":"), "\n")
      cat("Declared assembly:", r$variant$build[[1L]], "\n")
      cat("Successful samples:", sum(r$group_summary$analyzed_n), "; failed samples:", sum(r$group_summary$failed_n), "\n\n")
      print(r$provenance)
    })
    csv_export <- function(field, filename) downloadHandler(filename, function(file)
      utils::write.csv(result()[[field]], file, row.names = FALSE, na = "NA", fileEncoding = "UTF-8"))
    output$download_groups <- downloadHandler("genotype_sample_counts.csv", function(file)
      utils::write.csv(display_groups(), file, row.names = FALSE, na = "NA", fileEncoding = "UTF-8"))
    output$download_junctions <- csv_export("junctions", "genotype_junctions.csv")
    output$download_rna <- csv_export("rna_bases", "genotype_rna_bases.csv")
    output$download_samples <- csv_export("samples", "selected_rna_sample_audit.csv")
    output$download_genotypes <- csv_export("genotypes", "dna_genotype_calls.csv")
    output$download_depth <- downloadHandler("genotype_depth.tsv", function(file)
      utils::write.table(result()$depth, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"))
    output$download_audit <- downloadHandler("variant_parameters_and_provenance.json", function(file) {
      r <- result()
      sample_evidence <- lapply(seq_len(nrow(r$samples)), function(i) {
        z <- if (length(r$sample_results) >= i) r$sample_results[[i]] else NULL
        list(sample = as.list(r$samples[i, , drop = FALSE]),
          analysis = if (is.null(z)) NULL else list(config = z$config,
            candidate_alignments = z$candidate_alignments, retained_alignments = r$samples$retained_reads[[i]],
            span_only_excluded = z$span_only_excluded, native_audit_passed = z$native_audit_passed,
            versions = z$versions[setdiff(names(z$versions), "regtools_extract_help")],
            source_metadata = z$source_metadata, completed_at = z$completed_at,
            rna_bases = z$rna_bases, rna_base_exclusions = as.list(attr(z$rna_bases, "excluded")),
            warnings = z$warnings))
      })
      x <- list(variant = r$variant, config = launched_config(), controls = launched_controls(),
        group_summary = display_groups(), selected_samples = r$samples, provenance = r$provenance,
        sample_evidence = sample_evidence,
        warnings = r$warnings, definitions = list(selection = if (identical(result_mode(), "all")) "all linked RNA BAMs with complete DNA genotype calls" else "exploratory deterministic first n sorted sample IDs per observed genotype",
          dna_genotype = "VCF GT; raw phasing and missingness retained in DNA calls export",
          depth = "mean per successful sample at each position, including zero coverage",
          junction_mean = "retained alignment support / successful samples, including zero support",
          rna_bases = "retained RNA aligned-base observations at SNV; unavailable for other variant types",
          limitations = "descriptive RNA evidence; not library-size normalized, formal PSI or causal inference; failed and missing samples retain their states"))
      jsonlite::write_json(x, file, auto_unbox = TRUE, pretty = TRUE, na = "null")
    })
    invisible(list(result = result, lookup_task = lookup_task, analysis_task = analysis_task,
      results_stale = results_stale, selected_record = selected_record,
      active_mode = active_mode, active_job = active_job, job_state = job_state,
      open_saved_job = open_saved_job, full_result = full_result))
  })
}

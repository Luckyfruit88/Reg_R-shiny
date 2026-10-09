# Variant-centered interface inside Reg_Shiny. Preview workers share the app's
# future pool. Full comparisons use persistent SCC jobs and survive disconnects.

# Resolve each requested record independently, then issue a separate persistent
# job submission. This orchestrates metadata work only, never a shared RNA job.
variant_batch_submit <- function(api, resources, lines, settings, source_identity,
                                 backend_file, variant_file, job_source, job_root, controls) {
  outcomes <- list(); submitted_records <- character()
  for (i in seq_len(nrow(lines))) {
    row <- list(input_line = lines$input_line[[i]], query = lines$query[[i]], record_id = NA_character_,
      outcome = "LOOKUP_FAILED", job_id = NA_character_, job_dir = NA_character_, requested_cores = NA_integer_,
      candidates = "", message = "")
    records <- tryCatch(api$lookup_variants(resources, row$query), error = function(e) e)
    if (inherits(records, "error")) row$message <- conditionMessage(records) else if (!nrow(records)) {
      row$outcome <- "NOT_FOUND"; row$message <- "No exact VCF record; a reference genotype is not inferred."
    } else if (nrow(records) != 1L) {
      row$outcome <- "AMBIGUOUS_ALLELES"
      row$candidates <- paste(paste(records$chrom, records$pos1, records$ref, records$alt, sep = ":"), collapse = "; ")
      row$message <- "No job submitted. Specify one exact CHROM:POS:REF:ALT record."
    } else {
      variant <- records
      row$record_id <- as.character(variant$record_id[[1L]])
      if (row$record_id %in% submitted_records) {
        row$outcome <- "DUPLICATE_IN_BATCH"
        row$message <- "This exact record was already attempted in this batch; no duplicate submission."
      } else {
        # Mark attempted even when qsub is uncertain: never automatically retry
        # a possibly accepted submission later in the same batch.
        submitted_records <- c(submitted_records, row$record_id)
        attempt <- tryCatch({
          cfg <- settings
          flank <- cfg$flank; cfg$flank <- NULL
          cfg$chrom <- variant$chrom[[1L]]
          cfg$start1 <- max(1, variant$pos1[[1L]] - flank)
          cfg$end1 <- variant$pos1[[1L]] + flank
          if ("contig_length" %in% names(variant) && is.finite(variant$contig_length[[1L]]))
            cfg$end1 <- min(cfg$end1, variant$contig_length[[1L]])
          cfg <- api$validate_config(cfg)
          if (!is.null(source_identity)) {
            cfg$data_source_id <- source_identity$id; cfg$data_source_label <- source_identity$label
          }
          frozen <- controls
          frozen$query <- row$query; frozen$record_id <- row$record_id; frozen$sampling_mode <- "all"
          frozen$sample_cap <- NULL; frozen$batch_input_line <- row$input_line
          api$submit_variant_job(resources, variant, cfg, backend_file = backend_file,
            variant_file = variant_file, job_source = job_source, job_root = job_root, ui_controls = frozen)
        }, error = function(e) e)
        if (inherits(attempt, "error")) {
          row$outcome <- "SUBMISSION_FAILED"
          row$message <- conditionMessage(attempt)
          if (length(attempt$job_dir) == 1L) row$job_dir <- as.character(attempt$job_dir)
          if (length(attempt$job_id) == 1L) row$job_id <- as.character(attempt$job_id)
          if (length(attempt$status) == 1L) row$outcome <- as.character(attempt$status)
        } else {
          row$outcome <- "SUBMITTED"; row$job_id <- as.character(attempt$job_id)
          row$job_dir <- as.character(attempt$job_dir); row$requested_cores <- 16L
          row$message <- "Independent full-cohort job submitted; inspect its own status and logs."
        }
      }
    }
    outcomes[[i]] <- as.data.frame(row, stringsAsFactors = FALSE)
  }
  do.call(rbind, outcomes)
}

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
        helpText("All matched BAMs with complete DNA genotype calls are processed. The SCC job and sample checkpoints persist after the browser is closed; open the Jobs tab to see progress or results.")),
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
      actionButton(ns("show_jobs"), "Batch submissions and saved jobs"),
      helpText("Use Jobs to review batch resources, reopen a result, or resume an interrupted job."),
      helpText("DNA GT defines each group. Missing and partial calls remain in the audit and are excluded from RNA comparison. RNA coverage and junctions are descriptive evidence, not a formal splicing association test.")
    ),
    uiOutput(ns("view_context"), role = "status", `aria-live` = "polite"),
    uiOutput(ns("lookup_status"), role = "status", `aria-live` = "polite"),
    uiOutput(ns("submission_status")),
    uiOutput(ns("resume_status")),
    uiOutput(ns("job_progress")),
    uiOutput(ns("status"), role = "status", `aria-live` = "polite"),
    uiOutput(ns("selected_variant")),
    bslib::navset_card_tab(id = ns("result_tabs"),
      bslib::nav_panel("Genotype comparison",
        h5("Observed genotype groups and sample accounting"),
        DT::DTOutput(ns("groups")),
        tags$details(tags$summary("Sample accounting definitions and exclusions"),
          p("dna_n counts calls in the VCF; available_n counts calls linked to the RNA manifest. selected_n is the actual number requested, analyzed_n succeeded, and failed_n failed. Missing RNA links and incomplete DNA calls are separate, potentially overlapping exclusions. An unavailable BAM or index remains a recorded failure. A genotype absent from the VCF is not invented, and a missing call is not reference homozygous.")),
        bslib::card(bslib::card_header("Mean RNA read depth per successful sample"), plotOutput(ns("depth_plot"), height = 320)),
        bslib::card(bslib::card_header("Junction support by genotype · top 20 per group"), uiOutput(ns("junction_chart"))),
        annotation_track_ui(ns("annotation")),
        p("Depth includes zero-coverage positions. Junction means include successful samples with zero support; failed samples are excluded and shown above. These raw counts are not normalized for library size. Arcs are not transcript annotations.")),
      bslib::nav_panel("Jobs",
        bslib::layout_columns(col_widths = c(6, 6),
          bslib::card(bslib::card_header("Batch submission"),
            textAreaInput(ns("batch_queries"), "Batch variants (one per line)", rows = 5,
              placeholder = "chr21:14288395:G:A\nchr17:3156517:T:C"),
            bslib::input_task_button(ns("batch_submit"), "Review batch and resources"),
            uiOutput(ns("batch_review_status"), role = "status", `aria-live` = "polite"),
            helpText("Batch jobs always use all matched BAMs, even when the single-variant form is set to preview. The current interval and filters apply to every batch entry. Review the maximum resource request before confirming; no jobs are submitted by the review button.")),
          bslib::card(bslib::card_header("Saved full-comparison jobs"),
            selectInput(ns("saved_job"), "Your recent jobs", choices = character()),
            div(class = "regshiny-actions",
              actionButton(ns("refresh_jobs"), "Refresh jobs"),
              actionButton(ns("open_job"), "Open saved job"),
              bslib::input_task_button(ns("resume_job"), "Resume interrupted job")),
            helpText("Opening shows the job's frozen variant and parameters. Resume reuses verified checkpoints; editing the analysis form never changes a saved job."))),
        h5("Independent full-cohort jobs"),
        p("Every new full job requests 16 cores. Several jobs may run concurrently when the scheduler allocates their separate resources. Select a row and open that job's frozen result."),
        DT::DTOutput(ns("jobs_table")), actionButton(ns("open_dashboard_job"), "Open selected job"),
        h5("Latest batch outcomes"), uiOutput(ns("batch_status")), DT::DTOutput(ns("batch_outcomes")),
        downloadButton(ns("download_batch"), "Batch outcomes CSV"),
        helpText("Partial failures are retained per input line. Resolve ambiguous alleles before resubmitting only those lines. Inspect an uncertain submission's own receipt before retrying; successful jobs are never rolled back.")),
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
        hr(), uiOutput(ns("job_files")), verbatimTextOutput(ns("audit")),
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
    view_revision <- reactiveVal(0L)
    single_revision <- reactiveVal(0L)
    single_submission <- reactiveVal(NULL)
    single_receipt <- reactiveVal(NULL)
    observeEvent(is_active(), {
      if (!isTRUE(is_active())) {
        # Release loaded full results and stop polling; the scheduler job itself
        # remains under its original immutable dataset and is never cancelled.
        full_result(NULL); full_result_error(NULL); active_job(NULL); loaded_job(NULL)
      }
    })
    results_stale <- reactive(!identical(controls(), launched_controls()))
    single_submit_task <- shiny::ExtendedTask$new(function(res, variant, cfg, frozen_controls,
      source_file, variant_source, jobs_source, root, serial, view_token) {
      promises::future_promise({
        outcome <- tryCatch({
          e <- new.env(parent = globalenv())
          sys.source(source_file, envir = e); sys.source(variant_source, envir = e); sys.source(jobs_source, envir = e)
          receipt <- e$submit_variant_job(res, variant, cfg, backend_file = source_file,
            variant_file = variant_source, job_source = jobs_source, job_root = root, ui_controls = frozen_controls)
          if (length(receipt$job_dir) != 1L || !nzchar(receipt$job_dir)) stop("Submission did not return an independent job directory.")
          list(ok = TRUE, receipt = receipt)
        }, error = function(e) list(ok = FALSE, message = conditionMessage(e),
          job_dir = e$job_dir, status = e$status))
        c(outcome, list(variant = variant, cfg = cfg, controls = frozen_controls,
          revision = serial, view_token = view_token))
      }, seed = TRUE)
    })
    observeEvent(single_submit_task$status(), {
      req(is_active())
      status <- single_submit_task$status()
      if (!status %in% c("success", "error")) return()
      bslib::update_task_button("run", state = "ready", session = session)
      if (status == "error") return()
      value <- single_submit_task$result()
      if (!identical(value$revision, single_revision())) return()
      single_receipt(value)
      refresh_job_list()
      # An A receipt cannot replace a B query, a manually opened result, or a
      # newer source/view intent that changed while matching/qsub was running.
      if (isTRUE(value$ok) && identical(controls(), value$controls) && identical(view_revision(), value$view_token)) {
        launched_controls(value$controls); launched_config(value$cfg); restored_job(FALSE)
        active_mode("all"); active_job(value$receipt$job_dir); active_job_variant(value$variant)
        full_result(NULL); full_result_error(NULL); loaded_job(NULL)
      }
    }, ignoreInit = TRUE)
    output$submission_status <- renderUI({
      pending <- single_submission()
      if (is.null(pending)) return(NULL)
      s <- single_submit_task$status()
      label <- paste(pending$variant$chrom, pending$variant$pos1, pending$variant$ref, pending$variant$alt, sep = ":")
      if (s == "running") return(div(class = "alert alert-info", "Submitting ", label,
        " as an independent 16-core full job. Its allele record and filters are frozen; other lookups remain available."))
      if (s == "error") {
        message <- tryCatch(single_submit_task$result(), error = function(e) conditionMessage(e))
        return(div(class = "alert alert-warning", "Submission failed for ", label, ": ", as.character(message)))
      }
      receipt <- single_receipt()
      if (is.null(receipt)) return(NULL)
      if (!isTRUE(receipt$ok)) return(div(class = "alert alert-warning", "Submission failed for ", label, ": ", receipt$message,
        if (!is.null(receipt$job_dir)) paste0(" Inspect its own receipt/logs: ", receipt$job_dir)))
      div(class = "alert alert-success", "Submitted ", label, " · job ", receipt$receipt$job_id,
        " · 16 cores requested for this variant. Open it from Jobs if another query or result is now selected.")
    })
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
          single_revision(single_revision() + 1L); view_revision(view_revision() + 1L)
          single_submission(list(variant = variant, cfg = cfg, controls = controls(), revision = single_revision()))
          single_receipt(NULL)
          bslib::update_task_button("run", state = "busy", session = session)
          single_submit_task$invoke(resources, variant, cfg, controls(), backend_file, variant_file,
            job_source, job_root, single_revision(), view_revision())
        } else {
          view_revision(view_revision() + 1L)
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
    dashboard_revision <- reactiveVal(0L)
    dashboard_selected <- reactiveVal(NULL)
    selected_logs <- reactiveVal(NULL)
    batch_snapshot <- reactiveVal(NULL)
    batch_task <- shiny::ExtendedTask$new(function(res, lines, settings, identity, frozen_controls,
      source_file, variant_source, jobs_source, root) {
      promises::future_promise({
        e <- new.env(parent = globalenv())
        sys.source(source_file, envir = e); sys.source(variant_source, envir = e); sys.source(jobs_source, envir = e)
        variant_batch_submit(e, res, lines, settings, identity, source_file, variant_source, jobs_source, root, frozen_controls)
      }, seed = TRUE)
    }) |> bslib::bind_task_button("batch_submit")
    batch_review <- reactiveVal(NULL)
    batch_review_error <- reactiveVal(NULL)
    batch_inputs <- function() {
      if (!resource_ok) stop(conditionMessage(resources))
      if (is.null(job_api)) stop("Persistent SCC jobs are not configured.")
      raw <- strsplit(if (is.null(input$batch_queries)) "" else input$batch_queries, "\n", fixed = TRUE)[[1L]]
      rows <- data.frame(input_line = seq_along(raw), query = trimws(raw), stringsAsFactors = FALSE)
      rows <- rows[nzchar(rows$query), , drop = FALSE]
      if (!nrow(rows)) stop("Enter at least one variant coordinate, one per line.")
      flank <- backend$int_scalar(input$flank, "Flank", 50, 124999)
      settings <- list(demo = FALSE, bam = "", chrom = "validation", start1 = 1L, end1 = 101L,
        mapq = input$mapq, baseq = input$baseq, anchor = input$anchor, min_intron = input$min_intron,
        max_intron = input$max_intron, strand_mode = input$strand_mode,
        exclude_duplicates = input$exclude_dup, nh1_only = input$nh1)
      settings <- backend$validate_config(settings); settings$flank <- flank
      frozen <- controls(); frozen$sampling_mode <- "all"; frozen$sample_cap <- NULL
      list(lines = rows, settings = settings, controls = frozen)
    }
    observeEvent(input$show_jobs, {
      req(is_active())
      bslib::nav_select("result_tabs", "Jobs", session = session)
    })
    observeEvent(input$batch_submit, {
      req(is_active())
      if (batch_task$status() == "running") return()
      batch_review(NULL); batch_review_error(NULL)
      tryCatch({
        snapshot <- batch_inputs()
        batch_review(snapshot)
        maximum_jobs <- length(unique(snapshot$lines$query))
        bslib::nav_select("result_tabs", "Jobs", session = session)
        showModal(modalDialog(title = "Review full-cohort batch resource request", size = "l", easyClose = FALSE,
          p(strong("No jobs have been submitted by this review.")),
          p("Dataset: ", if (is.null(source_identity)) "Current configured dataset" else source_identity$label,
            " · assembly: ", resources$build),
          p(nrow(snapshot$lines), " non-empty input lines; ", maximum_jobs, " unique query strings."),
          div(class = "alert alert-warning",
            strong(paste0("Up to ", maximum_jobs, " independent jobs × 16 cores = ", 16 * maximum_jobs, " requested cores.")),
            p("This is an upper bound, not a validated variant count or a promise of simultaneous allocation. Exact VCF lookup, allele ambiguity, record-level deduplication and source checks run after confirmation; they can reduce the submitted job count.")),
          p("All matched BAMs with complete DNA calls. Analysis flank: ±", snapshot$settings$flank, " bp."),
          tags$details(tags$summary("Review query lines and frozen filters"),
            tags$pre(paste(snapshot$lines$query, collapse = "\n")),
            tags$pre(paste(capture.output(str(snapshot$settings)), collapse = "\n"))),
          p("Each accepted job persists independently of this browser. Previously submitted jobs are not cancelled or rolled back. This review does not detect duplicates in your earlier batches."),
          footer = tagList(actionButton(session$ns("batch_cancel"), "Back without submitting"),
            actionButton(session$ns("batch_confirm"), "Confirm and submit full jobs", class = "btn-primary"))))
      }, error = function(e) batch_review_error(conditionMessage(e)))
      bslib::update_task_button("batch_submit", state = "ready", session = session)
    })
    observeEvent(input$batch_cancel, {
      req(is_active())
      batch_review(NULL); batch_review_error(NULL)
      removeModal()
    })
    observeEvent(is_active(), {
      if (!isTRUE(is_active()) && !is.null(batch_review())) {
        batch_review(NULL)
        removeModal()
      }
    })
    observeEvent(input$batch_confirm, {
      req(is_active())
      # The server, not just the confirmation dialog, prevents unreviewed,
      # stale and repeated submissions. Consume approval before any side effect.
      if (batch_task$status() == "running") return()
      snapshot <- batch_review()
      batch_review(NULL)
      removeModal()
      problem <- tryCatch({
        if (is.null(snapshot)) stop("Review the batch before confirming a submission.")
        if (!identical(snapshot, batch_inputs()))
          stop("Inputs changed after review. Review the batch and resources again; no new jobs were submitted.")
        batch_review_error(NULL)
        batch_snapshot(snapshot)
        bslib::update_task_button("batch_submit", state = "busy", session = session)
        batch_task$invoke(resources, snapshot$lines, snapshot$settings, source_identity, snapshot$controls,
          backend_file, variant_file, job_source, job_root)
        NULL
      }, error = function(e) conditionMessage(e))
      if (!is.null(problem)) {
        batch_review_error(problem)
        bslib::update_task_button("batch_submit", state = "ready", session = session)
      }
    })
    output$batch_review_status <- renderUI({
      if (!is.null(batch_review_error())) return(div(class = "alert alert-warning", batch_review_error()))
      if (!is.null(batch_review())) return(p("Resource review is open. No new jobs have been submitted; confirm or go back."))
      p("Review first, then confirm. Each resolved variant requests its own 16 cores.")
    })
    observeEvent(batch_task$status(), {
      req(is_active())
      if (batch_task$status() == "success") refresh_job_list()
    }, ignoreInit = TRUE)
    batch_outcomes <- reactive({req(is_active()); batch_task$result()})
    output$batch_status <- renderUI({
      s <- batch_task$status()
      if (s == "initial") return(p("Enter multiple coordinates in Batch variants to submit independent full jobs."))
      if (s == "running") return(div(class = "alert alert-info", "Resolving and submitting the frozen batch. Each successful line becomes a separate 16-core job; jobs may start while later lines are being checked."))
      if (s == "error") {
        message <- tryCatch(batch_task$result(), error = function(e) conditionMessage(e))
        return(div(class = "alert alert-warning", "Batch orchestration failed: ", as.character(message), ". Inspect Jobs for any submissions already accepted before retrying."))
      }
      rows <- batch_outcomes(); ok <- sum(rows$outcome == "SUBMITTED"); duplicate <- sum(rows$outcome == "DUPLICATE_IN_BATCH")
      failed <- nrow(rows) - ok - duplicate
      div(class = if (failed) "alert alert-warning" else "alert alert-success",
        ok, " independent jobs submitted; ", duplicate, " duplicate records skipped; ", failed,
        " lines unresolved or failed. Each submitted variant has its own 16-core request and saved result.")
    })
    output$batch_outcomes <- DT::renderDT(batch_outcomes(), rownames = FALSE, selection = "none",
      options = list(pageLength = 10, scrollX = TRUE))
    output$download_batch <- downloadHandler("variant_batch_outcomes.csv", function(file)
      utils::write.csv(batch_outcomes(), file, row.names = FALSE, na = "NA"))
    refresh_job_list <- function() {
      if (is.null(job_api) || !isTRUE(isolate(is_active()))) return(invisible(NULL))
      found <- tryCatch(job_api$list_variant_jobs(job_root), error = function(e) e)
      if (inherits(found, "error")) {
        showNotification(paste("Cannot list saved jobs:", conditionMessage(found)), type = "error")
        return(invisible(NULL))
      }
      job_list(found)
      dashboard_revision(isolate(dashboard_revision()) + 1L)
      choices <- if (nrow(found)) stats::setNames(found$job_dir, found$label) else character()
      updateSelectInput(session, "saved_job", choices = choices,
        selected = if (!is.null(active_job()) && active_job() %in% found$job_dir) active_job() else NULL)
      invisible(found)
    }
    observeEvent(input$refresh_jobs, refresh_job_list())
    observeEvent(TRUE, refresh_job_list(), once = TRUE)
    dashboard_cache <- new.env(parent = emptyenv())
    dashboard_cache$value <- data.frame(job_dir = character(), label = character())
    job_dashboard <- reactivePoll(10000, session,
      checkFunc = function() {
        dashboard_revision()
        if (!isTRUE(is_active()) || is.null(job_api)) {
          dashboard_cache$value <- data.frame(job_dir = character(), label = character())
          return(NULL)
        }
        rows <- tryCatch({
          if (is.function(job_api$list_variant_job_statuses)) job_api$list_variant_job_statuses(job_root) else
            job_api$list_variant_jobs(job_root)
        }, error = function(e) data.frame(job_dir = character(), label = character()))
        dashboard_cache$value <- rows
        jsonlite::toJSON(rows, dataframe = "rows", na = "null")
      }, valueFunc = function() dashboard_cache$value)
    observeEvent(job_dashboard(), {
      req(is_active())
      rows <- job_dashboard()
      job_list(rows[, c("job_dir", "label"), drop = FALSE])
      current <- isolate(input$saved_job)
      selected <- if (length(current) == 1L && current %in% rows$job_dir) current else
        if (!is.null(active_job()) && active_job() %in% rows$job_dir) active_job() else NULL
      updateSelectInput(session, "saved_job", choices = stats::setNames(rows$job_dir, rows$label), selected = selected)
    }, ignoreNULL = TRUE)
    output$jobs_table <- DT::renderDT({
      req(is_active())
      rows <- job_dashboard()
      columns <- setdiff(names(rows), c("job_dir", "label", "source_id"))
      if (!length(columns)) columns <- "label"
      selected <- match(dashboard_selected(), rows$job_dir)
      DT::datatable(rows[, columns, drop = FALSE], rownames = FALSE,
        selection = list(mode = "single", selected = selected[!is.na(selected)]),
        options = list(pageLength = 10, scrollX = TRUE))
    })
    observeEvent(input$jobs_table_rows_selected, {
      req(is_active())
      index <- input$jobs_table_rows_selected; rows <- job_dashboard()
      if (length(index) == 1L && index >= 1L && index <= nrow(rows)) dashboard_selected(rows$job_dir[[index]])
    })
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
      view_revision(view_revision() + 1L)
      selected_logs(NULL)
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
    observeEvent(input$open_dashboard_job, {
      req(is_active())
      tryCatch({
        open_saved_job(dashboard_selected())
        bslib::nav_select("result_tabs", "Genotype comparison", session = session)
      }, error = function(e) showNotification(conditionMessage(e), type = "error", duration = 10))
    })
    resume_pending <- reactiveVal(NULL)
    resume_receipt <- reactiveVal(NULL)
    resume_task <- shiny::ExtendedTask$new(function(path, root, source_file, variant_source, jobs_source, identity, frozen_controls, token) {
      promises::future_promise({
        value <- tryCatch({
          e <- new.env(parent = globalenv())
          sys.source(source_file, envir = e); sys.source(variant_source, envir = e); sys.source(jobs_source, envir = e)
          state <- e$read_variant_job(path, job_root = root)
          if (!is.null(identity) && !identical(state$request$cfg$data_source_id, identity$id))
            stop("The saved job's source identity does not match this dataset.")
          list(ok = TRUE, receipt = e$resume_variant_job(path, job_root = root))
        }, error = function(e) list(ok = FALSE, message = conditionMessage(e)))
        c(value, list(job_dir = path, controls = frozen_controls, view_token = token))
      }, seed = TRUE)
    }) |> bslib::bind_task_button("resume_job")
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
        view_revision(view_revision() + 1L)
        resume_pending(path); resume_receipt(NULL)
        resume_task$invoke(path, job_root, backend_file, variant_file, job_source, source_identity, controls(), view_revision())
      }, error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 10)
        bslib::update_task_button("resume_job", state = "ready", session = session)
      })
    })
    observeEvent(resume_task$status(), {
      req(is_active())
      if (resume_task$status() != "success") return()
      value <- resume_task$result(); resume_receipt(value)
      refresh_job_list()
      if (isTRUE(value$ok) && identical(view_revision(), value$view_token) && identical(controls(), value$controls))
        tryCatch(open_saved_job(value$job_dir), error = function(e) showNotification(conditionMessage(e), type = "error", duration = 10))
    }, ignoreInit = TRUE)
    output$resume_status <- renderUI({
      if (is.null(resume_pending())) return(NULL)
      state <- resume_task$status()
      if (state == "running") return(div(class = "alert alert-info", "Requesting a resume for ", basename(resume_pending()), ". Other queries remain available."))
      if (state == "error") {
        message <- tryCatch(resume_task$result(), error = function(e) conditionMessage(e))
        return(div(class = "alert alert-warning", "Resume failed: ", as.character(message)))
      }
      value <- resume_receipt(); if (is.null(value)) return(NULL)
      div(class = if (isTRUE(value$ok)) "alert alert-success" else "alert alert-warning",
        if (isTRUE(value$ok)) paste0("Resume submitted for ", basename(value$job_dir), "; its independent attempt is listed in Jobs.") else paste0("Resume failed: ", value$message))
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
      status_class <- if (!is.null(state$error) || !is.null(full_result_error())) "danger" else
        if (isTRUE(state$result_ready)) {
          if (!is.null(p$failed) && is.finite(p$failed) && p$failed > 0) "warning" else "success"
        } else if (grepl("FAIL|INTERRUPT|UNAVAILABLE|UNKNOWN", state$status)) "warning" else "info"
      tagList(div(class = paste("alert alert-", status_class, sep = ""),
        strong(paste0("All matched BAMs · ", state$status)),
        if (!is.null(state$job_id)) paste0(" · SCC job ", state$job_id),
        if (!is.null(state$request$cores)) paste0(" · ", state$request$cores, " dedicated cores; ", state$request$workers, " sample workers"),
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
    output$job_files <- renderUI({
      if (!identical(active_mode(), "all") || is.null(active_job())) return(NULL)
      tagList(h5("Selected job's independent files"), p(code(active_job())),
        p("Matching/preparation logs, each submission attempt's stdout/stderr, checkpoints and final results belong to this directory only."),
        actionButton(session$ns("refresh_job_logs"), "Read selected job log tails"),
        verbatimTextOutput(session$ns("job_log_tail")))
    })
    observeEvent(input$refresh_job_logs, {
      req(is_active(), !is.null(active_job()))
      path <- active_job()
      logs <- tryCatch({
        if (!is.function(job_api$read_variant_job_logs)) stop("Bounded log reading is unavailable in this backend; use the selected job directory above.")
        job_api$read_variant_job_logs(path, job_root = job_root, max_lines = 100L)
      }, error = function(e) list(job_dir = path, error = conditionMessage(e)))
      selected_logs(logs)
    })
    output$job_log_tail <- renderPrint({
      req(is_active())
      logs <- selected_logs()
      if (is.null(logs) || !identical(logs$job_dir, active_job())) {
        cat("Use Read selected job log tails to inspect this job's latest attempt.\n"); return(invisible(NULL))
      }
      if (!is.null(logs$error)) {cat(logs$error, "\n"); return(invisible(NULL))}
      cat("Job directory:", logs$job_dir, "\nAttempt:", logs$attempt, "\n")
      for (name in names(logs$logs)) {
        entry <- logs$logs[[name]]
        cat("\n", name, if (isTRUE(entry$truncated)) " (bounded tail)" else "", "\n", sep = "")
        if (!isTRUE(entry$exists)) cat("Not written for this attempt.\n") else cat(paste(entry$lines, collapse = "\n"), "\n")
      }
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
    output$view_context <- renderUI({
      req(is_active())
      mode <- active_mode()
      if (is.null(mode)) return(div(class = "regshiny-context",
        strong("No analysis result selected"),
        p("Find a variant and compare matched RNA samples, or open an existing result in Jobs.")))
      cfg <- launched_config()
      variant <- if (identical(mode, "all")) active_job_variant() else
        tryCatch(result()$variant, error = function(e) NULL)
      label <- if (!is.null(variant)) paste(variant$chrom[[1L]], variant$pos1[[1L]],
        variant$ref[[1L]], variant$alt[[1L]], sep = ":") else
        if (!is.null(launched_controls())) launched_controls()$query else "Submitted analysis"
      state_label <- if (identical(mode, "all")) {
        state <- job_state()
        if (is.null(state)) "Reading job state" else state$status
      } else analysis_task$status()
      div(class = "regshiny-context",
        strong(paste("Viewing:", label)),
        p(if (!is.null(source_identity)) paste0(source_identity$label, " · "),
          if (!is.null(variant)) paste0(variant$build[[1L]], " · "),
          if (identical(mode, "all")) "Full cohort" else "Exploratory preview",
          " · ", state_label,
          if (!is.null(active_job())) paste0(" · ", basename(active_job()))),
        if (!is.null(cfg)) p("Frozen analysis interval: ", cfg$chrom, ":", cfg$start1, "–", cfg$end1),
        if (restored_job()) p("Saved snapshot. The form does not edit this job; plots and downloads belong to this snapshot.") else
          if (results_stale()) p(class = "regshiny-pending", "Unapplied form changes. The viewed analysis retains its submitted variant and parameters."))
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
      single_submit_task = single_submit_task, resume_task = resume_task, batch_task = batch_task, batch_outcomes = batch_outcomes,
      job_dashboard = job_dashboard,
      open_saved_job = open_saved_job, full_result = full_result))
  })
}

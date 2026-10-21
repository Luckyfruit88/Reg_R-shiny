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

# The only full runtime initialization: dependency checks, SCC module paths and
# thread limits. run_app.R and RStudio's Run App button both arrive here.
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
source("single_bam_module.R", local = TRUE)
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

ui <- page_navbar(
  title = "Reg_Shiny", id = "workspace", selected = "Data sources", fillable = FALSE,
  header = tagList(tags$head(tags$link(rel = "stylesheet", href = "frontend.css")),
    uiOutput("active_data_banner", role = "status", `aria-live` = "polite")),
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
    div(class = "regshiny-context",
      strong(paste("Connected dataset:", b$label)),
      if (!inherits(b$variant_resources, "error")) span(" · ", b$variant_resources$build),
      tags$details(tags$summary("Source identity"), code(b$id)))
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

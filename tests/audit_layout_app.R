# Browser-only synthetic fixture. Render the production WGS UI and production
# audit-table expressions. No app.R startup, native tools, cohort data or jobs.
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(bslib))
source("annotation_module.R")
source("splice_evidence_module.R")
source("variant_module.R")

# Reuse the production table renderer, rather than a simplified HTML imitation.
module_call <- Filter(function(x) is.call(x) &&
  identical(x[[1L]], quote(shiny::moduleServer)), as.list(body(variant_server))[-1L])
stopifnot(length(module_call) == 1L)
module_body <- module_call[[1L]][[3L]][[3L]]
renderers <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("for")) &&
  identical(x[[2L]], as.name("name")) &&
  identical(x[[3L]], quote(c("junctions", "samples", "genotypes", "rna_bases"))),
  as.list(module_body)[-1L])
stopifnot(length(renderers) == 1L)
audit_renderer <- renderers[[1L]]

n <- 75L
calls <- data.frame(
  vcf_sample = sprintf("synthetic-%03d", seq_len(n)),
  raw_gt = rep(c("0/0", "0/1", "1/1", "./."), length.out = n),
  genotype = rep(c("0/0", "0/1", "1/1", "NO_CALL"), length.out = n),
  call_status = rep(c("CALLED", "CALLED", "CALLED", "NO_CALL"), length.out = n),
  ploidy = 2L, phased = FALSE, linked = rep(c(TRUE, TRUE, FALSE), length.out = n),
  bam = paste0("/synthetic/not-a-real-dataset/", strrep("long-directory/", 14),
               sprintf("sample-%03d.bam", seq_len(n))),
  stringsAsFactors = FALSE
)
samples <- calls
samples$status <- rep(c("SUCCESS", "FAILED"), length.out = n)
samples$retained_reads <- ifelse(samples$status == "SUCCESS", 100, NA_real_)
samples$error <- ifelse(samples$status == "FAILED", "Synthetic test failure", "")
fixture <- list(samples = samples, genotypes = calls,
                junctions = data.frame(), rna_bases = data.frame())

shiny::addResourcePath("audit-fixture-assets", normalizePath("www"))
ui <- function(request) {
  baseline <- identical(shiny::parseQueryString(request$QUERY_STRING)$baseline, "1")
  bslib::page_navbar(
    title = "Synthetic audit layout check", fillable = FALSE,
    header = if (!baseline) tags$head(tags$link(rel = "stylesheet", href = "audit-fixture-assets/frontend.css")),
    nav_panel("WGS genotype comparison", variant_ui("audit"))
  )
}
server <- function(input, output, session) {
  moduleServer("audit", function(input, output, session) {
    result <- reactive(fixture)
    eval(audit_renderer, envir = environment())
  })
}
message("Audit fixture package versions: ", paste(vapply(c("shiny", "bslib", "DT", "htmlwidgets"),
  function(p) paste0(p, "=", packageVersion(p)), character(1)), collapse = "; "))
shiny::runApp(shinyApp(ui, server), host = "127.0.0.1", port = 8765, launch.browser = FALSE)

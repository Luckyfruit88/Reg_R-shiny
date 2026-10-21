# Native-free structure and dynamic-input checks. No app startup or SCC access.
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(bslib))
source("annotation_module.R")
source("splice_evidence_module.R")
source("variant_module.R")
source("single_bam_module.R")
for (file in c("app.R", "variant_module.R", "single_bam_module.R", "tests/test_variant_ui.R")) parse(file)
ui_text <- as.character(variant_ui("contract"))
for (id in c("query", "batch_queries", "batch_submit", "saved_job", "open_job", "resume_job", "view_context")) {
  needle <- paste0('id="contract-', id, '"')
  matches <- gregexpr(needle, ui_text, fixed = TRUE)[[1L]]
  stopifnot(sum(matches > 0) == 1L)
}
stopifnot(grepl('aria-live="polite"', ui_text, fixed = TRUE))
stopifnot(grepl('Review batch and resources', ui_text, fixed = TRUE))
# Exercise the production display-window wiring in a small module. This isolates
# the UI binding from BAM analysis, whose native test is unchanged.
display_test_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    result <- reactive(list(config = list(start1 = 101, end1 = 450), denominator_reads = 4L))
    view_range <- single_bam_display_server(input, output, session, result)
  })
}
shiny::testServer(display_test_server, {
  rendered <- output$display_control
  html <- if (is.list(rendered)) as.character(rendered$html) else as.character(rendered)
  expected_id <- paste0('id="', session$ns("display"), '"')
  stopifnot(grepl(expected_id, html, fixed = TRUE))
  stopifnot(!grepl('id="display"', html, fixed = TRUE))
  stopifnot(identical(isolate(view_range()), c(101, 450)))
  session$setInputs(display = c(150, 320))
  stopifnot(identical(isolate(view_range()), c(150, 320)))
  stopifnot(isolate(result())$config$start1 == 101, isolate(result())$config$end1 == 450)
  stopifnot(isolate(result())$denominator_reads == 4L)
  session$setInputs(display = c(200, 200))
  stopifnot(identical(isolate(view_range()), c(101, 450)))
  session$setInputs(display = c(NA_real_, 320))
  stopifnot(identical(isolate(view_range()), c(101, 450)))
  session$setInputs(display = c(1, 900))
  stopifnot(identical(isolate(view_range()), c(101, 450)))
})
cat("FRONTEND_CONTRACT_TESTS_OK\n")

# Native-free structure and dynamic-input checks. No app startup or SCC access.
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(bslib))
source("annotation_module.R")
source("splice_evidence_module.R")
source("variant_module.R")
for (file in c("app.R", "variant_module.R", "tests/test_variant_ui.R")) parse(file)
ui_text <- as.character(variant_ui("contract"))
for (id in c("query", "batch_queries", "batch_submit", "saved_job", "open_job", "resume_job", "view_context")) {
  needle <- paste0('id="contract-', id, '"')
  matches <- gregexpr(needle, ui_text, fixed = TRUE)[[1L]]
  stopifnot(sum(matches > 0) == 1L)
}
stopifnot(grepl('aria-live="polite"', ui_text, fixed = TRUE))
stopifnot(grepl('Review batch and resources', ui_text, fixed = TRUE))
# Extract production functions without sourcing startup settings or native tools.
for (expression in parse("app.R")) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.symbol(expression[[2L]]) &&
      as.character(expression[[2L]]) %in% c("single_bam_ui", "single_bam_server")) eval(expression)
}
# Exercise the exact production display-control expressions in a small module.
# This isolates the UI binding from BAM analysis, whose native test is unchanged.
find_assignment <- function(expr, field) {
  if (!is.call(expr)) return(NULL)
  if (identical(expr[[1L]], as.name("<-")) &&
      identical(expr[[2L]], field)) return(expr[[3L]])
  for (part in as.list(expr)[-1L]) {
    found <- find_assignment(part, field)
    if (!is.null(found)) return(found)
  }
  NULL
}
display_expr <- find_assignment(body(single_bam_server), quote(output$display_control))
range_expr <- find_assignment(body(single_bam_server), quote(view_range))
stopifnot(!is.null(display_expr), !is.null(range_expr))
display_test_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    result <- reactive(list(config = list(start1 = 101, end1 = 450), denominator_reads = 4L))
    output$display_control <- eval(display_expr)
    view_range <- eval(range_expr)
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

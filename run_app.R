# In the RStudio Console: source("run_app.R"); run_reg_shiny()
.regshiny_launcher_dir <- local({
  paths <- vapply(sys.frames(), function(frame) {
    value <- frame$ofile
    if (is.character(value) && length(value) == 1L) value else ""
  }, character(1))
  paths <- paths[nzchar(paths)]
  normalizePath(if (length(paths)) dirname(tail(paths, 1L)) else getwd(), mustWork = TRUE)
})

run_reg_shiny <- function(app_dir = .regshiny_launcher_dir, host = "127.0.0.1",
                          port = getOption("shiny.port"),
                          launch.browser = getOption("shiny.launch.browser", interactive())) {
  app_dir <- normalizePath(path.expand(app_dir), mustWork = TRUE)
  if (!file.exists(file.path(app_dir, "app.R"))) stop("app_dir must be the Reg_R-shiny checkout.")
  runtime <- new.env(parent = globalenv())
  sys.source(file.path(app_dir, "runtime_setup.R"), envir = runtime)
  runtime$regshiny_runtime_setup()
  # RStudio provides its own viewer/proxy through the standard Shiny launch option.
  shiny::runApp(appDir = app_dir, host = host, port = port, launch.browser = launch.browser)
}

# Optional only: SCC R 4.5.2 already provides these shared packages.
# Run from the clone root. The package list lives only in runtime_setup.R.
deps <- new.env()
sys.source("src/runtime_setup.R", envir = deps)
packages <- deps$regshiny_required_packages()
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (requireNamespace("shiny", quietly = TRUE) && utils::packageVersion("shiny") < deps$regshiny_min_shiny_version)
  missing <- union(missing, "shiny")
if (!length(missing)) {
  message("Reg_Shiny R dependencies are already available.")
} else {
  short_version <- paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1L]][[1L]], sep = ".")
  user_library <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(user_library)) user_library <- file.path("~", "R", paste0(R.version$platform, "-library"), short_version)
  user_library <- strsplit(user_library, .Platform$path.sep, fixed = TRUE)[[1L]][[1L]]
  user_library <- gsub("%p", R.version$platform, user_library, fixed = TRUE)
  user_library <- gsub("%v", short_version, user_library, fixed = TRUE)
  user_library <- gsub("%V", paste(R.version$major, R.version$minor, sep = "."), user_library, fixed = TRUE)
  user_library <- path.expand(user_library)
  dir.create(user_library, recursive = TRUE, showWarnings = FALSE)
  if (file.access(user_library, 2) != 0L) stop("Choose a writable personal R_LIBS_USER before installing.")
  .libPaths(c(user_library, .libPaths()))
  install.packages(missing, lib = user_library, repos = "https://cloud.r-project.org")
}

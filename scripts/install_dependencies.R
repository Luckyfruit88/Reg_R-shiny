# Optional only: SCC R 4.5.2 already provides these shared packages.
packages <- c("shiny", "bslib", "DT", "data.table", "processx", "future", "promises", "jsonlite")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (requireNamespace("shiny", quietly = TRUE) && utils::packageVersion("shiny") < "1.8.1")
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

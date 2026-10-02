# Shared runtime bootstrap for RStudio's Run App button and run_app.R.
# Native module imports affect this R process and its children, not other sessions.
regshiny_runtime_setup <- function(strict_native = FALSE, import_scc_modules = TRUE) {
  required <- c("shiny", "bslib", "DT", "data.table", "processx", "future", "promises", "jsonlite")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "),
    ". In SCC RStudio select R 4.5.2, or run scripts/install_dependencies.R in your own library.")
  if (utils::packageVersion("shiny") < "1.8.1") stop("Reg_Shiny requires shiny >= 1.8.1.")
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  native <- c("samtools", "regtools", "bcftools", "python3")
  scheduler <- c("qsub", "qstat", "qacct")
  absent <- function(names) names[!nzchar(Sys.which(names))]
  is_scc <- dir.exists("/share/module.8") && file.exists("/bin/bash")
  imported <- FALSE; module_error <- NULL
  if (import_scc_modules && is_scc && length(absent(c(native, scheduler)))) {
    # Do not print or import the complete login environment: only executable,
    # native-library and Grid Engine locations are needed by the application.
    keys <- c("PATH", "LD_LIBRARY_PATH", "SGE_ROOT", "SGE_CELL", "SGE_QMASTER_PORT", "SGE_EXECD_PORT")
    py <- "import json,os; print(json.dumps({k:os.environ.get(k,'') for k in ['PATH','LD_LIBRARY_PATH','SGE_ROOT','SGE_CELL','SGE_QMASTER_PORT','SGE_EXECD_PORT']}))"
    shell <- paste("module load samtools/1.23 regtools/1.0.0 bcftools/1.23 && python3 -c", shQuote(py))
    result <- tryCatch(processx::run("/bin/bash", c("-lc", shell), timeout = 60,
      error_on_status = FALSE, cleanup_tree = FALSE), error = function(e) e)
    if (inherits(result, "error")) module_error <- conditionMessage(result) else if (result$status != 0L) {
      module_error <- trimws(result$stderr)
    } else {
      lines <- strsplit(result$stdout, "\n", fixed = TRUE)[[1L]]
      json <- lines[startsWith(lines, "{") & endsWith(lines, "}")]
      imported_values <- if (length(json)) tryCatch(jsonlite::fromJSON(tail(json, 1L)), error = function(e) NULL) else NULL
      if (is.list(imported_values) && all(keys %in% names(imported_values))) {
        values <- imported_values[keys]
        # Empty optional SGE values must not erase a working caller setting.
        values <- values[vapply(values, function(x) is.character(x) && length(x) == 1L && nzchar(x), logical(1))]
        if (length(values)) do.call(Sys.setenv, values)
        imported <- TRUE
      } else module_error <- "The SCC login shell did not return a usable module environment."
    }
  }
  missing_native <- absent(native)
  if (length(missing_native)) {
    message <- paste0("Native tools unavailable: ", paste(missing_native, collapse = ", "),
      ". Configure the SCC RStudio session with samtools/1.23 regtools/1.0.0 bcftools/1.23.")
    if (strict_native) stop(message, if (!is.null(module_error)) paste0(" ", module_error))
    warning(message, call. = FALSE)
  }
  tmp <- tempfile("regshiny-runtime-")
  if (!dir.create(tmp, mode = "0700")) stop("A private temporary directory cannot be created.")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  if (.Platform$OS.type == "unix" && bitwAnd(as.integer(file.info(tmp)$mode), 63L) != 0L)
    stop("The temporary directory is not private to this user.")
  invisible(list(scc = is_scc, modules_imported = imported, native = Sys.which(native),
    scheduler = Sys.which(scheduler), missing_native = missing_native,
    missing_scheduler = absent(scheduler), R_version = R.version.string,
    package_versions = vapply(required, function(x) as.character(utils::packageVersion(x)), character(1)),
    private_tempdir = TRUE))
}

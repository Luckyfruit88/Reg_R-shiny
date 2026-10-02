args <- commandArgs(trailingOnly = TRUE)
app <- normalizePath(args[[1L]], mustWork = TRUE)
expressions <- parse(app)
stopifnot(identical(expressions[[1L]][[1L]], as.name("local")),
          identical(expressions[[2L]][[2L]], as.name("needed")))
bootstrap <- expressions[[1L]]
fixture <- tempfile("regshiny_autoload_"); dir.create(fixture)
oldwd <- getwd(); setwd(fixture)
cleanup <- function() {setwd(oldwd); unlink(fixture, recursive = TRUE)}
tryCatch({
  Sys.unsetenv(c("R_ENVIRON_USER", "REGSHINY_VCF_MANIFEST", "REGSHINY_SAMPLE_MANIFEST", "REGSHINY_FULL_CORES"))
  # No config: demo and externally preconfigured resources must remain usable.
  Sys.setenv(REGSHINY_FULL_CORES = "4")
  eval(bootstrap)
  stopifnot(Sys.getenv("REGSHINY_FULL_CORES") == "4")
  dir.create("config")
  writeLines(c('REGSHINY_VCF_MANIFEST="/synthetic path/vcf.tsv"',
               'REGSHINY_SAMPLE_MANIFEST=/synthetic/samples.tsv',
               'REGSHINY_FULL_CORES=8'), "config/runtime.Renviron")
  eval(bootstrap)
  stopifnot(Sys.getenv("REGSHINY_VCF_MANIFEST") == "/synthetic path/vcf.tsv",
            Sys.getenv("REGSHINY_SAMPLE_MANIFEST") == "/synthetic/samples.tsv",
            Sys.getenv("REGSHINY_FULL_CORES") == "8")
  # Existing launchers retain explicit override priority (including relative paths).
  writeLines('REGSHINY_FULL_CORES=6', "custom.Renviron")
  Sys.setenv(R_ENVIRON_USER = "custom.Renviron")
  eval(bootstrap)
  stopifnot(Sys.getenv("REGSHINY_FULL_CORES") == "6")
  Sys.setenv(R_ENVIRON_USER = normalizePath("custom.Renviron"))
  eval(bootstrap)
  stopifnot(Sys.getenv("REGSHINY_FULL_CORES") == "6")
  # An explicit mistake must not quietly fall back to the project config.
  Sys.setenv(R_ENVIRON_USER = "missing.Renviron")
  err <- tryCatch({eval(bootstrap); NULL}, error = conditionMessage)
  stopifnot(!is.null(err), grepl("does not exist", err), Sys.getenv("REGSHINY_FULL_CORES") == "6")
  Sys.setenv(R_ENVIRON_USER = "config")
  err <- tryCatch({eval(bootstrap); NULL}, error = conditionMessage)
  stopifnot(!is.null(err), grepl("regular file", err))
  writeLines('REGSHINY_FULL_CORES no_equals', "malformed.Renviron")
  Sys.setenv(R_ENVIRON_USER = "malformed.Renviron")
  err <- tryCatch({eval(bootstrap); NULL}, error = conditionMessage)
  stopifnot(!is.null(err), grepl("Invalid Reg_Shiny runtime configuration syntax", err))
  # Child R processes inherit loaded app settings without R_ENVIRON_USER.
  Sys.unsetenv("R_ENVIRON_USER"); eval(bootstrap)
  writeLines('stopifnot(Sys.getenv("REGSHINY_FULL_CORES") == "8"); cat("CHILD_ENV_PASS\\n")', "child.R")
  output <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "child.R"), stdout = TRUE, stderr = TRUE)
  stopifnot(is.null(attr(output, "status")), any(grepl("CHILD_ENV_PASS", output)))
  cat("PASS: parse/order, absent default, default load, quoting, explicit relative/absolute override, invalid/malformed/directory override, child inheritance\n")
}, finally = cleanup())

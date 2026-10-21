# Worktree guard of tests/run_tests.R, in a throwaway git repository.
# Status codes alone cannot show that a dirty file was rewritten or removed.
source("tests/run_tests.R")
git <- function(dir, ...) stopifnot(identical(system2("git", c("-C", shQuote(dir), ...),
  stdout = FALSE, stderr = FALSE), 0L))
repo <- tempfile("regshiny_worktree_"); dir.create(repo)
git(repo, "init", "-q")
writeLines("x <- 1", file.path(repo, "tracked.R"))
git(repo, "add", "tracked.R")
git(repo, "-c", "user.name=test", "-c", "user.email=test@example.invalid", "commit", "-q", "-m", "init")
verdict <- function(change, require_clean = FALSE) {
  before <- worktree_snapshot(repo); change(); worktree_verdict(before, worktree_snapshot(repo), require_clean)
}
fails_on <- function(v, path) stopifnot(!v$ok, path %in% v$changes)

v <- verdict(function() NULL)
stopifnot(v$ok, v$message == "worktree clean before and after")
stopifnot(verdict(function() NULL, require_clean = TRUE)$ok)
fails_on(verdict(function() writeLines("new", file.path(repo, "new-output.txt"))), "new-output.txt")

# Dirty before the run: allowed, but the tests must leave it exactly as it was.
writeLines("x <- 2", file.path(repo, "tracked.R"))
writeLines("keep", file.path(repo, "scratch.txt"))
v <- verdict(function() NULL)
stopifnot(v$ok, grepl("did not alter", v$message))
fails_on(verdict(function() writeLines("x <- 3", file.path(repo, "tracked.R"))), "tracked.R")
fails_on(verdict(function() unlink(file.path(repo, "scratch.txt"))), "scratch.txt")
stopifnot(!verdict(function() NULL, require_clean = TRUE)$ok)

# A failed git query is unverified, never a pass.
snapshot <- worktree_snapshot(repo)
stopifnot(!worktree_verdict(snapshot, NULL)$ok, grepl("UNVERIFIED", worktree_verdict(snapshot, NULL)$message))
stopifnot(!worktree_verdict(NULL, NULL, require_clean = TRUE)$ok)
nongit <- tempfile("regshiny_nongit_"); dir.create(nongit)
Sys.setenv(GIT_CEILING_DIRECTORIES = dirname(nongit))
stopifnot(is.null(worktree_snapshot(nongit)))
unlink(c(repo, nongit), recursive = TRUE)
cat("Worktree guard tests: PASS\n")

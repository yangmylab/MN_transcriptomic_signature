# Paper 1 - shared library bootstrap. Source this at the top of every R script.
# The system/default R library is read-only, so we prepend a workspace-local library
# while KEEPING the default libs (which carry BiocManager, remotes, GEOquery, Seurat...).
#
# The project root is resolved from THIS file's own location (falling back to the running
# script, then to the working directory), so the project can be moved between machines or
# drives without editing every script. Set PAPER1_PROJ to override.

.paper1_self <- function() {
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(f) || !nzchar(f)) {
    a <- commandArgs(FALSE)
    f <- sub("^--file=", "", a[startsWith(a, "--file=")])
    if (!length(f)) f <- getwd()
  }
  normalizePath(f[1], mustWork = FALSE)
}

PROJ <- Sys.getenv("PAPER1_PROJ", unset = "")
if (!nzchar(PROJ)) PROJ <- dirname(dirname(.paper1_self()))
PROJ <- normalizePath(PROJ, mustWork = FALSE)

WS_LIB <- file.path(PROJ, ".Rlib")
if (!dir.exists(WS_LIB)) dir.create(WS_LIB, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(WS_LIB, .libPaths()))
options(Ncpus = max(1, parallel::detectCores() - 1), timeout = 1800)

# Data root: the GEO inputs live in the SHARED MN dataset collection (../Datasets), not in the
# project folder, because Paper_2 (single-cell/spatial) and Paper_3 (MR drug targets) reuse several
# of the same accessions. Moved there 2026-10-19; each accession is stored once as
#   <DATA>/<GSE>/          extracted (raw/, processed/, series matrices)
#   <DATA>/<GSE>_RAW.tar   the original archive
# Override with PAPER1_DATA; fall back to a project-local Datasets/ if the shared one is absent.
DATA <- Sys.getenv("PAPER1_DATA", unset = "")

if (!nzchar(DATA)) {
  shared_data <- file.path(dirname(PROJ), "Datasets")
  DATA <- if (dir.exists(shared_data)) shared_data else file.path(PROJ, "Datasets")
}

DATA <- normalizePath(DATA, mustWork = FALSE)

OUT  <- file.path(PROJ, "output")

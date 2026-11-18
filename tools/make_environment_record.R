#!/usr/bin/env Rscript
# Paper 1 - reproducibility record.
#
# Writes an environment snapshot that can rebuild the analysis on another machine:
#   environment/renv.lock                    renv-compatible package lockfile (CRAN/Bioconductor)
#   environment/sessionInfo.txt              R version, platform, BLAS, loaded namespaces
#   environment/packages_<lib>.csv           every installed package + version per library
#   environment/README.md                    reconstruction recipe
#
# The project deliberately does NOT activate renv at run time (the analysis resolves its library in
# R/_lib.R); the lockfile is written as a record and for `renv::restore()` on a fresh machine.

# tools/ sits one level below the project root, so the library bootstrap is ../R/_lib.R
local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[startsWith(a, "--file=")])
  here <- if (length(f)) dirname(normalizePath(f[1])) else getwd()
  source(file.path(dirname(here), "R", "_lib.R"))
})

env_dir <- file.path(PROJ, "environment")
dir.create(env_dir, showWarnings = FALSE, recursive = TRUE)

sink(file.path(env_dir, "sessionInfo.txt"))
cat("Project:", PROJ, "\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")
print(sessionInfo())
cat("\n--- library paths ---\n"); print(.libPaths())
cat("\n--- BLAS ---\n"); print(extSoftVersion()["BLAS"])
sink()

libs <- .libPaths()
libs <- libs[dir.exists(libs)]
for (l in libs) {
  ip <- installed.packages(lib.loc = l)
  if (!nrow(ip)) next
  tag <- switch(basename(l),
                ".Rlib" = "workspace",
                "4.6" = "user",
                gsub("[^A-Za-z0-9]+", "_", basename(l)))
  f <- file.path(env_dir, paste0("packages_", tag, ".csv"))
  write.csv(ip[, c("Package", "Version", "LibPath")], f, row.names = FALSE)
  cat("wrote", f, "(", nrow(ip), "packages )\n")
}

# renv lockfile: the workspace library holds locally built Brainarray CDF packages that cannot be
# restored from a repository, so they are marked Local and are rebuilt by R/00d_install_cdf.R and
# R/02a_build_cdf.R.
build_lockfile <- function() {
  ws <- installed.packages(lib.loc = WS_LIB)
  us <- installed.packages(lib.loc = .libPaths()[.libPaths() != WS_LIB])
  us <- us[!rownames(us) %in% rownames(ws), , drop = FALSE]
  mk <- function(ip, source) {
    lapply(seq_len(nrow(ip)), function(i) {
      list(Package = ip[i, "Package"], Version = ip[i, "Version"],
           Source = source,
           Repository = if (identical(source, "Bioconductor")) "Bioconductor" else "CRAN")
    })
  }
  entries <- c(mk(us, "Repository"), mk(ws, "Repository"))
  names(entries) <- vapply(entries, function(e) e$Package, character(1))
  local_cdf <- c("hugene21sthsentrezgcdf", "hgu133plus2hsentrezgcdf",
                 "hgu133plus2hsentrezg2cdf", "hgu133ahsentrezgcdf", "hgu133ahsentrezg2cdf",
                 "hgu133ahsentrezgcdfb", "hgu133plus2hsentrezgcdfb")
  for (nm in intersect(local_cdf, names(entries))) {
    entries[[nm]]$Source <- "Local"
    entries[[nm]]$Repository <- NULL
    entries[[nm]]$Path <- file.path("work", "cdf", nm)
  }
  list(R = list(Version = paste(R.version$major, R.version$minor, sep = "."),
                Repositories = list(CRAN = list(URL = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"),
                                    BioCsoft = list(URL = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor"))),
       Packages = entries)
}

lock <- build_lockfile()
lockfile <- file.path(env_dir, "renv.lock")
if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(lock, lockfile, pretty = TRUE, auto_unbox = TRUE, null = "null")
} else {
  cat("jsonlite unavailable; writing a package CSV only\n")
}
cat("wrote", lockfile, "\n")

readme <- c(
  "# Environment record",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("R: ", R.version.string, " on ", R.version$platform),
  "",
  "## Rebuild on another machine",
  "",
  "1. Install R 4.6.x and Bioconductor 3.23.",
  "2. Restore the repository packages from the lockfile:",
  "",
  "   ```r",
  "   install.packages(\"renv\")",
  "   renv::restore(lockfile = \"environment/renv.lock\")",
  "   ```",
  "",
  "3. The locally built Brainarray CDF packages are NOT in any repository. Rebuild them from the",
  "   `.cdf.gz` files inside the GEO archives (no network needed):",
  "",
  "   ```",
  "   Rscript R/00d_install_cdf.R      # makecdfenv into ./.Rlib",
  "   Rscript R/02a_build_cdf.R        # build the 5 CDF source packages from $DATA/*/raw/",
  "   for p in work/cdf/*; do R CMD INSTALL -l .Rlib \"$p\"; done",
  "   ```",
  "",
  "4. Point the data root at the shared dataset collection (or copy it next to the project):",
  "",
  "   ```",
  "   export PAPER1_DATA=/path/to/MN/Datasets",
  "   ```",
  "",
  "5. Verify: all six verifiers must exit 0.",
  "",
  "   ```",
  "   for v in phase_a phase_b downstream scrna phase_d phase_e; do Rscript R/99_verify_$v.R; done",
  "   ```",
  "",
  "## Notes",
  "",
  "- The analysis library is `./.Rlib` (prepended by `R/_lib.R`); the user library is",
  "  `~/R/x86_64-pc-linux-gnu-library/4.6`. Both are recorded in `packages_*.csv`.",
  "- `renv` is intentionally NOT activated at run time; `renv.lock` is a reproducibility record.",
  "- Packages with `Source = Local` are the Brainarray CDF builds listed above.")
writeLines(readme, file.path(env_dir, "README.md"))
cat("wrote", file.path(env_dir, "README.md"), "\n")

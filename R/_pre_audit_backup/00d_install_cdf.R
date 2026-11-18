#!/usr/bin/env Rscript
# Paper 1 - Phase A step 0d: install the one package Phase A needs that is not already present.
# The Brainarray .cdf.gz files are shipped inside the GEO RAW tars, but converting them into an
# installable CDF package needs makecdfenv. Everything else (affy, affyio, Biobase, preprocessCore,
# oligo, limma, edgeR, sva, GSVA, ...) is already installed in the workspace .Rlib or in the
# user library. Packages go into the workspace library, never into the user's personal one.

WS_LIB <- file.path(dirname(dirname(normalizePath(sub("^--file=", "",
  commandArgs(FALSE)[startsWith(commandArgs(FALSE), "--file=")][1]))), ".Rlib")

.libPaths(c(WS_LIB, .libPaths()))

options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

cat("library paths:\n"); print(.libPaths())

need <- c("makecdfenv")

for (p in need) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat("present :", p, "\n")
  } else {
    cat("installing:", p, "\n")
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = WS_LIB)
    BiocManager::install(p, lib = WS_LIB, ask = FALSE, update = FALSE)
  }
}

cat("\n--- result ---\n")

for (p in need) cat(sprintf("%-14s %s\n", p, if (requireNamespace(p, quietly = TRUE)) "OK" else "FAILED"))

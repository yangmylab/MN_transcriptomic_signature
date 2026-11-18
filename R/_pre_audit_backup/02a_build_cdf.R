#!/usr/bin/env Rscript
# Paper 1 - Phase A step 2a: turn the Brainarray .cdf.gz files (shipped inside the GEO RAW tars)
# into installable CDF packages, so every array cohort can be RMA-normalized against the same
# custom CDF the original authors used. Packages are installed into the workspace .Rlib only.

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages(library(makecdfenv))

cdf_specs <- data.frame(
  gz = c(
    file.path(DATA, "GSE200828", "raw", "GPL19983_hugene21st_Hs_ENTREZG.cdf.gz"),
    file.path(DATA, "GSE104948", "raw", "GPL22945_HGU133Plus2_Hs_ENTREZG.cdf.gz"),
    file.path(DATA, "GSE99339",  "raw", "GPL19109_HGU133Plus2_Hs_ENTREZG.cdf.gz"),
    file.path(DATA, "GSE104948", "raw", "GPL24120_HGU133A_Hs_ENTREZG.cdf.gz"),
    file.path(DATA, "GSE99339",  "raw", "GPL19184_HGU133A_Hs_ENTREZG.cdf.gz")
  ),
  # affy's cleancdfname() normalises a CDF package name by appending "cdf" unless the name already
  # ends in it, so both variants must be named accordingly (the ERCB re-release and the NEPTUNE
  # release ship different Brainarray revisions of the same array)
  pkg = c("hugene21sthsentrezgcdf", "hgu133plus2hsentrezgcdf", "hgu133plus2hsentrezg2cdf",
          "hgu133ahsentrezgcdf", "hgu133ahsentrezg2cdf"),
  array = c("HuGene-2_1-st", "HG-U133-Plus2", "HG-U133-Plus2", "HG-U133A", "HG-U133A"),
  stringsAsFactors = FALSE
)

work <- file.path(PROJ, "work", "cdf")

dir.create(work, showWarnings = FALSE, recursive = TRUE)

old_wd <- setwd(work)

on.exit(setwd(old_wd))

# binary-safe gzip decompression (the CDF is an XDA binary file, not text)

gunzip_bin <- function(src, dst) {
  ci <- gzfile(src, "rb")
  co <- file(dst, "wb")
  repeat {
    chunk <- readBin(ci, "raw", n = 1048576)
    if (length(chunk) == 0) break
    writeBin(chunk, co)
  }
  close(ci)
  close(co)
}

walk <- function(i) {
  pkg <- cdf_specs$pkg[i]
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("[", i, "]", pkg, "- already installed\n")
    return(invisible(NULL))
  }
  # make.cdf.package resolves its input as file.path(getwd(), filename), so pass a bare name
  cdf <- paste0(pkg, ".cdf")
  cat("[", i, "]", pkg, "- decompressing", basename(cdf_specs$gz[i]), "\n")
  gunzip_bin(cdf_specs$gz[i], cdf)
  cat("[", i, "]", pkg, "- building package from", file.size(cdf), "byte CDF\n")
  make.cdf.package(cdf, packagename = pkg, species = "Homo_sapiens", compress = FALSE)
  ok <- requireNamespace(pkg, quietly = TRUE)
  cat("[", i, "]", pkg, if (ok) "OK" else "FAILED", "\n")
  invisible(NULL)
}

cat(".libPaths():\n"); print(.libPaths())

for (i in seq_len(nrow(cdf_specs))) walk(i)

cat("\n------------ CDF PACKAGE STATUS ------------\n")

for (p in cdf_specs$pkg) {
  cat(sprintf("%-26s %s\n", p, if (requireNamespace(p, quietly = TRUE)) "OK" else "FAILED"))
}

writeLines(c(
  paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 02a_build_cdf.R"),
  paste0("  ", cdf_specs$pkg, " (", cdf_specs$array, ") ",
         ifelse(sapply(cdf_specs$pkg, function(p) requireNamespace(p, quietly = TRUE)), "OK", "FAILED"))
), file.path(PROJ, "logs", "build_cdf.log"))

cat("\nDone.\n")

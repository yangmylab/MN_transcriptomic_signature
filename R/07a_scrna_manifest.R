#!/usr/bin/env Rscript
# Paper 1 - scRNA arm, step 07a: validate the user-supplied single-cell archives, extract them,
# and inventory samples. This script NEVER downloads anything.
#
# GEO RAW.tar files come in two layouts, both handled here:
#   1. 10x triplets, flat in one folder:  <GSM>_<sample>_{barcodes,features,matrix}.*
#      -> organised into extracted/<GSE>/samples/<sample_id>/{barcodes,features,matrix}.tsv/mtx.gz
#   2. per-sample dense matrices:         <GSM>_<title>.expression_matrix.txt.gz
#      -> used in place (07b parses orientation)
#
# Inputs : scrna/manifest.tsv, scrna/raw/<GSE>/*.tar, scrna/sample_metadata.tsv (07a0)
# Outputs: scrna/extracted/<GSE>/..., scrna/readiness.csv, scrna/samples.csv

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

scrna_dir <- file.path(PROJ, "scrna")
ext_dir <- file.path(scrna_dir, "extracted")
meta_file <- file.path(scrna_dir, "sample_metadata.tsv")

dir.create(ext_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(file.path(scrna_dir, "manifest.tsv"))) stop("missing scrna/manifest.tsv")

man <- read.delim(file.path(scrna_dir, "manifest.tsv"), stringsAsFactors = FALSE)

cat("scRNA manifest:", nrow(man), "datasets\n")

# --- helpers -----------------------------------------------------------------

# 10x triplets arrive flat in one folder: <prefix>_{barcodes,features,matrix}.*
# We deliberately do NOT build one folder of links per sample any more: this volume is **exfat**,
# which supports neither symlinks nor hardlinks, so a link layout silently produces empty
# directories (R's file.symlink() just returns FALSE) and a copy layout doubles the disk use -
# that is exactly what happened to the 3.5 GB of dereferenced trios found on 2026-10-19.
# Instead we record the shared directory + the filename prefix, and 07b reads the triple directly
# with Seurat::ReadMtx().
inventory_10x <- function(dest) {

  files <- list.files(dest, recursive = TRUE, full.names = TRUE)

  trio <- grepl("_(barcodes|features|genes|matrix)\\.(tsv|mtx)\\.gz$", basename(files))

  if (!any(trio)) return(NULL)

  prefix <- sub("_(barcodes|features|genes|matrix)\\.(tsv|mtx)\\.gz$", "", basename(files[trio]))

  dirs <- tapply(dirname(files[trio]), prefix, function(d) unique(d)[1])

  data.frame(prefix = names(dirs), dir = as.character(dirs), stringsAsFactors = FALSE)
}

# --- inspect each dataset -----------------------------------------------------

rows <- list()
sample_rows <- list()

for (i in seq_len(nrow(man))) {

  acc <- man$accession[i]
  arch <- file.path(PROJ, man$archive[i])

  status <- "MISSING"
  note <- ""
  found <- character(0)
  fmt <- if ("format" %in% names(man)) man$format[i] else "10x"

  if (!file.exists(arch)) {

    note <- "archive not found - see scrna/DATA_REQUEST.md"

  } else if (grepl("\\.rds\\.gz$", arch)) {

    status <- "READY"
    found <- arch
    fmt <- "rds"
    note <- "pre-annotated Seurat object"

  } else {

    dest <- file.path(ext_dir, acc)

    if (!dir.exists(dest) || !length(list.files(dest))) {
      dir.create(dest, showWarnings = FALSE, recursive = TRUE)
      cat("  extracting", basename(arch), "->", dest, "\n")
      utils::untar(arch, exdir = dest)
    }

    files <- list.files(dest, recursive = TRUE, full.names = TRUE)

    if (any(grepl("\\.expression_matrix\\.txt\\.gz$", basename(files)))) {
      found <- files[grepl("\\.expression_matrix\\.txt\\.gz$", basename(files))]
      fmt <- "expression_matrix"
      inv <- NULL
    } else {
      inv <- inventory_10x(dest)
      found <- if (is.null(inv)) character(0) else inv$prefix
      fmt <- "10x"
    }

    status <- if (length(found)) "READY" else "NO_10X_FOUND"

    note <- paste(length(found), "samples of format", fmt)
  }

  rows[[acc]] <- data.frame(
    accession = acc, dataset = man$dataset[i], layer = man$layer[i], priority = man$priority[i],
    expected_samples = man$expected_samples[i], detected_samples = length(found),
    format = fmt, status = status, archive = man$archive[i], note = note,
    stringsAsFactors = FALSE)

  if (length(found)) {
    if (identical(fmt, "expression_matrix")) {
      sample_rows[[acc]] <- data.frame(
        accession = acc, dataset = man$dataset[i], layer = man$layer[i], format = fmt,
        # one dense matrix per sample, used in place
        path = sub(paste0("^", PROJ, "/?"), "", found),
        prefix = NA_character_,
        sample_id = sub("\\.expression_matrix\\.txt\\.gz$", "", basename(found)),
        stringsAsFactors = FALSE)
    } else {
      sample_rows[[acc]] <- data.frame(
        accession = acc, dataset = man$dataset[i], layer = man$layer[i], format = fmt,
        # all samples of a 10x dataset share one untarred directory; prefix selects the triple
        path = sub(paste0("^", PROJ, "/?"), "", inv$dir),
        prefix = inv$prefix,
        sample_id = inv$prefix, stringsAsFactors = FALSE)
    }
  }
}

readiness <- do.call(rbind, rows)

write.csv(readiness, file.path(scrna_dir, "readiness.csv"), row.names = FALSE)

if (length(sample_rows)) {

  samples <- do.call(rbind, sample_rows)

  if (file.exists(meta_file)) {

    meta <- read.delim(meta_file, stringsAsFactors = FALSE)

    samples <- merge(samples, meta[, c("sample_id", "gsm", "title", "patient_id",
                                       "disease_group", "source_accession")],
                     by = "sample_id", all.x = TRUE, sort = FALSE)
  } else {
    samples$gsm <- NA_character_
    samples$title <- NA_character_
    samples$patient_id <- NA_character_
    samples$disease_group <- NA_character_
    samples$source_accession <- samples$accession
  }

  write.csv(samples, file.path(scrna_dir, "samples.csv"), row.names = FALSE)

  cat("\nsamples inventoried:", nrow(samples),
      "| with a disease label:", sum(!is.na(samples$disease_group)), "\n")

  print(table(samples$accession, samples$disease_group, useNA = "ifany"))
} else {
  cat("\nno samples detected yet - nothing written to scrna/samples.csv\n")
}

cat("\n================ MANIFEST READINESS ================\n")
print(readiness[, c("accession", "layer", "format", "priority", "expected_samples",
                    "detected_samples", "status")], row.names = FALSE)

tier1 <- readiness[readiness$priority == "tier1", ]

cat("\ntier-1 datasets ready:", sum(tier1$status == "READY"), "of", nrow(tier1), "\n")

if (all(tier1$status == "READY")) {
  cat("scRNA arm: READY for 07b\n")
} else {
  cat("scRNA arm: NOT READY\n")
}

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 07a_scrna_manifest.R"),
             capture.output(print(readiness, row.names = FALSE))),
           file.path(PROJ, "logs", "scrna_manifest.log"))

cat("\nDone.\n")

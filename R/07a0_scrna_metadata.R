#!/usr/bin/env Rscript
# Paper 1 - scRNA arm, step 07a0: build scrna/sample_metadata.tsv from the GEO series matrices.
#
# The RAW.tar archives carry no phenotype table, and the Series Matrix files are metadata only
# (sample titles + characteristics), not expression data. They were downloaded into
# scrna/metadata/ and are parsed here so every 10x sample gets an explicit disease group.
#
# Output: scrna/sample_metadata.tsv
#   sample_id      prefix shared by the archive files for that sample
#   gsm, title, patient_id, disease_group, source_accession, raw_characteristic

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

scrna_dir <- file.path(PROJ, "scrna")
meta_dir <- file.path(scrna_dir, "metadata")

accs <- c("GSE271188", "GSE171458", "GSE241302")

canonical <- function(x) {
  x <- tolower(x)
  out <- rep("Other", length(x))
  out[grepl("membranous|\\bimn\\b|\\bmn\\b", x)] <- "MN"
  out[grepl("healthy|control|normal", x)] <- "Control"
  out[grepl("iga", x)] <- "IgAN"
  out[grepl("minimal change|\\bmcd\\b", x)] <- "MCD"
  out[grepl("diabetic|\\bdkd\\b", x)] <- "DKD"
  out[grepl("lupus|\\bsle\\b", x)] <- "LN"
  out
}

parse_sm <- function(f) {

  lines <- readLines(gzfile(f), warn = FALSE)

  take <- function(tag) {
    l <- lines[startsWith(lines, tag)]
    if (!length(l)) return(NULL)
    lapply(l, function(x) gsub('^"|"$', "", strsplit(x, "\t", fixed = TRUE)[[1]][-1]))
  }

  gsm <- take("!Sample_geo_accession")[[1]]
  title <- take("!Sample_title")[[1]]
  chars <- take("!Sample_characteristics_ch1")

  n <- length(gsm)

  raw <- rep(NA_character_, n)

  for (j in seq_len(n)) {
    vals <- vapply(chars, function(v) if (length(v) >= j) v[j] else "", "")
    vals <- vals[nzchar(vals)]
    hit <- vals[grepl("diagnosis|disease state|subject status|group", vals, ignore.case = TRUE)]
    if (length(hit)) raw[j] <- hit[1]
  }

  data.frame(gsm = gsm, title = title, raw_characteristic = raw, stringsAsFactors = FALSE)
}

rows <- list()

for (acc in accs) {

  f <- file.path(meta_dir, paste0(acc, "_series_matrix.txt.gz"))

  if (!file.exists(f)) {
    warning("missing series matrix for ", acc)
    next
  }

  sm <- parse_sm(f)

  # sample_id = the archive file prefix (everything before the 10x suffix)
  arch <- file.path(scrna_dir, "raw", acc, paste0(acc, "_RAW.tar"))
  members <- if (file.exists(arch)) untar(arch, list = TRUE) else character(0)

  prefix <- unique(sub("_(barcodes|features|matrix)\\.(tsv|mtx)\\.gz$", "", basename(members)))
  prefix <- prefix[!grepl("\\.expression_matrix\\.txt\\.gz$", prefix)]
  prefix <- prefix[nzchar(prefix)]

  # drop entries that are not sample prefixes (no GSM at the start)
  prefix <- prefix[grepl("^GSM[0-9]+", prefix)]

  if (!length(prefix)) {

    # per-sample expression-matrix layout: the prefix is the GSM + title
    prefix <- sub("\\.expression_matrix\\.txt\\.gz$", "", basename(members))
    prefix <- prefix[grepl("^GSM[0-9]+", prefix)]
  }

  sid_gsm <- sub("_.*$", "", prefix)

  idx <- match(sid_gsm, sm$gsm)

  rows[[acc]] <- data.frame(
    sample_id = prefix, gsm = sid_gsm, title = sm$title[idx],
    patient_id = sub("^(GSM[0-9]+)_", "", prefix),
    disease_group = canonical(sm$raw_characteristic[idx]),
    source_accession = acc, raw_characteristic = sm$raw_characteristic[idx],
    stringsAsFactors = FALSE)

  cat(acc, ":", nrow(rows[[acc]]), "samples |",
      paste(names(table(rows[[acc]]$disease_group)), table(rows[[acc]]$disease_group),
            sep = "=", collapse = ", "), "\n")
}

meta <- do.call(rbind, rows)

write.table(meta, file.path(scrna_dir, "sample_metadata.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

cat("\ntotal samples with metadata:", nrow(meta), "\n")
print(table(meta$source_accession, meta$disease_group))

cat("\nunlabelled samples (raw characteristic missing):", sum(is.na(meta$disease_group)), "\n")

#!/usr/bin/env Rscript
# Paper 1 - Phase A step 2c: normalise the three RNA-seq cohorts.
#
#   GSE182380 (TI,  NEPTUNE, symbol-keyed)   raw counts -> filterByExpr + TMM + voom
#   GSE197307 (Glom, NEPTUNE, Ensembl-keyed) raw counts -> filterByExpr + TMM + voom
#   GSE216841 (Glom, LCM-FFPE, Ensembl)      values are already normalised decimals, so they are
#                                            log2-transformed rather than re-counted
#
# Every cohort is normalised on its own and never merged with the microarrays (two-layer design).

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(edgeR); library(limma) })

out_dir <- file.path(PROJ, "work", "rnaseq")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_counts <- function(path) {
  d <- read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
  ids <- as.character(d[[1]])
  d <- d[, -1, drop = FALSE]
  # some files carry annotation columns beside the gene id (GSE216841 ships hgnc symbol and
  # description), so keep only the columns that really are numeric sample columns
  numeric_col <- vapply(d, function(x) {
    y <- suppressWarnings(as.numeric(x))
    sum(!is.na(y)) > 0.9 * length(y)
  }, logical(1))
  dropped <- names(d)[!numeric_col]
  if (length(dropped)) cat("  dropping non-numeric columns:", paste(dropped, collapse = ", "), "\n")
  d <- d[, numeric_col, drop = FALSE]
  m <- as.matrix(d)
  mode(m) <- "numeric"
  rownames(m) <- ids
  m
}

save_mat <- function(m, tag) {
  f <- file.path(out_dir, paste0(tag, "_log2expr.txt.gz"))
  gz <- gzfile(f, "w")
  write.table(data.frame(gene = rownames(m), m, check.names = FALSE), gz, sep = "\t",
              quote = FALSE, row.names = FALSE)
  close(gz)
  saveRDS(m, file.path(out_dir, paste0(tag, "_log2expr.rds")))
  cat(sprintf("  %-16s %6d genes x %4d samples -> %s\n", tag, nrow(m), ncol(m), basename(f)))
}

# ---- integer-count cohorts: filterByExpr -> TMM -> voom -----------------------------------------

counts_to_log2 <- function(m, tag) {
  cat(tag, ": raw counts", nrow(m), "x", ncol(m),
      "| integer-looking:", all(abs(m - round(m)) < 1e-6), "\n")
  dge <- DGEList(counts = round(m))
  keep <- filterByExpr(dge)
  cat("  filterByExpr keeps", sum(keep), "of", nrow(m), "genes\n")
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  cat("  TMM norm factors range:", paste(round(range(dge$samples$norm.factors), 3), collapse = " .. "), "\n")
  design <- model.matrix(~1, data = data.frame(row.names = colnames(dge)))
  v <- voom(dge, design, plot = FALSE)
  cat("  voom log2-CPM range:", paste(round(range(v$E), 2), collapse = " .. "), "\n")
  saveRDS(dge, file.path(out_dir, paste0(tag, "_dge.rds")))
  v$E
}

# ---- GSE182380 (tubulointerstitium, NEPTUNE) ----------------------------------------------------

cat("=== RNA-seq normalisation ===\n")

m182 <- read_counts(file.path(DATA, "GSE182380", "processed", "GSE182380_counts_symbol.txt.gz"))

save_mat(counts_to_log2(m182, "GSE182380"), "GSE182380")

# ---- GSE197307 (glomeruli, NEPTUNE) -------------------------------------------------------------

m197 <- read_counts(file.path(DATA, "GSE197307", "processed", "GSE197307_counts_ensg.txt.gz"))

save_mat(counts_to_log2(m197, "GSE197307"), "GSE197307")

# ---- GSE216841 (glomeruli, LCM-FFPE) ------------------------------------------------------------

m216 <- read_counts(file.path(DATA, "GSE216841", "processed", "GSE216841_counts_ensg.txt.gz"))

cat("GSE216841: values", nrow(m216), "x", ncol(m216),
    "| integer-looking:", all(abs(m216 - round(m216)) < 1e-6),
    "| col-sum range:", paste(round(range(colSums(m216, na.rm = TRUE))), collapse = " .. "), "\n")

# the file is named "ncounts" but the values are normalised decimals, so treat it as a normalised
# matrix: drop all-zero genes, then log2.

m216 <- m216[rowSums(m216, na.rm = TRUE) > 0, ]

save_mat(log2(m216 + 1), "GSE216841")

cat("\n=== RNA-seq outputs ===\n")

print(list.files(out_dir))

writeLines(c(
  paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 02c_rnaseq_normalize.R"),
  paste0("  ", list.files(out_dir))
), file.path(PROJ, "logs", "rnaseq_normalize.log"))

cat("\nDone.\n")

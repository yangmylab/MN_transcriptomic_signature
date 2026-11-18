#!/usr/bin/env Rscript
# Paper 1 - scRNA arm, step 07b: per-sample import, QC, doublet removal and normalisation.
#
# Inputs : scrna/samples.csv (07a), one row per sample (10x triplet dir or dense expression matrix)
# Outputs: work/scrna/qc/<accession>__<sample>.rds   filtered, log-normalised Seurat object
#          work/scrna/qc/<accession>__<sample>_qc.rds per-sample QC row (makes the step resumable)
#          output/tables/scrna_qc.csv
#          output/scrna/07b_*.{pdf,jpeg}
#
# Filtering: MAD-based (3 MAD) on log1p(nCount), log1p(nFeature) and percent.mt, plus a hard
# percent.mt ceiling of 20 %. Doublets: scDblFinder, parallelised inside each sample.
# The step is resumable: a sample whose object and QC row already exist is skipped, which matters
# because 49 samples at ~2 min each do not fit in one foreground call.

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(Seurat); library(SingleCellExperiment); library(scDblFinder);
                                 library(BiocParallel); library(ggplot2) })

scrna_dir <- file.path(PROJ, "scrna")
qc_dir <- file.path(PROJ, "work", "scrna", "qc")
tab_dir <- file.path(OUT, "tables")
fig_dir <- file.path(OUT, "scrna")

dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, stem, width = 9, height = 5.5) {
  if (interactive()) print(p)
  ggsave(paste0(stem, ".pdf"), plot = p, device = "pdf", width = width, height = height)
  ggsave(paste0(stem, ".jpeg"), plot = p, device = "jpeg", quality = 95, dpi = 300,
         bg = "white", width = width, height = height)
}

if (!file.exists(file.path(scrna_dir, "samples.csv"))) {
  stop("scrna/samples.csv missing - place the archives and run R/07a_scrna_manifest.R first")
}

samples <- read.csv(file.path(scrna_dir, "samples.csv"), stringsAsFactors = FALSE)

if (!nrow(samples)) stop("scrna/samples.csv has no rows")

cat("samples to process:", nrow(samples), "\n")

# --- helpers -----------------------------------------------------------------

mad_bounds <- function(x, nmads = 3) {
  m <- median(x, na.rm = TRUE)
  s <- mad(x, na.rm = TRUE)
  c(lower = m - nmads * s, upper = m + nmads * s)
}

# samples.csv stores paths relative to PROJ; older manifests carry absolute paths from a previous
# mount point, so fall back to rebasing anything under .../extracted/ onto this project.
resolve_path <- function(p) {
  if (file.exists(p)) return(p)
  alt <- file.path(PROJ, p)
  if (file.exists(alt)) return(alt)
  rebased <- file.path(PROJ, "scrna", sub("^.*?/extracted/", "extracted/", p))
  if (file.exists(rebased)) return(rebased)
  p
}

read_sample <- function(path, format = "10x", prefix = NA_character_) {
  if (identical(format, "expression_matrix")) {
    d <- read.delim(gzfile(path), header = TRUE, row.names = 1, check.names = FALSE)
    m <- as.matrix(d)
    # GEO dense matrices are genes x cells; transpose only if the shape says otherwise
    if (ncol(m) > nrow(m)) m <- t(m)
    if (any(m < 0, na.rm = TRUE)) stop("negative values in ", path)
    if (!isTRUE(all.equal(m, round(m)))) {
      warning("non-integer expression matrix in ", path, " - rounded for count-based QC")
      m <- round(m)
    }
    return(Matrix::Matrix(m, sparse = TRUE))
  }
  if (grepl("\\.h5$", path)) return(Read10X_h5(path))
  # 10x: the whole dataset is untarred into ONE directory and `prefix` selects the triplet.
  # This replaced the old per-sample folder of symlinks, which cannot exist on an exfat volume
  # (see P-34) - 07a now records dir + prefix instead of building links.
  if (!is.na(prefix) && nzchar(prefix)) {
    f <- list.files(path, pattern = paste0("^", prefix, "_.*\\.gz$"), full.names = TRUE)
    mtx <- f[grepl("_matrix\\.mtx\\.gz$", f)]
    bar <- f[grepl("_barcodes\\.tsv\\.gz$", f)]
    fea <- f[grepl("_(features|genes)\\.tsv\\.gz$", f)]
    if (length(mtx) == 1 && length(bar) == 1 && length(fea) == 1) {
      return(ReadMtx(mtx = mtx, cells = bar, features = fea, feature.column = 2))
    }
    stop("expected exactly one barcodes/features/matrix file for prefix '", prefix, "' in ", path)
  }
  Read10X(path)
}

tag_of <- function(i) {
  paste0(samples$accession[i], "__", gsub("[^A-Za-z0-9]+", "_", samples$sample_id[i]))
}

process_sample <- function(i) {

  acc <- samples$accession[i]
  sid <- samples$sample_id[i]
  tag <- tag_of(i)

  f_obj <- file.path(qc_dir, paste0(tag, ".rds"))
  f_qc <- file.path(qc_dir, paste0(tag, "_qc.rds"))

  if (file.exists(f_obj) && file.exists(f_qc)) {
    cat("  [resume]", tag, "\n")
    return(readRDS(f_qc))
  }

  fmt <- if ("format" %in% names(samples)) samples$format[i] else "10x"

  # 07a records the untarred directory plus the filename prefix for 10x datasets (no per-sample
  # link folders exist on an exfat volume); older samples.csv files have no prefix column.
  pfx <- if ("prefix" %in% names(samples)) samples$prefix[i] else NA_character_

  mat <- read_sample(resolve_path(samples$path[i]), fmt, pfx)

  if (is.list(mat)) stop("multi-modal 10x directory for ", sid)

  obj <- CreateSeuratObject(counts = mat, project = tag, min.cells = 3, min.features = 200)

  rm(mat)
  invisible(gc())

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  n0 <- ncol(obj)

  b_count <- mad_bounds(log1p(obj$nCount_RNA))
  b_feat <- mad_bounds(log1p(obj$nFeature_RNA))

  keep <- obj$nCount_RNA >= expm1(b_count["lower"]) &
    obj$nCount_RNA <= expm1(b_count["upper"]) &
    obj$nFeature_RNA >= expm1(b_feat["lower"]) &
    obj$nFeature_RNA <= expm1(b_feat["upper"]) &
    obj$percent.mt < 20

  obj <- obj[, keep]

  n_after_qc <- ncol(obj)

  # Doublets are called in the droplet-based scRNA cohorts only. In single-nucleus preparations the
  # droplet-doublet expectation is not calibrated (nuclei are not encapsulated at the same rate), so
  # snRNA nuclei are retained without a doublet filter and the choice is recorded per sample.
  lay <- if ("layer" %in% names(samples)) samples$layer[i] else "kidney_scRNA"

  if (identical(lay, "kidney_snRNA")) {
    obj$doublet <- FALSE
    n_doublet <- 0L
    dbl_method <- "not applied (single-nucleus)"
  } else {
    sce <- as.SingleCellExperiment(obj)
    sce <- scDblFinder(sce, dbr = 0.04, verbose = FALSE,
                       BPPARAM = MulticoreParam(workers = 6, progressbar = FALSE))
    obj$doublet <- sce$scDblFinder.class == "doublet"
    n_doublet <- sum(obj$doublet)
    dbl_method <- "scDblFinder dbr=0.04"
  }

  obj <- obj[, !obj$doublet]

  n_kept <- ncol(obj)

  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)

  obj$sample_id <- sid
  obj$accession <- acc
  obj$dataset <- samples$dataset[i]
  obj$layer <- samples$layer[i]
  obj$donor <- if ("patient_id" %in% names(samples)) samples$patient_id[i] else sid
  obj$disease_group <- if ("disease_group" %in% names(samples)) samples$disease_group[i] else NA

  qc_row <- data.frame(
    accession = acc, sample_id = sid, cells_raw = n0,
    removed_count_or_feature = n0 - n_after_qc,
    removed_doublet = n_doublet, cells_kept = n_kept, doublet_method = dbl_method,
    median_nFeature = round(median(obj$nFeature_RNA), 1),
    median_nCount = round(median(obj$nCount_RNA), 1),
    median_percent_mt = round(median(obj$percent.mt), 2),
    count_lower = round(expm1(b_count["lower"]), 1), count_upper = round(expm1(b_count["upper"]), 1),
    feature_lower = round(expm1(b_feat["lower"]), 1), feature_upper = round(expm1(b_feat["upper"]), 1),
    percent_mt_cut = 20, stringsAsFactors = FALSE)

  saveRDS(obj, f_obj)
  saveRDS(qc_row, f_qc)

  cat("  ", tag, ": kept", n_kept, "of", n0, "cells (doublets", n_doublet, ")\n")

  qc_row
}

# --- run ---------------------------------------------------------------------

# samples run sequentially: two concurrent Seurat+SCE objects plus their doublet workers exceeded
# the 14 GB budget and a worker was killed (SIGPIPE). Parallelism is inside scDblFinder instead.
res <- lapply(seq_len(nrow(samples)), function(i) {
  out <- tryCatch(process_sample(i), error = function(e) conditionMessage(e))
  invisible(gc())
  out
})

errs <- !vapply(res, is.data.frame, logical(1))

if (any(errs)) {
  cat("\nsamples that failed:\n")
  for (i in which(errs)) cat("  ", samples$sample_id[i], ":", as.character(res[[i]]), "\n")
  stop(sum(errs), " sample(s) failed - fix before continuing")
}

qc <- do.call(rbind, res)

write.csv(qc, file.path(tab_dir, "scrna_qc.csv"), row.names = FALSE)

# cells table: rebuilt from the saved objects so resumed samples are included
cells <- do.call(rbind, lapply(seq_len(nrow(samples)), function(i) {
  o <- readRDS(file.path(qc_dir, paste0(tag_of(i), ".rds")))
  data.frame(cell = colnames(o), sample_id = o$sample_id[1], donor = o$donor[1],
             accession = o$accession[1], disease_group = o$disease_group[1],
             celltype = NA_character_, stringsAsFactors = FALSE)
}))

write.csv(cells, file.path(PROJ, "work", "scrna", "cells_qc.csv"), row.names = FALSE)

# --- QC figure ---------------------------------------------------------------

long <- data.frame(
  sample = rep(qc$sample_id, 3),
  metric = rep(c("cells kept", "median nFeature", "median percent.mt"), each = nrow(qc)),
  value = c(qc$cells_kept, qc$median_nFeature, qc$median_percent_mt))

p <- ggplot(long, aes(value, reorder(sample, value))) +
  geom_col(fill = "#2471A3") +
  facet_wrap(~ metric, scales = "free_x") +
  labs(title = "Per-sample QC (after filtering)", x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 6))

save_plot(p, file.path(fig_dir, "07b_sample_qc"), width = 11, height = 9)

cat("\n================ scRNA QC ================\n")
print(qc, row.names = FALSE)

cat("\ntotal cells kept:", sum(qc$cells_kept), "\n")

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 07b_scrna_qc.R"),
             capture.output(print(qc, row.names = FALSE))),
           file.path(PROJ, "logs", "scrna_qc.log"))

cat("\nDone.\n")

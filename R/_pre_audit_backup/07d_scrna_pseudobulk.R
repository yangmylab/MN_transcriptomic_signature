#!/usr/bin/env Rscript
# Paper 1 - scRNA arm, step 07d: pseudobulk per donor x cell type, differential expression and
# composition.
#
# Why pseudobulk: with 1-6 donors per group, cell-level tests are anticonservative. Counts are summed
# per donor and cell type, then modelled with limma-voom; the donor is the unit of replication.
#
# Inputs : work/scrna/qc/*.rds (07b, all genes), work/scrna/atlas_cells.csv (07c)
# Outputs: work/scrna/pseudobulk/<celltype>.rds        genes x donors count matrix + donor table
#          output/tables/scrna_pseudobulk_de.csv        long table of all per-cell-type contrasts
#          output/tables/scrna_detection_fraction.csv  gene x cell type detection fraction
#          output/tables/scrna_composition.csv          cell-type proportions per donor
#          output/tables/scrna_composition_tests.csv    CLR + Wilcoxon MN vs control
#          output/scrna/07d_*.{pdf,jpeg}

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(Seurat); library(edgeR); library(limma); library(ggplot2) })

qc_dir <- file.path(PROJ, "work", "scrna", "qc")
pb_dir <- file.path(PROJ, "work", "scrna", "pseudobulk")
tab_dir <- file.path(OUT, "tables")
fig_dir <- file.path(OUT, "scrna")

dir.create(pb_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, stem, width = 9, height = 5.5) {
  if (interactive()) print(p)
  ggsave(paste0(stem, ".pdf"), plot = p, device = "pdf", width = width, height = height)
  ggsave(paste0(stem, ".jpeg"), plot = p, device = "jpeg", quality = 95, dpi = 300,
         bg = "white", width = width, height = height)
}

cells <- read.csv(file.path(PROJ, "work", "scrna", "atlas_cells.csv"), stringsAsFactors = FALSE)

if (!nrow(cells)) stop("atlas_cells.csv is empty - run 07c first")

# 07b/07c label healthy donors "Control" (capital C), while every comparison below is written
# against the lower-case "control" level. Normalise once, here, instead of silently skipping
# every contrast that needs the control group.
cells$disease_group[cells$disease_group == "Control"] <- "control"

MIN_CELLS <- 10

files <- list.files(qc_dir, pattern = "\\.rds$", full.names = TRUE)

# the per-sample *_qc.rds rows live in the same directory and are not Seurat objects
files <- files[!grepl("_qc\\.rds$", files)]

# --- aggregate counts per donor x cell type ----------------------------------

agg <- list()        # celltype -> list(genes, matrix, donors)
detect <- list()     # celltype -> named numeric detection fraction

for (f in files) {

  obj <- readRDS(f)

  # atlas cell names are the sample tag (accession__sample_id) + "__" + the original barcode: 07c
  # renames them so that barcodes shared between samples stay traceable. Rebuild the same key here.
  uid <- paste0(sub("\\.rds$", "", basename(f)), "__", colnames(obj))

  keep_cells <- which(uid %in% cells$cell)

  if (!length(keep_cells)) { rm(obj); next }

  cm <- cells[match(uid[keep_cells], cells$cell), ]

  counts <- LayerData(obj, assay = "RNA", layer = "counts")[, keep_cells, drop = FALSE]

  donor <- cm$donor[1]
  acc <- cm$accession[1]
  disease <- cm$disease_group[1]

  for (ct in unique(cm$celltype)) {

    sel <- which(cm$celltype == ct)

    if (length(sel) < MIN_CELLS) next

    s <- Matrix::rowSums(counts[, sel, drop = FALSE])

    agg[[ct]][[donor]] <- s

    detect[[ct]][[donor]] <- Matrix::rowMeans(counts[, sel, drop = FALSE] > 0)

  }

  rm(obj, counts)
  invisible(gc())

  cat("  aggregated", basename(f), "\n")
}

# --- build matrices and run DE ----------------------------------------------

all_ct <- names(agg)

donor_info <- unique(cells[, c("donor", "accession", "disease_group")])

de_rows <- list()
detect_rows <- list()

for (ct in all_ct) {

  donors <- names(agg[[ct]])

  genes <- Reduce(intersect, lapply(agg[[ct]], names))

  mat <- do.call(cbind, lapply(agg[[ct]], function(v) v[genes]))

  colnames(mat) <- donors

  info <- donor_info[match(donors, donor_info$donor), ]

  saveRDS(list(counts = mat, donors = info), file.path(pb_dir, paste0(ct, ".rds")))

  detect_rows[[ct]] <- data.frame(
    celltype = ct, gene = genes,
    detection_fraction = rowMeans(do.call(cbind, lapply(detect[[ct]], function(v) v[genes])),
                                  na.rm = TRUE),
    n_donors = length(donors), stringsAsFactors = FALSE)

  # contrast: each disease vs control, inside the cohort that has controls
  grp <- info$disease_group

  if (!any(grp == "control", na.rm = TRUE)) next

  grp[is.na(grp)] <- "unknown"

  levels_needed <- intersect(c("control", "MN", "IgAN", "MCD", "DKD"), unique(grp))

  if (length(levels_needed) < 2) next

  keep <- grp %in% levels_needed

  mat_f <- mat[, keep, drop = FALSE]
  grp_f <- factor(grp[keep], levels = levels_needed)

  if (min(table(grp_f)) < 2) next

  y <- DGEList(mat_f)

  keep_g <- filterByExpr(y, group = grp_f)

  y <- y[keep_g, ]

  y <- calcNormFactors(y)

  design <- model.matrix(~ 0 + grp_f)
  colnames(design) <- levels(grp_f)

  v <- voom(y, design, plot = FALSE)

  fit <- eBayes(lmFit(v, design), trend = TRUE, robust = TRUE)

  for (g in setdiff(levels(grp_f), "control")) {

    cm <- makeContrasts(contrasts = paste(g, "- control"), levels = design)

    fitc <- eBayes(contrasts.fit(fit, cm), trend = TRUE, robust = TRUE)

    tt <- topTable(fitc, coef = 1, number = Inf, sort.by = "P")

    de_rows[[paste(ct, g)]] <- data.frame(
      celltype = ct, contrast = paste0(g, "_vs_control"), n_case = sum(grp_f == g),
      n_control = sum(grp_f == "control"), gene = rownames(tt), logFC = tt$logFC,
      P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val, AveExpr = tt$AveExpr,
      stringsAsFactors = FALSE)
  }
}

if (length(de_rows)) {

  de <- do.call(rbind, de_rows)
  write.csv(de, file.path(tab_dir, "scrna_pseudobulk_de.csv"), row.names = FALSE)

  cat("\npseudobulk DE rows:", nrow(de), "| cell types:", length(unique(de$celltype)), "\n")
  print(table(de$celltype, de$contrast))
} else {
  cat("\nno cell type had >= 2 donors per group - pseudobulk DE not run\n")
}

if (length(detect_rows)) {
  write.csv(do.call(rbind, detect_rows), file.path(tab_dir, "scrna_detection_fraction.csv"),
            row.names = FALSE)
}

# --- composition -------------------------------------------------------------

comp <- as.data.frame(table(donor = cells$donor, celltype = cells$celltype))

tot <- aggregate(Freq ~ donor, data = comp, sum)

comp$proportion <- comp$Freq / tot$Freq[match(comp$donor, tot$donor)]

comp <- merge(comp, unique(cells[, c("donor", "disease_group")]), by = "donor")

write.csv(comp, file.path(tab_dir, "scrna_composition.csv"), row.names = FALSE)

comp_wide <- reshape(comp[, c("donor", "celltype", "proportion", "disease_group")],
                     idvar = c("donor", "disease_group"), timevar = "celltype",
                     direction = "wide")

props <- as.matrix(comp_wide[, grep("^proportion\\.", names(comp_wide)), drop = FALSE])
props[is.na(props)] <- 0

clr <- log(props + 1e-6)
clr <- clr - rowMeans(clr)

test_rows <- list()

if (any(comp$disease_group == "control", na.rm = TRUE)) {

  for (ct in unique(comp$celltype)) {

    x <- comp$proportion[comp$celltype == ct & comp$disease_group == "MN"]
    y <- comp$proportion[comp$celltype == ct & comp$disease_group == "control"]

    if (length(x) < 2 || length(y) < 2) next

    test_rows[[ct]] <- data.frame(
      celltype = ct, median_MN = round(median(x), 4), median_control = round(median(y), 4),
      p_value = round(wilcox.test(x, y)$p.value, 4), n_MN = length(x), n_control = length(y),
      stringsAsFactors = FALSE)
  }
}

if (length(test_rows)) {
  tests <- do.call(rbind, test_rows)
  tests$p_adj <- round(p.adjust(tests$p_value, method = "BH"), 4)
  write.csv(tests, file.path(tab_dir, "scrna_composition_tests.csv"), row.names = FALSE)
  print(tests, row.names = FALSE)
}

p <- ggplot(comp[!is.na(comp$disease_group), ], aes(donor, proportion, fill = celltype)) +
  geom_col(width = 0.85) +
  facet_grid(~ disease_group, scales = "free_x", space = "free_x") +
  labs(title = "Cell-type composition per donor", x = NULL, y = "proportion of cells") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p, file.path(fig_dir, "07d_composition"), width = 12, height = 6)

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 07d_scrna_pseudobulk.R"),
             paste("cell types aggregated:", length(all_ct))),
           file.path(PROJ, "logs", "scrna_pseudobulk.log"))

cat("\nDone.\n")

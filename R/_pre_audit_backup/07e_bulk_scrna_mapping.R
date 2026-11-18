#!/usr/bin/env Rscript
# Paper 1 - scRNA arm, step 07e: the bridge between the bulk signature and the single-cell atlas.
#
# Three products:
#   1. signature_matrix.csv               cell-type mean expression (CPM), CIBERSORTx-compatible
#   2. bulk_deg_celltype_attribution.csv  top cell type + specificity (tau) + detection for every
#                                         meta-DEG and every gate feature
#   3. bulk_vs_pseudobulk_concordance.csv per cell type, bulk MN-vs-control logFC vs pseudobulk logFC
#
# Inputs : work/scrna/pseudobulk/*.rds (07d), output/tables/meta/*_metaDEG.csv,
#          output/tables/separability_selected_genes.csv, output/tables/deg/*.csv (04b)
# Outputs: output/tables/signature_matrix.csv, bulk_deg_celltype_attribution.csv,
#          panel_celltype_specificity.csv, bulk_vs_pseudobulk_concordance.csv,
#          output/scrna/07e_*.{pdf,jpeg}

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(edgeR); library(ggplot2) })

pb_dir <- file.path(PROJ, "work", "scrna", "pseudobulk")
tab_dir <- file.path(OUT, "tables")
fig_dir <- file.path(OUT, "scrna")
set_dir <- file.path(PROJ, "work", "sets")

save_plot <- function(p, stem, width = 9, height = 5.5) {
  if (interactive()) print(p)
  ggsave(paste0(stem, ".pdf"), plot = p, device = "pdf", width = width, height = height)
  ggsave(paste0(stem, ".jpeg"), plot = p, device = "jpeg", quality = 95, dpi = 300,
         bg = "white", width = width, height = height)
}

files <- list.files(pb_dir, pattern = "\\.rds$", full.names = TRUE)

if (!length(files)) stop("no pseudobulk objects - run 07d first")

# --- 1. signature matrix -----------------------------------------------------

cpm_list <- list()

for (f in files) {

  ct <- sub("\\.rds$", "", basename(f))
  pb <- readRDS(f)

  y <- calcNormFactors(DGEList(pb$counts))

  cpm <- cpm(y, log = FALSE)

  cpm_list[[ct]] <- rowMeans(cpm)
}

# The cell types do NOT share a gene universe: each pseudobulk keeps only the genes detected in the
# samples where that cell type was found (13 868 - 28 342 genes here). do.call(cbind, ...) would
# then align the vectors BY POSITION and silently scramble every column after the first length
# mismatch (it warns "number of rows of result is not a multiple of vector length"), so align on the
# union of gene names and fill 0 for a gene a cell type never had.
all_genes <- sort(unique(unlist(lapply(cpm_list, names))))

sig <- vapply(cpm_list, function(v) {
  out <- numeric(length(all_genes))
  names(out) <- all_genes
  out[names(v)] <- v
  out
}, numeric(length(all_genes)))

rownames(sig) <- all_genes

cat("aligned signature:", nrow(sig), "genes x", ncol(sig), "cell types\n")

# genes to keep: those measurable in the bulk cohorts (symbol keys)
bulk_genes <- rownames(readRDS(file.path(set_dir, "disc_glom_array_expr_symbol.rds")))

sig <- sig[intersect(rownames(sig), bulk_genes), , drop = FALSE]

sig <- sig[rowSums(sig) > 0, , drop = FALSE]

write.csv(round(as.data.frame(sig), 4), file.path(tab_dir, "signature_matrix.csv"))

cat("signature matrix:", nrow(sig), "genes x", ncol(sig), "cell types\n")

# --- 2. cell-type attribution ------------------------------------------------

tau <- function(x) {
  n <- length(x)
  if (n < 2 || all(x == 0)) return(0)
  sum(1 - x / max(x)) / (n - 1)
}

attribute <- function(genes) {
  g <- intersect(genes, rownames(sig))
  if (!length(g)) return(NULL)
  m <- sig[g, , drop = FALSE]
  top <- colnames(m)[max.col(m, ties.method = "first")]
  data.frame(gene = g, top_celltype = top,
             tau = round(apply(m, 1, tau), 4),
             max_cpm = round(apply(m, 1, max), 2),
             second_celltype = colnames(m)[apply(m, 1, function(z) order(z, decreasing = TRUE)[2])],
             stringsAsFactors = FALSE)
}

attr_rows <- list()

for (f in list.files(file.path(tab_dir, "meta"), pattern = "_metaDEG\\.csv$", full.names = TRUE)) {

  nm <- sub("_metaDEG\\.csv$", "", basename(f))
  d <- read.csv(f, stringsAsFactors = FALSE)
  a <- attribute(d$gene)
  if (is.null(a)) next
  a$set <- nm
  a$pooled_logFC <- d$logFC[match(a$gene, d$gene)]
  attr_rows[[nm]] <- a
}

if (length(attr_rows)) {
  attr <- do.call(rbind, attr_rows)
  write.csv(attr, file.path(tab_dir, "bulk_deg_celltype_attribution.csv"), row.names = FALSE)

  cat("\nmeta-DEGs attributed:", nrow(attr), "rows\n")
  print(head(sort(table(attr$set, attr$top_celltype), decreasing = TRUE), 20))
}

gate <- read.csv(file.path(tab_dir, "separability_selected_genes.csv"), stringsAsFactors = FALSE)
gate <- gate[gate$K == "all", ]

if (nrow(gate)) {
  pg <- attribute(gate$gene)
  if (!is.null(pg)) {
    write.csv(pg, file.path(tab_dir, "panel_celltype_specificity.csv"), row.names = FALSE)
    cat("\ngate features attributed:", nrow(pg), "\n")
    print(sort(table(pg$top_celltype), decreasing = TRUE))
  }
}

# --- 3. bulk vs pseudobulk concordance ---------------------------------------

bulk_deg <- read.csv(file.path(tab_dir, "deg", "disc_glom_array__MN_vs_Control.csv"),
                     stringsAsFactors = FALSE)

pb_de <- read.csv(file.path(tab_dir, "scrna_pseudobulk_de.csv"), stringsAsFactors = FALSE)
pb_de <- pb_de[pb_de$contrast == "MN_vs_control", ]

conc_rows <- list()

if (nrow(pb_de)) {

  for (ct in unique(pb_de$celltype)) {

    d <- pb_de[pb_de$celltype == ct, ]

    common <- intersect(d$gene, bulk_deg$gene)

    if (length(common) < 100) next

    x <- bulk_deg$logFC[match(common, bulk_deg$gene)]
    y <- d$logFC[match(common, d$gene)]

    conc_rows[[ct]] <- data.frame(
      celltype = ct, genes = length(common),
      pearson = round(cor(x, y), 3), spearman = round(cor(x, y, method = "spearman"), 3),
      direction_agreement = round(100 * mean(sign(x) == sign(y)), 1),
      stringsAsFactors = FALSE)
  }
}

if (length(conc_rows)) {
  conc <- do.call(rbind, conc_rows)
  conc <- conc[order(-conc$pearson), ]
  write.csv(conc, file.path(tab_dir, "bulk_vs_pseudobulk_concordance.csv"), row.names = FALSE)
  print(conc, row.names = FALSE)

  dd <- pb_de[pb_de$celltype == conc$celltype[1], ]
  common <- intersect(dd$gene, bulk_deg$gene)

  plot_df <- data.frame(
    bulk_logFC = bulk_deg$logFC[match(common, bulk_deg$gene)],
    pseudobulk_logFC = dd$logFC[match(common, dd$gene)],
    sig = bulk_deg$deg[match(common, bulk_deg$gene)])

  p <- ggplot(plot_df, aes(bulk_logFC, pseudobulk_logFC, colour = sig)) +
    geom_point(size = 0.7, alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "grey70"),
                        labels = c(`TRUE` = "bulk DEG", `FALSE` = "not a bulk DEG"),
                        name = NULL) +
    labs(title = paste0("Bulk vs pseudobulk (", conc$celltype[1], ", MN vs control)"),
         subtitle = paste0("Pearson r = ", conc$pearson[1],
                           "; direction agreement ", conc$direction_agreement[1], "%"),
         x = "bulk log2FC", y = "pseudobulk log2FC") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  save_plot(p, file.path(fig_dir, "07e_bulk_vs_pseudobulk"), width = 6.5, height = 5.5)
}

# attribution figure: top cell type counts for the glomerular control meta-DEGs
if (length(attr_rows) && "Glom_MN_vs_Control" %in% names(attr_rows)) {

  a <- attr_rows[["Glom_MN_vs_Control"]]

  p <- ggplot(a, aes(reorder(top_celltype, -tau), tau)) +
    geom_boxplot(outlier.shape = NA, fill = "#D6EAF8") +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.6, colour = "#2471A3") +
    labs(title = "Cell-type attribution of glomerular MN-vs-control meta-DEGs",
         subtitle = "tau = cell-type specificity of the top cell type (1 = exclusive)",
         x = NULL, y = "tau") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  save_plot(p, file.path(fig_dir, "07e_metaDEG_attribution"), width = 8, height = 5)
}

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 07e_bulk_scrna_mapping.R"),
             paste("signature matrix:", nrow(sig), "x", ncol(sig))),
           file.path(PROJ, "logs", "scrna_mapping.log"))

cat("\nDone.\n")

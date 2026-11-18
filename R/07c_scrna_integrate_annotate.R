#!/usr/bin/env Rscript
# Paper 1 - scRNA arm, step 07c: integrate the per-sample objects and annotate cell types.
#
# Inputs : work/scrna/qc/*.rds (07b)
# Outputs: work/scrna/atlas.rds                 Harmony-integrated object on the integrated HVG set
#          work/scrna/atlas_cells.csv           cell -> sample / donor / disease / cell type / confidence
#          output/tables/scrna_celltype_counts.csv
#          output/scrna/07c_*.{pdf,jpeg}
#
# Memory strategy: the per-sample objects keep ALL genes (07d needs them for pseudobulk), but the
# atlas is built on a reduced gene set - genes that are variable in >= 10 % of samples (max 3000)
# plus the marker panel. That keeps a ~350 k-cell object inside the 14 GB budget.

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(Seurat); library(harmony); library(ggplot2) })

qc_dir <- file.path(PROJ, "work", "scrna", "qc")
out_dir <- file.path(PROJ, "work", "scrna")
tab_dir <- file.path(OUT, "tables")
fig_dir <- file.path(OUT, "scrna")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, stem, width = 9, height = 6) {
  if (interactive()) print(p)
  ggsave(paste0(stem, ".pdf"), plot = p, device = "pdf", width = width, height = height)
  ggsave(paste0(stem, ".jpeg"), plot = p, device = "jpeg", quality = 95, dpi = 300,
         bg = "white", width = width, height = height)
}

files <- list.files(qc_dir, pattern = "\\.rds$", full.names = TRUE)

# the per-sample *_qc.rds rows live in the same directory and must not be loaded as objects
files <- files[!grepl("_qc\\.rds$", files)]

if (!length(files)) stop("no QC objects in work/scrna/qc - run 07b first")

cat("QC objects:", length(files), "\n")

markers <- list(
  Podocyte = c("NPHS1", "NPHS2", "PODXL", "PLA2R1", "WT1"),
  PEC = c("PAX8", "CLDN1", "KRT19", "TACSTD2"),
  PT = c("LRP2", "CUBN", "SLC34A1", "SLC5A12"),
  TAL = c("UMOD", "SLC12A1", "CLDN16"),
  DCT = c("SLC12A3", "TRPM6", "CALB1"),
  CD_PC = c("AQP2", "AQP3", "FXYD4"),
  Endothelium = c("PECAM1", "EMCN", "FLT1", "CDH5"),
  Mesangium = c("PDGFRB", "ITGA8", "ACTA2", "NOTCH3"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM"),
  T_NK = c("CD3D", "CD3E", "IL7R", "NKG7", "GNLY"),
  B = c("MS4A1", "CD79A", "CD19"),
  Plasma = c("MZB1", "JCHAIN", "SDC1"),
  Myeloid = c("LYZ", "CD68", "C1QA", "CD14"),
  pDC = c("LILRA4", "CLEC4C"),
  Mast = c("TPSAB1", "CPA3", "KIT"),
  Neutrophil = c("FCGR3B", "S100A8", "S100A9"))

# --- choose the integrated gene set ------------------------------------------

# Memory: the 49 per-sample objects are never all held at once. First pass reads only the variable
# genes; second pass keeps a capped number of cells per sample on a reduced gene set, because
# ScaleData needs a dense matrix. The caps were originally 3 000 cells / 1 500 genes to fit the
# 14 GB laptop; the current workstation has 61 GB, so they were raised (user decision 2026-10-19)
# to 6 000 cells / 2 000 genes => ~200 k cells x ~2 050 genes, still ~3.5 GB dense.
MAX_CELLS <- 6000
MAX_GENES <- 2000

hvg_tab <- table(unlist(lapply(files, function(f) VariableFeatures(readRDS(f)))))

min_hits <- max(2, ceiling(0.1 * length(files)))

hvg <- names(hvg_tab)[hvg_tab >= min_hits]

hvg <- hvg[order(-as.integer(hvg_tab[hvg]))][1:min(MAX_GENES, length(hvg))]

cat("HVG selected:", length(hvg), "genes\n")

set.seed(20261018)

sub <- list()

for (f in files) {
  o <- readRDS(f)
  cells <- colnames(o)
  if (length(cells) > MAX_CELLS) cells <- sample(cells, MAX_CELLS)
  g <- intersect(unique(c(hvg, unlist(markers))), rownames(o))
  nm <- sub("\\.rds$", "", basename(f))
  # 10x barcodes repeat across samples, so merge() would append _1 _2 ... to the same barcode from
  # different samples and the resulting names could no longer be traced back to their sample.
  # Prefix every barcode with the sample tag (accession__sample_id) instead: the atlas cell names
  # are then unique AND 07d can rebuild the same key to recover the full-gene counts.
  sub[[nm]] <- RenameCells(subset(o, cells = cells, features = g),
                           new.names = paste0(nm, "__", cells))
  rm(o)
  invisible(gc())
}

marker_present <- intersect(unique(unlist(markers)), rownames(sub[[1]]))

genes <- union(hvg, marker_present)

cat("integrated gene set:", length(genes), "genes (", length(hvg), "HVG +",
    length(marker_present), "markers ) | cells:",
    sum(vapply(sub, function(z) as.integer(ncol(z)), integer(1))), "\n")

# --- merge, scale, PCA, Harmony ----------------------------------------------

atlas <- merge(sub[[1]], y = sub[-1], project = "MN_scrna_atlas")

rm(sub)
invisible(gc())

# Seurat 5 keeps one `counts.<sample>` / `data.<sample>` layer per merged object, and each layer
# carries only that sample's genes. AddModuleScore iterates EVERY layer returned by
# Layers(search = "data") and dies with "feature lists do not have enough features present" as soon
# as one sample lacks a marker, so the layers must be joined into one before annotation.
# JoinLayers fills 0 for a gene a given sample never detected, which is what we want.
atlas <- JoinLayers(atlas)

atlas$donor_key <- paste(atlas$accession, atlas$donor, sep = "_")

atlas <- ScaleData(atlas, features = intersect(genes, rownames(atlas)), verbose = FALSE)
atlas <- RunPCA(atlas, features = intersect(hvg, rownames(atlas)), npcs = 30, verbose = FALSE)

atlas <- RunHarmony(atlas, group.by.vars = "donor_key", reduction = "pca",
                    reduction.save = "harmony", verbose = FALSE)

atlas <- FindNeighbors(atlas, reduction = "harmony", dims = 1:30, verbose = FALSE)
atlas <- FindClusters(atlas, resolution = 0.8, verbose = FALSE)
atlas <- RunUMAP(atlas, reduction = "harmony", dims = 1:30, verbose = FALSE)

# --- annotate with marker modules --------------------------------------------

score_sets <- lapply(markers, function(g) intersect(g, rownames(atlas)))

score_sets <- score_sets[vapply(score_sets, length, integer(1)) >= 2]

# The atlas carries only the integrated HVG + marker set, so the default control pool
# (ctrl = 100 genes per bin) is larger than a bin of 2 019 genes / 24 bins and AddModuleScore
# fails in sample.int(). Seurat draws `ctrl` control genes from EVERY expression bin, so the cap
# is the SMALLEST bin, not the average one - compute the bins exactly as Seurat does and keep a
# 40 % margin.
N_BIN <- 24L

avg_expr <- Matrix::rowMeans(LayerData(atlas, assay = "RNA", layer = "data"))

expr_bin <- ggplot2::cut_number(avg_expr + rnorm(length(avg_expr)) / 1e30, n = N_BIN,
                                labels = FALSE, right = FALSE)

min_bin <- min(tabulate(expr_bin))

ctrl_use <- as.integer(max(5L, floor(0.6 * min_bin)))

ctrl_use <- min(100L, ctrl_use)

cat("AddModuleScore control pool:", ctrl_use, "genes per bin x", N_BIN, "bins (smallest bin",
    min_bin, "genes)\n")

for (nm in names(score_sets)) {
  atlas <- AddModuleScore(atlas, features = list(score_sets[[nm]]), name = nm,
                          nbin = N_BIN, ctrl = ctrl_use, seed = 20261018)
}

score_cols <- paste0(names(score_sets), "1")

scores <- atlas@meta.data[, score_cols, drop = FALSE]
colnames(scores) <- names(score_sets)

best <- colnames(scores)[max.col(scores, ties.method = "first")]
second <- apply(scores, 1, function(z) sort(z, decreasing = TRUE)[2])

atlas$celltype <- best
atlas$celltype_confidence <- round(scores[cbind(seq_len(nrow(scores)), max.col(scores, ties.method = "first"))] -
                                     second, 4)

cell_meta <- atlas@meta.data[, c("celltype", "celltype_confidence", "sample_id", "donor",
                                  "accession", "disease_group", "seurat_clusters")]

cell_meta$cell <- rownames(cell_meta)

write.csv(cell_meta[, c("cell", "sample_id", "donor", "accession", "disease_group",
                        "celltype", "celltype_confidence", "seurat_clusters")],
          file.path(out_dir, "atlas_cells.csv"), row.names = FALSE)

saveRDS(atlas, file.path(out_dir, "atlas.rds"))

counts <- as.data.frame(table(sample = atlas$sample_id, celltype = atlas$celltype))
write.csv(counts, file.path(tab_dir, "scrna_celltype_counts.csv"), row.names = FALSE)

# --- figures -----------------------------------------------------------------

p <- DimPlot(atlas, group.by = "celltype", label = TRUE, repel = TRUE, raster = TRUE) +
  labs(title = "Integrated kidney atlas - cell types") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

save_plot(p, file.path(fig_dir, "07c_umap_celltype"), width = 8.5, height = 6.5)

p <- DimPlot(atlas, group.by = "disease_group", raster = TRUE) +
  labs(title = "Integrated kidney atlas - disease group") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

save_plot(p, file.path(fig_dir, "07c_umap_disease"), width = 8.5, height = 6.5)

p <- DotPlot(atlas, features = unique(unlist(score_sets)), group.by = "celltype") +
  RotatedAxis() +
  labs(title = "Canonical markers by cell type", x = NULL, y = NULL) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

save_plot(p, file.path(fig_dir, "07c_marker_dotplot"), width = 14, height = 6)

cat("\n================ CELL TYPES ================\n")
print(table(atlas$celltype))

cat("\nlow-confidence cells (score margin < 0.05):",
    sum(atlas$celltype_confidence < 0.05), "of", ncol(atlas), "\n")

cat("\nclusters:", length(unique(atlas$seurat_clusters)), "\n")

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 07c_scrna_integrate_annotate.R"),
             capture.output(print(table(atlas$celltype)))),
           file.path(PROJ, "logs", "scrna_atlas.log"))

cat("\nDone.\n")

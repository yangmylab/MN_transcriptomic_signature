#!/usr/bin/env Rscript
# Paper 1 - Phase F (exploratory): molecular subtypes within glomerular MN.
#
# Power is the binding constraint: glomerular MN samples are 51 (NEPTUNE array), 62 (NEPTUNE
# RNA-seq) and 21 (ERCB). IMPORTANT: the two NEPTUNE cohorts CANNOT be linked at patient level from
# the downloaded metadata - GSE197307/GSE182380 sample titles are local sample indices and the
# series matrices carry only disease + tissue, no participant id. Subtypes are therefore discovered
# WITHIN each platform and compared as *signatures*, not as labels on shared patients.
#
# Design
#   discovery-1 : disc_glom_array, MN only (n = 51)
#   discovery-2 : cross_glom_rnaseq, MN only (n = 62)
#   mapping     : val_glom_array, MN only (n = 21)   [nearest-centroid transfer, ERCB]
#   clustering  : ConsensusClusterPlus, ward.D2, Euclidean, top 2000 variable genes,
#                 z-scored within platform, k = 2..6, 1000 resamples (pItem 0.9, pFeature 0.8)
#   stability   : PAC (proportion of ambiguous clustering); PAC < 0.1 = stable partition
#   comparison  : Jaccard overlap of subtype-defining genes, cross-platform centroid correlation,
#                 and agreement of the two frameworks when transporting ERCB samples
#
# Everything here is labelled EXPLORATORY until clinical outcomes exist.
#
# Outputs
#   output/tables/subtype_stability.csv
#   output/tables/subtype_clusters.csv
#   output/tables/subtype_genes_<platform>.csv
#   output/tables/subtype_signature_concordance.csv
#   output/tables/subtype_ercb_assignment.csv
#   output/subtypes/06a_*.{pdf,jpeg}

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(ConsensusClusterPlus); library(limma); library(pheatmap);
                                 library(ggplot2); library(cluster) })

set_dir <- file.path(PROJ, "work", "sets")
tab_dir <- file.path(OUT, "tables")
fig_dir <- file.path(OUT, "subtypes")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

MAXK <- 6
REPS <- 1000
TOP_GENES <- 2000
PAC_STABLE <- 0.1
SEED <- 20261018

save_device_plot <- function(fn, stem, width = 7, height = 6) {
  while (!is.null(dev.list())) dev.off()
  dev.list()
  pdf(paste0(stem, ".pdf"), width = width, height = height)
  fn()
  dev.off()
  jpeg(paste0(stem, ".jpeg"), width = width, height = height, units = "in", res = 300,
       quality = 95, bg = "white")
  fn()
  dev.off()
}

save_plot <- function(p, stem, width = 7, height = 5.5) {
  if (interactive()) print(p)
  ggsave(paste0(stem, ".pdf"), plot = p, device = "pdf", width = width, height = height)
  ggsave(paste0(stem, ".jpeg"), plot = p, device = "jpeg", quality = 95, dpi = 300,
         bg = "white", width = width, height = height)
}

ari <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  if (n < 2) return(NA_real_)
  s_ij <- sum(choose(tab, 2))
  s_a <- sum(choose(rowSums(tab), 2))
  s_b <- sum(choose(colSums(tab), 2))
  exp <- s_a * s_b / choose(n, 2)
  mx <- (s_a + s_b) / 2
  (s_ij - exp) / (mx - exp)
}

mn_only <- function(set) {
  ms <- readRDS(file.path(set_dir, paste0(set, "_expr_symbol.rds")))
  ann <- read.csv(file.path(set_dir, paste0(set, "_samples.csv")), stringsAsFactors = FALSE)
  keep <- ann$disease == "MN"
  ms <- ms[, keep, drop = FALSE]
  ms <- ms[rowSums(is.na(ms)) == 0, , drop = FALSE]
  list(ms = ms, ann = ann[keep, , drop = FALSE])
}

prep <- function(ms, top = TOP_GENES) {
  v <- apply(ms, 1, var)
  g <- names(sort(v, decreasing = TRUE))[seq_len(min(top, length(v)))]
  z <- t(scale(t(ms[g, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  z
}

pac <- function(cons, lo = 0.1, hi = 0.9) mean(cons > lo & cons < hi)

run_ccp <- function(z, label) {

  set.seed(SEED)

  # plot = NULL is not honoured by this ConsensusClusterPlus build: it opens the default device
  # and leaves an Rplots.pdf in the project root. Park it on a throw-away device.
  tmp_dev <- tempfile(fileext = ".pdf")
  pdf(tmp_dev)

  res <- ConsensusClusterPlus(z, maxK = MAXK, reps = REPS, pItem = 0.9, pFeature = 0.8,
                              clusterAlg = "hc", distance = "euclidean",
                              innerLinkage = "ward.D2", finalLinkage = "ward.D2",
                              plot = NULL, verbose = FALSE)

  dev.off()
  unlink(tmp_dev)

  d <- dist(t(z))

  stab <- data.frame(
    platform = label, k = 2:MAXK,
    pac = vapply(res[2:MAXK], function(r) pac(r$consensusMatrix), numeric(1)),
    silhouette = vapply(2:MAXK, function(kk) {
      s <- silhouette(res[[kk]]$consensusClass, d)
      mean(s[, "sil_width"])
    }, numeric(1)),
    stringsAsFactors = FALSE)

  stab$stable <- stab$pac < PAC_STABLE

  list(res = res, stab = stab, dist = d)
}

disc1 <- mn_only("disc_glom_array")
disc2 <- mn_only("cross_glom_rnaseq")
ercb <- mn_only("val_glom_array")

cat("MN samples: array", ncol(disc1$ms), "| RNA-seq", ncol(disc2$ms), "| ERCB", ncol(ercb$ms), "\n")

z1 <- prep(disc1$ms)
z2 <- prep(disc2$ms)

c1 <- run_ccp(z1, "disc_glom_array (n=51)")
c2 <- run_ccp(z2, "cross_glom_rnaseq (n=62)")

stab <- rbind(c1$stab, c2$stab)

write.csv(stab, file.path(tab_dir, "subtype_stability.csv"), row.names = FALSE)

# k = smallest stable partition (PAC < 0.1); if none is stable, the best k is reported but flagged
pick_k <- function(s) {
  st <- s[s$stable, ]
  if (nrow(st)) return(min(st$k))
  s$k[which.min(s$pac)]
}

k1 <- pick_k(c1$stab)
k2 <- pick_k(c2$stab)

stable1 <- any(c1$stab$stable)
stable2 <- any(c2$stab$stable)

cat("chosen k: array", k1, "(stable:", stable1, ") | RNA-seq", k2, "(stable:", stable2, ")\n")

# --- consensus matrices and stability curves ---------------------------------

for (pl in list(list(nm = "array", obj = c1), list(nm = "rnaseq", obj = c2))) {

  for (kk in 2:min(4, MAXK)) {

    cm <- pl$obj$res[[kk]]$consensusMatrix

    save_device_plot(function() pheatmap(cm, cluster_rows = TRUE, cluster_cols = TRUE,
                                         show_rownames = FALSE, show_colnames = FALSE,
                                         color = colorRampPalette(c("white", "#2471A3", "#C0392B"))(101),
                                         main = paste0(pl$nm, " consensus matrix, k = ", kk)),
                     file.path(fig_dir, paste0("06a_", pl$nm, "_consensus_k", kk)),
                     width = 6.5, height = 6)
  }
}

p <- ggplot(stab, aes(k, pac, colour = platform)) +
  geom_line() +
  geom_point(aes(shape = stable), size = 2.6) +
  geom_hline(yintercept = PAC_STABLE, linetype = "dashed", colour = "grey50") +
  labs(title = "Subtype stability (proportion of ambiguous clustering)",
       subtitle = paste0("dashed line = PAC ", PAC_STABLE,
                         "; filled points are stable partitions"),
       x = "k", y = "PAC", colour = NULL, shape = "stable") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

save_plot(p, file.path(fig_dir, "06a_stability_pac"), width = 7, height = 4.6)

# --- subtype-defining genes per platform -------------------------------------

sub_genes <- function(z, cl, label) {

  d <- model.matrix(~ 0 + factor(cl))
  colnames(d) <- paste0("C", levels(factor(cl)))

  fit <- eBayes(lmFit(z, d), trend = TRUE, robust = TRUE)

  out <- list()

  for (i in seq_len(ncol(d))) {
    tt <- topTable(fit, coef = i, number = 200, sort.by = "t")
    out[[i]] <- data.frame(platform = label, cluster = colnames(d)[i], gene = rownames(tt),
                           logFC = tt$logFC, adj.P.Val = tt$adj.P.Val, stringsAsFactors = FALSE)
  }

  do.call(rbind, out)
}

cl1 <- c1$res[[k1]]$consensusClass
cl2 <- c2$res[[k2]]$consensusClass

g1 <- sub_genes(z1, cl1, "array")
g2 <- sub_genes(z2, cl2, "rnaseq")

write.csv(g1, file.path(tab_dir, "subtype_genes_array.csv"), row.names = FALSE)
write.csv(g2, file.path(tab_dir, "subtype_genes_rnaseq.csv"), row.names = FALSE)

for (pl in list(list(nm = "array", z = z1, cl = cl1, g = g1),
                list(nm = "rnaseq", z = z2, cl = cl2, g = g2))) {

  top <- unique(pl$g$gene[pl$g$adj.P.Val < 0.05])[1:min(40, sum(pl$g$adj.P.Val < 0.05))]
  top <- top[!is.na(top)]

  if (length(top) < 3) next

  ann_col <- data.frame(subtype = factor(pl$cl), row.names = colnames(pl$z))

  save_device_plot(function() pheatmap(pl$z[top, , drop = FALSE],
                                       annotation_col = ann_col,
                                       cluster_cols = TRUE, cluster_rows = TRUE,
                                       show_colnames = FALSE, fontsize_row = 7,
                                       color = colorRampPalette(c("#2471A3", "white", "#C0392B"))(101),
                                       main = paste0(pl$nm, " MN subtypes: top ", length(top),
                                                     " defining genes")),
                   file.path(fig_dir, paste0("06a_", pl$nm, "_subtype_heatmap")),
                   width = 7, height = 8)
}

# --- is the array subtype split a quality artefact? --------------------------
# No clinical or demographic covariate is available for these samples, so the only confounder that
# can be checked is technical: per-sample median correlation and mean expression.

qc1 <- read.csv(file.path(OUT, "qc", "03c_disc_glom_array_sample_qc.csv"), stringsAsFactors = FALSE)

qc1 <- qc1[match(disc1$ann$gsm, qc1$gsm), ]

qc1$subtype <- factor(paste0("C", cl1))
qc1$mean_expr <- colMeans(disc1$ms)

qc_tab <- data.frame(
  subtype = names(tapply(qc1$median_corr, qc1$subtype, median)),
  n = as.integer(table(qc1$subtype)),
  median_corr = as.numeric(tapply(qc1$median_corr, qc1$subtype, median)),
  mean_expr = as.numeric(tapply(qc1$mean_expr, qc1$subtype, median)),
  p_median_corr = round(wilcox.test(median_corr ~ subtype, data = qc1)$p.value, 4),
  p_mean_expr = round(wilcox.test(mean_expr ~ subtype, data = qc1)$p.value, 4),
  stringsAsFactors = FALSE)

write.csv(qc_tab, file.path(tab_dir, "subtype_qc_check.csv"), row.names = FALSE)

# --- cross-platform signature concordance (no patient linkage available) -----

sig1 <- unique(g1$gene[g1$adj.P.Val < 0.05])
sig2 <- unique(g2$gene[g2$adj.P.Val < 0.05])

jacc <- if (length(sig1) && length(sig2)) {
  length(intersect(sig1, sig2)) / length(union(sig1, sig2))
} else NA_real_

centroids <- function(z, cl) {
  cls <- sort(unique(cl))
  cen <- sapply(cls, function(c) rowMeans(z[, cl == c, drop = FALSE]))
  colnames(cen) <- paste0("C", cls)
  cen
}

cen1 <- centroids(z1, cl1)
cen2 <- centroids(z2, cl2)

common <- intersect(rownames(cen1), rownames(cen2))
cc <- cor(cen1[common, , drop = FALSE], cen2[common, , drop = FALSE])

conc <- data.frame(
  metric = c("subtype-defining genes: array", "subtype-defining genes: RNA-seq",
             "shared defining genes", "Jaccard of defining genes",
             "max cross-platform centroid correlation",
             "array stable partition", "RNA-seq stable partition"),
  value = c(length(sig1), length(sig2), length(intersect(sig1, sig2)), round(jacc, 4),
            round(max(cc, na.rm = TRUE), 4), stable1, stable2),
  stringsAsFactors = FALSE)

write.csv(conc, file.path(tab_dir, "subtype_signature_concordance.csv"), row.names = FALSE)

cc_df <- as.data.frame(as.table(cc))
names(cc_df) <- c("array_subtype", "rnaseq_subtype", "correlation")

p <- ggplot(cc_df, aes(array_subtype, rnaseq_subtype, fill = correlation)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(correlation, 2)), size = 4) +
  scale_fill_gradient2(low = "#2471A3", mid = "white", high = "#C0392B", midpoint = 0) +
  labs(title = "Cross-platform centroid correlation",
       subtitle = "same genes, different cohorts; labels are platform-specific clusters",
       x = "array subtype centroid", y = "RNA-seq subtype centroid") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

save_plot(p, file.path(fig_dir, "06a_centroid_concordance"), width = 6, height = 5)

# --- transfer to ERCB with both frameworks -----------------------------------

map_ercb <- function(z_ref, cl_ref, ms_ercb) {

  genes <- intersect(rownames(z_ref), rownames(ms_ercb))
  z_ref <- z_ref[genes, , drop = FALSE]
  z_e <- t(scale(t(ms_ercb[genes, , drop = FALSE])))
  z_e[!is.finite(z_e)] <- 0

  cls <- sort(unique(cl_ref))
  cents <- sapply(cls, function(c) rowMeans(z_ref[, cl_ref == c, drop = FALSE]))
  colnames(cents) <- paste0("C", cls)

  cc <- cor(z_e, cents)

  data.frame(sample = colnames(z_e), assigned = colnames(cc)[apply(cc, 1, which.max)],
             max_correlation = round(apply(cc, 1, max), 3), stringsAsFactors = FALSE)
}

e1 <- map_ercb(z1, cl1, ercb$ms)
e2 <- map_ercb(z2, cl2, ercb$ms)

ercb_map <- data.frame(sample = ercb$ann$gsm, patient = ercb$ann$patient_id,
                       array_framework = e1$assigned, array_cor = e1$max_correlation,
                       rnaseq_framework = e2$assigned, rnaseq_cor = e2$max_correlation,
                       stringsAsFactors = FALSE)

write.csv(ercb_map, file.path(tab_dir, "subtype_ercb_assignment.csv"), row.names = FALSE)

p <- ggplot(ercb_map, aes(array_framework, array_cor)) +
  geom_boxplot(outlier.shape = NA, fill = "#D6EAF8") +
  geom_jitter(width = 0.12, size = 1.8, colour = "#2471A3") +
  labs(title = "ERCB MN samples assigned to discovery subtypes",
       subtitle = "nearest centroid on shared genes, array-based framework",
       x = "assigned subtype", y = "centroid correlation") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

save_plot(p, file.path(fig_dir, "06a_ercb_assignment"), width = 5.6, height = 4.6)

# --- report ------------------------------------------------------------------

cat("\n================ STABILITY ================\n")
print(stab, row.names = FALSE)

cat("\narray cluster sizes (k =", k1, "):\n")
print(table(cl1))

cat("\nRNA-seq cluster sizes (k =", k2, "):\n")
print(table(cl2))

cat("\n---- array subtype vs technical covariates ----\n")
print(qc_tab, row.names = FALSE)

cat("\n================ SIGNATURE CONCORDANCE ================\n")
print(conc, row.names = FALSE)

cat("\nERCB assignment with the array framework:\n")
print(table(e1$assigned))

cat("\nERCB assignment with the RNA-seq framework:\n")
print(table(e2$assigned))

cat("\nagreement of the two frameworks on ERCB samples (ARI):",
    round(ari(e1$assigned, e2$assigned), 3), "\n")

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 06a_subtypes.R"),
             capture.output(print(stab, row.names = FALSE)),
             capture.output(print(conc, row.names = FALSE)),
             paste("ERCB framework agreement ARI:", round(ari(e1$assigned, e2$assigned), 3))),
           file.path(PROJ, "logs", "subtypes.log"))

cat("\nDone.\n")

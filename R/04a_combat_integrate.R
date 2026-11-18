#!/usr/bin/env Rscript
# Paper 1 - Phase B step 6: within-platform batch integration of the discovery cohorts.
#
# The plan (steps 6) asks for "merge GPL19983 cohorts with ComBat on log-intensities, study as
# batch, disease as biological covariate; verify batch removal by PCA silhouette".
#
# After the Phase 1 de-duplication only ONE discovery release remains per compartment except for
# the tubulointerstitial arm, where patient de-duplication kept 244 GSE133288 samples plus the 25
# GSE200818 samples that had no GSE133288 counterpart. Both are HuGene-2.1-ST (GPL19983) but they
# were normalised as separate series, so accession is a genuine technical batch. GSE200828 (glom)
# is a single series and needs no correction - that is recorded, not silently skipped.
#
# Outputs
#   work/combat/disc_ti_array_expr_combat.rds         (Entrez-keyed, batch-adjusted)
#   work/combat/disc_ti_array_expr_combat_symbol.rds  (symbol-keyed, batch-adjusted)
#   output/tables/combat_summary.csv
#   output/qc/04a_disc_ti_array_pca_before.{pdf,jpeg}
#   output/qc/04a_disc_ti_array_pca_after.{pdf,jpeg}

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(ggplot2); library(sva) })

set_dir <- file.path(PROJ, "work", "sets")
combat_dir <- file.path(PROJ, "work", "combat")
qc_dir <- file.path(OUT, "qc")
tab_dir <- file.path(OUT, "tables")

dir.create(combat_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, stem, width = 7.5, height = 5.5) {
  if (interactive()) print(p)
  ggsave(paste0(stem, ".pdf"), plot = p, device = "pdf", width = width, height = height)
  ggsave(paste0(stem, ".jpeg"), plot = p, device = "jpeg", quality = 95, dpi = 300,
         bg = "white", width = width, height = height)
}

# --- a group vocabulary shared by every Phase B script -----------------------
# MN / FSGS / MCD are kept apart, both healthy-control classes collapse to Control, and every other
# diagnosis (IgAN, SLE, DN, HT, RPGN, TMD, FSGS_MCD, NEPTUNE "Other") becomes Other.
harmonise <- function(x) {
  out <- ifelse(x %in% c("Control-LD", "Control-TN"), "Control", x)
  out[out %in% c("MN", "FSGS", "MCD", "Control")] <- out[out %in% c("MN", "FSGS", "MCD", "Control")]
  out[!(out %in% c("MN", "FSGS", "MCD", "Control"))] <- "Other"
  out
}

pca_frame <- function(ms, meta, group) {
  v <- apply(ms, 1, sd)
  top <- names(sort(v, decreasing = TRUE))[seq_len(min(2000, length(v)))]
  m <- ms[top, , drop = FALSE]
  m <- m[apply(m, 1, function(z) all(is.finite(z))), , drop = FALSE]
  pc <- prcomp(t(m), scale. = TRUE)
  list(pc = pc, scaled = t(scale(t(m))), top = rownames(m))
}

batch_silhouette <- function(scaled, labels) {
  labels <- factor(labels)
  if (nlevels(labels) < 2) return(NA_real_)
  d <- dist(t(scaled))
  s <- cluster::silhouette(as.integer(labels), d)
  mean(s[, "sil_width"])
}

batch_r2 <- function(pc, labels) {
  y <- pc$x[, seq_len(min(5, ncol(pc$x))), drop = FALSE]
  if (length(unique(labels)) < 2) return(NA_real_)
  r2 <- sapply(seq_len(ncol(y)), function(j) {
    summary(lm(y[, j] ~ factor(labels)))$r.squared
  })
  mean(r2)
}

# --- ComBat on the tubulointerstitial discovery set --------------------------

s <- "disc_ti_array"

expr <- readRDS(file.path(set_dir, paste0(s, "_expr.rds")))
expr_sym <- readRDS(file.path(set_dir, paste0(s, "_expr_symbol.rds")))
ann <- read.csv(file.path(set_dir, paste0(s, "_samples.csv")), stringsAsFactors = FALSE)

stopifnot(identical(colnames(expr), ann$gsm), identical(colnames(expr_sym), ann$gsm))

grp <- harmonise(ann$disease)
batch <- ann$accession

cat("ComBat on", s, "\n")
cat("  samples:", ncol(expr), "| batches:", paste(names(table(batch)), table(batch), sep = "=", collapse = ", "), "\n")
cat("  genes:", nrow(expr), "| NA:", sum(is.na(expr)), "\n")
print(table(batch, grp))

mod <- model.matrix(~ grp)

# guard the classic ComBat failure mode: a group present in one batch only is confounded with it
conf <- table(batch, grp)
cat("  batch x group cells with 0 samples:", sum(conf == 0), "of", length(conf), "\n")

pb <- pca_frame(expr, ann, grp)

pc_before <- pb$pc
sil_before <- batch_silhouette(pb$scaled, batch)
r2_before <- batch_r2(pc_before, batch)

adjusted <- ComBat(dat = expr, batch = batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)
adjusted_sym <- ComBat(dat = expr_sym, batch = batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)

saveRDS(adjusted, file.path(combat_dir, paste0(s, "_expr_combat.rds")))
saveRDS(adjusted_sym, file.path(combat_dir, paste0(s, "_expr_combat_symbol.rds")))

pa <- pca_frame(adjusted, ann, grp)

pc_after <- pa$pc
sil_after <- batch_silhouette(pa$scaled, batch)
r2_after <- batch_r2(pc_after, batch)

# the biological covariate must survive: batch-free disease silhouette should stay well above 0
sil_disease_before <- batch_silhouette(pb$scaled, grp)
sil_disease_after <- batch_silhouette(pa$scaled, grp)

ve_b <- round(100 * summary(pc_before)$importance[2, 1:2], 1)
ve_a <- round(100 * summary(pc_after)$importance[2, 1:2], 1)

df <- rbind(
  data.frame(PC1 = pc_before$x[, 1], PC2 = pc_before$x[, 2], batch = batch, disease = grp,
             stage = "before ComBat", check.names = FALSE),
  data.frame(PC1 = pc_after$x[, 1], PC2 = pc_after$x[, 2], batch = batch, disease = grp,
             stage = "after ComBat", check.names = FALSE))

df$stage <- factor(df$stage, levels = c("before ComBat", "after ComBat"))

p <- ggplot(df, aes(PC1, PC2, colour = batch, shape = disease)) +
  geom_point(size = 2.1, alpha = 0.85) +
  facet_wrap(~ stage, scales = "free") +
  scale_shape_manual(values = c(Control = 1, FSGS = 2, MCD = 3, MN = 16, Other = 4)) +
  labs(title = paste0(s, " - batch structure before and after ComBat"),
       subtitle = paste0("PC1 variance ", ve_b[1], "% -> ", ve_a[1], "%"),
       x = "PC1", y = "PC2") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

save_plot(p, file.path(qc_dir, "04a_disc_ti_array_pca_before_after"), width = 10)

# separate single-panel plots keep the before/after comparison readable at print size
for (st in c("before", "after")) {

  dd <- df[df$stage == paste(st, "ComBat"), ]
  ve <- if (st == "before") ve_b else ve_a

  pp <- ggplot(dd, aes(PC1, PC2, colour = batch, shape = disease)) +
    geom_point(size = 2.2, alpha = 0.85) +
    scale_shape_manual(values = c(Control = 1, FSGS = 2, MCD = 3, MN = 16, Other = 4)) +
    labs(title = paste0(s, " - PCA ", st, " ComBat"),
         x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  save_plot(pp, file.path(qc_dir, paste0("04a_disc_ti_array_pca_", st)))
}

rows <- data.frame(
  set = c("disc_ti_array", "disc_glom_array"),
  batches = c("GSE133288, GSE200818", "GSE200828"),
  samples = c(ncol(expr), 153),
  batch_silhouette_before = c(round(sil_before, 4), NA),
  batch_silhouette_after = c(round(sil_after, 4), NA),
  batch_r2_pc1to5_before = c(round(r2_before, 4), NA),
  batch_r2_pc1to5_after = c(round(r2_after, 4), NA),
  disease_silhouette_before = c(round(sil_disease_before, 4), NA),
  disease_silhouette_after = c(round(sil_disease_after, 4), NA),
  pc1_var_before = c(ve_b[1], NA),
  pc1_var_after = c(ve_a[1], NA),
  action = c("ComBat applied (accession as batch, group as covariate)",
             "no ComBat: single series, no technical batch to remove"),
  matrix_out = c(file.path("work/combat", "disc_ti_array_expr_combat.rds"), "n/a"),
  stringsAsFactors = FALSE)

write.csv(rows, file.path(tab_dir, "combat_summary.csv"), row.names = FALSE)

cat("\n================ COMBAT SUMMARY ================\n")

print(rows, row.names = FALSE)

cat("\nbatch silhouette:", round(sil_before, 3), "->", round(sil_after, 3),
    "| disease silhouette:", round(sil_disease_before, 3), "->", round(sil_disease_after, 3), "\n")

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 04a_combat_integrate.R"),
             capture.output(print(rows, row.names = FALSE))),
           file.path(PROJ, "logs", "combat.log"))

cat("\nDone.\n")

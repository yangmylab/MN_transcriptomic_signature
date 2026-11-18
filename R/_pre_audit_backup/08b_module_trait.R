#!/usr/bin/env Rscript
# Paper 1 - Phase D, step 08b: module eigengenes vs traits.
#
# Traits are MN status, a non-MN glomerulopathy indicator, and the 14 reliable deconvolution
# fractions (Neutrophil and pDC excluded - handoff B8). Clinical outcomes (eGFR/UPCR/remission)
# are still unavailable (B2); when they arrive they are appended to wgcna_traits.csv and this
# script re-runs unchanged.
#
# Two association layers, deliberately kept separate
#   1. PRIMARY  : bicor(module eigengene, trait) with Student p, BH across the whole module x trait
#                 matrix. This is the protocol's module-trait mapping.
#   2. ADJUSTED : limma lmFit(MEs ~ MN + Podocyte + Fibroblast + Myeloid). A module that tracks MN
#                 only because MN samples are more podocyte-rich would lose significance here, so
#                 this layer is the guard against a pure composition artefact.
#   Module significance (MS) = mean gene-level bicor within the module, for each trait.
#
# Outputs
#   output/tables/wgcna_module_trait_<comp>.csv             primary correlations
#   output/tables/wgcna_module_trait_adjusted_<comp>.csv    composition-adjusted MN effect
#   output/tables/wgcna_module_significance_<comp>.csv      gene-level mean correlation
#   output/wgcna/08b_<comp>_{module_trait,ms_bar,me_box}.{pdf,jpeg}

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(WGCNA); library(limma); library(ggplot2) })

tab_dir <- file.path(OUT, "tables")
wk_dir  <- file.path(PROJ, "work", "wgcna")
fig_dir <- file.path(OUT, "wgcna")

TRAIT_CELLS <- c("Podocyte", "PT", "TAL", "DCT", "CD_PC", "PEC", "Endothelium",
                 "Fibroblast", "Mesangium", "Myeloid", "B", "T_NK", "Plasma", "Mast")
ADJ_COVARS  <- c("Podocyte", "Fibroblast", "Myeloid")
FDR_CUT     <- 0.05

fig <- function(stem, width, height, draw) {
  pdf(paste0(stem, ".pdf"), width = width, height = height)
  draw()
  dev.off()
  jpeg(paste0(stem, ".jpeg"), width = width, height = height, units = "in", res = 300,
       quality = 95, bg = "white")
  draw()
  dev.off()
}

write_tab <- function(x, name) {
  f <- file.path(tab_dir, name)
  write.csv(x, f, row.names = FALSE)
  cat("  wrote", f, "(", nrow(x), "rows )\n")
}

strip_me <- function(x) sub("^ME", "", x)

for (comp in c("glom", "ti")) {
  cat("\n================ ", toupper(comp), " ================\n", sep = "")

  obj <- readRDS(file.path(wk_dir, paste0("net_", comp, ".rds")))
  mes <- obj$mes$eigengenes
  tr  <- obj$traits
  dat <- obj$dat
  stopifnot(identical(rownames(mes), tr$gsm))

  # module labels are numeric inside the eigengene names; carry the colour everywhere too
  lab2col <- setNames(obj$module_table$module, obj$module_table$module_label)

  traits_num <- as.matrix(tr[, c("MN", "NonMN_GN", TRAIT_CELLS)])
  cat("modules:", ncol(mes), "| samples:", nrow(mes), "| traits:", ncol(traits_num), "\n")

  # bicor needs a non-zero MAD; traits that are > 50% exact zeros (every tubule-specific
  # fraction in glomerular tissue, and MN itself at 33% prevalence) silently fall back to
  # Pearson. Record which ones did so instead of letting the reader assume bicor throughout.
  tr_mad  <- apply(traits_num, 2, stats::mad)
  tr_zero <- colMeans(traits_num == 0)
  write_tab(data.frame(compartment = comp, trait = names(tr_mad),
                       mad = round(tr_mad, 5), frac_exact_zero = round(tr_zero, 4),
                       method = ifelse(tr_mad > 0, "bicor", "pearson (zero MAD)"),
                       stringsAsFactors = FALSE),
            paste0("wgcna_trait_diagnostics_", comp, ".csv"))

  # --- 1. primary module-trait correlations ---------------------------------
  cp <- bicorAndPvalue(mes, traits_num, use = "p")
  # NB bicorAndPvalue returns $bicor (not $cor); zero-MAD traits silently fall back to Pearson
  long <- data.frame(
    compartment = comp,
    module = strip_me(rownames(cp$bicor)),
    module_color = unname(lab2col[strip_me(rownames(cp$bicor))]),
    trait = rep(strip_me(colnames(cp$bicor)), each = nrow(cp$bicor)),
    cor = round(as.vector(cp$bicor), 4),
    p = as.vector(cp$p),
    stringsAsFactors = FALSE)
  long$padj <- signif(p.adjust(long$p, method = "BH"), 4)
  long$trait_mad <- round(tr_mad[long$trait], 5)
  long$cor_method <- ifelse(long$trait_mad > 0, "bicor", "pearson (zero MAD)")
  long$signif <- ifelse(long$padj < 0.001, "***", ifelse(long$padj < 0.01, "**",
                  ifelse(long$padj < 0.05, "*", "")))
  long <- long[order(long$p), ]
  write_tab(long, paste0("wgcna_module_trait_", comp, ".csv"))

  sig_mn  <- long$module[long$trait == "MN" & long$padj < FDR_CUT]
  sig_pod <- long$module[long$trait == "Podocyte" & long$padj < FDR_CUT]
  cat("modules FDR<", FDR_CUT, " with MN:", if (length(sig_mn)) paste(sig_mn, collapse = ",") else "none",
      "| with Podocyte:", if (length(sig_pod)) paste(sig_pod, collapse = ",") else "none", "\n")

  # --- 2. composition-adjusted MN effect -----------------------------------
  covars <- intersect(ADJ_COVARS, colnames(tr))
  adj_df <- data.frame(MN = tr$MN, tr[, covars, drop = FALSE])
  design <- model.matrix(as.formula(paste("~", paste(c("MN", covars), collapse = " + "))),
                         data = adj_df)
  fit  <- eBayes(lmFit(t(mes), design))
  tt   <- topTable(fit, coef = "MN", number = Inf, sort.by = "P")
  adj  <- data.frame(compartment = comp, module = strip_me(rownames(tt)),
                     module_color = unname(lab2col[strip_me(rownames(tt))]),
                     adj_logFC_MN = round(tt$logFC, 4), t = round(tt$t, 3),
                     p = tt$P.Value, padj = signif(tt$adj.P.Val, 4),
                     covars_adjusted = paste(covars, collapse = "+"), stringsAsFactors = FALSE)
  adj  <- adj[order(adj$p), ]
  write_tab(adj, paste0("wgcna_module_trait_adjusted_", comp, ".csv"))

  # --- 3. module significance (gene-level mean correlation) ----------------
  gene_sig <- bicorAndPvalue(dat, traits_num, use = "p")$bicor   # genes x traits
  stopifnot(identical(rownames(gene_sig), colnames(dat)))
  mod_of_gene <- setNames(obj$module_table$module[match(colnames(dat), obj$module_table$gene)],
                          colnames(dat))
  ms <- do.call(rbind, lapply(setdiff(unique(mod_of_gene), "grey"), function(m) {
    g <- names(mod_of_gene)[mod_of_gene == m]
    data.frame(compartment = comp, module = m, n_genes = length(g),
               module_color = unname(lab2col[m]),
               t(round(colMeans(gene_sig[g, , drop = FALSE]), 4)),
               check.names = FALSE, stringsAsFactors = FALSE)
  }))
  ms <- ms[order(-abs(ms$MN)), ]
  write_tab(ms, paste0("wgcna_module_significance_", comp, ".csv"))

  # --- figures -------------------------------------------------------------
  cor_mat <- cp$bicor
  rownames(cor_mat) <- unname(lab2col[strip_me(rownames(cor_mat))])
  colnames(cor_mat) <- strip_me(colnames(cor_mat))
  if (anyNA(rownames(cor_mat))) rownames(cor_mat) <- strip_me(rownames(cp$bicor))
  # order modules by their MN association so the heatmap reads top-to-bottom
  ord <- order(abs(cor_mat[, "MN"]), decreasing = TRUE)
  cor_mat <- cor_mat[ord, , drop = FALSE]
  star_mat <- matrix("", nrow(cor_mat), ncol(cor_mat), dimnames = dimnames(cor_mat))
  for (i in seq_len(nrow(cor_mat))) for (j in seq_len(ncol(cor_mat))) {
    r <- long[long$module_color == rownames(cor_mat)[i] & long$trait == colnames(cor_mat)[j], ]
    if (nrow(r)) star_mat[i, j] <- r$signif[1]
  }

  fig(file.path(fig_dir, paste0("08b_", comp, "_module_trait")), 9, 6.5, function() {
    par(mar = c(7, 6, 3, 1))
    labeledHeatmap(Matrix = cor_mat,
                   xLabels = colnames(cor_mat), yLabels = rownames(cor_mat),
                   ySymbols = rownames(cor_mat), colorLabels = FALSE,
                   colors = blueWhiteRed(50), textMatrix = star_mat,
                   setStdMargins = FALSE, cex.text = 0.8, zlim = c(-1, 1),
                   main = paste0("Module-trait relationships - ", comp,
                                 " (bicor, BH * p<0.05 ** <0.01 *** <0.001)"))
  })

  fig(file.path(fig_dir, paste0("08b_", comp, "_ms_bar")), 8, 5, function() {
    par(mar = c(4.5, 7, 3, 1))
    ms_o <- ms[order(ms$MN), ]
    cols <- ifelse(ms_o$MN > 0, "firebrick", "steelblue")
    barplot(ms_o$MN, names.arg = ms_o$module_color, horiz = TRUE, las = 1, col = cols,
            border = NA, xlab = "module significance (mean gene bicor with MN status)",
            main = paste0("Module significance for MN status - ", comp))
    abline(v = 0, col = "grey30")
  })

  # eigengene by diagnosis for the modules most associated with MN
  top_mods <- head(strip_me(rownames(cp$bicor))[order(-abs(cp$bicor[, "MN"]))], 4)
  box_df <- do.call(rbind, lapply(top_mods, function(m) {
    data.frame(module = unname(lab2col[m]), group = tr$group, ME = mes[, paste0("ME", m)],
               stringsAsFactors = FALSE)
  }))
  box_df$group <- factor(box_df$group, levels = c("Control", "FSGS", "MCD", "Other", "MN"))

  fig(file.path(fig_dir, paste0("08b_", comp, "_me_box")), 9, 6, function() {
    print(ggplot(box_df, aes(x = group, y = ME, fill = group)) +
            geom_boxplot(outlier.shape = NA, width = 0.62, colour = "grey30") +
            geom_jitter(width = 0.16, size = 0.9, alpha = 0.55, colour = "grey20") +
            facet_wrap(~ module, scales = "free_y", nrow = 2) +
            scale_fill_manual(values = c(Control = "grey75", FSGS = "#8DA0CB", MCD = "#A6D854",
                                         Other = "#FFD92F", MN = "#E41A1C")) +
            labs(title = paste0("Module eigengenes by diagnosis - ", comp),
                 x = NULL, y = "module eigengene") +
            theme_bw(base_size = 11) +
            theme(panel.grid = element_blank(), legend.position = "none",
                  strip.background = element_rect(fill = "grey92", colour = NA),
                  plot.title = element_text(face = "bold")))
  })

  ic <- cp$bicor[, "MN"]
  nm <- unname(lab2col[strip_me(names(ic))])
  cat("MN association ranked:", paste(nm[order(-abs(ic))][1:3], collapse = " > "), "\n")
}

cat("\n08b done\n")

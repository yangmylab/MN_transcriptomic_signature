#!/usr/bin/env Rscript
# Paper 1 - Phase D, step 08a: signed WGCNA networks per compartment.
#
# The separability gate (05a/05b) established that MN can be told apart from other
# glomerulopathies (external AUC 0.87-0.92), so a co-expression network is the right next
# step: it turns the 584/323 thresholded meta-DEGs into modules whose eigengenes can be
# related to traits, and whose preserved hub genes become the Phase E candidate space.
#
# Design (per the protocol, Phase D steps 12-14)
#   network    : signed, bicor, TOMType = signed, dynamic tree cut, mergeCutHeight 0.25
#   gene space : detected in >= 25% of samples (> 5 log2 intensity, the 04b filter), complete,
#                non-constant, then the top 8000 by median absolute deviation (single block, so
#                the TOM stays ~0.5 GB)
#   soft power : smallest power with scale-free fit R2 >= 0.85; if unreachable within the grid
#                the fallback (>= 0.80) is used and the rule is written into the table - never
#                silently relaxed
#   samples    : all retained; the 3 (glom) / 29 (TI) QC-flagged samples are FLAGGED, not dropped,
#                per the frozen QC policy (03d). A sensitivity refit without them is reported.
#   traits     : MN status + group indicators + the 14 reliable deconvolution fractions.
#                Neutrophil and pDC are EXCLUDED - their atlas columns rest on 535 / 97 cells and
#                carry ambient signal (handoff B8), so they are not trustworthy traits.
#   order      : glomerular first (Phase B, the gate and Phase C all put the signal there),
#                then tubulointerstitial with identical parameters.
#
# Outputs
#   output/tables/wgcna_traits.csv                    trait matrix, both compartments
#   output/tables/wgcna_softpower_<comp>.csv          scale-free fit per candidate power
#   output/tables/wgcna_modules_<comp>.csv            gene, module, kME, p.kME
#   output/tables/wgcna_module_sizes_<comp>.csv       module sizes + proportion
#   output/tables/wgcna_outlier_sensitivity.csv       refit without QC-flagged samples
#   output/tables/wgcna_module_stability_<comp>.csv   per-module Jaccard conservation vs the refit
#   output/tables/wgcna_network_summary.csv           one row per compartment (power, modules, ...)
#   work/wgcna/net_<comp>.rds                         blockwiseModules object
#   work/wgcna/me_<comp>.csv, datExpr_<comp>.rds      module eigengenes, filtered expression
#   output/wgcna/08a_<comp>_{sample_dendro,softpower,module_dendro,me_cor}.{pdf,jpeg}

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(WGCNA) })

set_dir  <- file.path(PROJ, "work", "sets")
tab_dir  <- file.path(OUT, "tables")
wk_dir   <- file.path(PROJ, "work", "wgcna")
fig_dir  <- file.path(OUT, "wgcna")

dir.create(wk_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

options(stringsAsFactors = FALSE)

n_threads <- max(1L, parallel::detectCores() - 1L)
allowWGCNAThreads(nThreads = n_threads)

cat("WGCNA", as.character(packageVersion("WGCNA")), "| threads:", n_threads, "\n")

# --- parameters --------------------------------------------------------------

SEED        <- 20261019L
MAX_GENES   <- 8000L
MIN_MODULE  <- 30L
MERGE_CUT   <- 0.25
DEEP_SPLIT  <- 2L
DETECT_FRAC <- 0.25
DETECT_THR  <- 5
R2_CUT      <- 0.85
R2_FALLBACK <- 0.80
POWERS      <- c(1:10, 12, 14, 16, 18, 20)
MAX_POUT    <- 0.05

# 14 reliable atlas cell types (Neutrophil and pDC deliberately excluded, handoff B8)
TRAIT_CELLS <- c("Podocyte", "PT", "TAL", "DCT", "CD_PC", "PEC", "Endothelium",
                 "Fibroblast", "Mesangium", "Myeloid", "B", "T_NK", "Plasma", "Mast")

COMPS <- list(
  glom = list(set = "disc_glom_array",
              expr = file.path(set_dir, "disc_glom_array_expr_symbol.rds"),
              label = "glomerular"),
  ti   = list(set = "disc_ti_array",
              # the only valid TI discovery matrix: GSE133288 + GSE200818 are separate series
              expr = file.path(PROJ, "work", "combat", "disc_ti_array_expr_combat_symbol.rds"),
              label = "tubulointerstitial")
)

# --- helpers -----------------------------------------------------------------

# draw to PDF and JPEG, then close both devices (no grid lines anywhere)
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

# traits: MN status + group + the 14 reliable deconvolution fractions
build_traits <- function(set) {
  ann <- read.csv(file.path(set_dir, paste0(set, "_samples.csv")), stringsAsFactors = FALSE)
  dec <- read.csv(file.path(tab_dir, "deconvolution_fractions.csv"), stringsAsFactors = FALSE)
  dec <- dec[dec$set == set, c("sample", TRAIT_CELLS)]
  # key on (set, sample): cross_glom_rnaseq / cross_ti_rnaseq reuse local sample indices
  stopifnot(!anyDuplicated(dec$sample))
  tr <- merge(ann[, c("gsm", "patient_key", "disease")], dec, by.x = "gsm", by.y = "sample",
              all.x = TRUE, sort = FALSE)
  tr <- tr[match(ann$gsm, tr$gsm), ]
  stopifnot(identical(tr$gsm, ann$gsm))
  stopifnot(!anyNA(tr[, TRAIT_CELLS]))
  grp <- ifelse(tr$disease %in% c("Control-LD", "Control-TN"), "Control", tr$disease)
  grp[!(grp %in% c("MN", "FSGS", "MCD", "Control"))] <- "Other"
  tr$group      <- factor(grp, levels = c("Control", "FSGS", "MCD", "Other", "MN"))
  tr$MN         <- as.integer(tr$group == "MN")
  tr$NonMN_GN   <- as.integer(tr$group %in% c("FSGS", "MCD", "Other"))
  tr
}

# gene filter: complete, detected, non-constant, then top MAX_GENES by MAD
select_genes <- function(ms, max_genes = MAX_GENES) {
  complete <- rowSums(is.na(ms)) == 0
  ms_c     <- ms[complete, , drop = FALSE]
  detected <- rowMeans(ms_c > DETECT_THR) >= DETECT_FRAC
  ms_d     <- ms_c[detected, , drop = FALSE]
  mads     <- if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowMads(ms_d)
  } else {
    apply(ms_d, 1, mad)
  }
  keep <- mads > 0
  ms_v <- ms_d[keep, , drop = FALSE]
  mads <- mads[keep]
  ord  <- order(mads, decreasing = TRUE)
  n    <- min(max_genes, length(ord))
  list(expr  = ms_v[ord[1:n], , drop = FALSE],
       n_complete = sum(complete), n_detected = sum(detected),
       n_variable = length(ord), n_used = n)
}

# smallest power meeting the scale-free criterion, with an explicit fallback rule
pick_power <- function(fi) {
  r2 <- fi$SFT.R.sq
  ok <- !is.na(r2) & r2 >= R2_CUT
  if (any(ok)) {
    p <- fi$Power[which(ok)[1]]
    rule <- sprintf("smallest power with scale-free R2 >= %.2f", R2_CUT)
  } else {
    ok2 <- !is.na(r2) & r2 >= R2_FALLBACK
    if (any(ok2)) {
      p <- fi$Power[which(ok2)[1]]
      rule <- sprintf("FALLBACK: smallest power with R2 >= %.2f (%.2f unreachable)", R2_FALLBACK, R2_CUT)
    } else {
      p <- fi$Power[which.max(r2)]
      rule <- "FALLBACK: power with the best achievable R2"
    }
  }
  list(power = p, rule = rule, r2 = r2[fi$Power == p])
}

# adjusted Rand index without extra dependencies (module agreement, full vs sensitivity fit)
ari <- function(a, b) {
  tab <- table(a, b)
  n   <- sum(tab)
  s   <- sum(choose(tab, 2))
  sa  <- sum(choose(rowSums(tab), 2))
  sb  <- sum(choose(colSums(tab), 2))
  exp <- sa * sb / choose(n, 2)
  mx  <- (sa + sb) / 2
  if (mx == exp) return(1)
  (s - exp) / (mx - exp)
}

fit_network <- function(dat, power, tag) {
  set.seed(SEED)
  blockwiseModules(dat, power = power,
                   networkType = "signed", TOMType = "signed", corType = "bicor",
                   maxPOutliers = MAX_POUT,
                   maxBlockSize = MAX_GENES, minModuleSize = MIN_MODULE,
                   deepSplit = DEEP_SPLIT, mergeCutHeight = MERGE_CUT,
                   numericLabels = TRUE, pamRespectsDendro = FALSE,
                   saveTOMs = FALSE, randomSeed = SEED, nThreads = n_threads, verbose = 2)
}

qc_flagged <- read.csv(file.path(tab_dir, "qc_flagged_samples_retained.csv"), stringsAsFactors = FALSE)

all_traits <- list()
summary_rows <- list()
sens_rows <- list()

# --- per compartment --------------------------------------------------------

for (comp in names(COMPS)) {
  cfg <- COMPS[[comp]]
  cat("\n================ ", toupper(comp), " (", cfg$label, ") ================\n", sep = "")

  ms  <- readRDS(cfg$expr)
  tr  <- build_traits(cfg$set)
  stopifnot(identical(colnames(ms), tr$gsm))

  cat("input:", nrow(ms), "genes x", ncol(ms), "samples\n")
  print(table(tr$group))

  all_traits[[comp]] <- cbind(compartment = comp, set = cfg$set, tr)

  num_tr <- tr[, c("MN", "NonMN_GN", TRAIT_CELLS)]
  rownames(num_tr) <- tr$gsm
  write.csv(num_tr, file.path(wk_dir, paste0("traits_", comp, ".csv")))

  # --- gene space ----------------------------------------------------------
  sel <- select_genes(ms)
  cat(sprintf("genes: complete %d -> detected %d -> variable %d -> used %d (cap %d)\n",
              sel$n_complete, sel$n_detected, sel$n_variable, sel$n_used, MAX_GENES))

  dat <- t(sel$expr)
  gsg <- goodSamplesGenes(dat, verbose = 0)
  cat("goodSamplesGenes:", gsg$allOK, "| bad genes:", sum(!gsg$goodGenes),
      "| bad samples:", sum(!gsg$goodSamples), "\n")
  if (!gsg$allOK) dat <- dat[gsg$goodSamples, gsg$goodGenes]
  saveRDS(dat, file.path(wk_dir, paste0("datExpr_", comp, ".rds")))

  # --- sample dendrogram + trait heatmap -----------------------------------
  sample_tree <- hclust(dist(dat), method = "average")
  tr_plot <- data.frame(MN = tr$MN, Podocyte = tr$Podocyte, Fibroblast = tr$Fibroblast,
                        Myeloid = tr$Myeloid)
  rownames(tr_plot) <- tr$gsm
  tr_plot <- tr_plot[rownames(dat), , drop = FALSE]
  flagged  <- tr$gsm %in% qc_flagged$gsm[qc_flagged$set == cfg$set]
  tr_cols  <- numbers2colors(as.matrix(tr_plot), signed = FALSE)
  tr_cols  <- cbind(tr_cols, QC_flagged = ifelse(flagged, "red", "grey80"))

  fig(file.path(fig_dir, paste0("08a_", comp, "_sample_dendro")), 9, 6, function() {
    par(mar = c(1, 5, 3, 1))
    plotDendroAndColors(sample_tree, tr_cols,
                        groupLabels = c(colnames(tr_plot), "QC_flagged"),
                        dendroLabels = FALSE, hang = 0.03, addGuide = FALSE, guideHang = 0.05,
                        main = paste0("Sample clustering - ", cfg$label, " (", cfg$set, ")"),
                        marAll = c(1, 6, 3, 1))
  })

  # --- soft threshold ------------------------------------------------------
  sft <- pickSoftThreshold(dat, powerVector = POWERS, networkType = "signed",
                           RsquaredCut = R2_CUT, verbose = 0)
  fi <- sft$fitIndices
  fi_out <- data.frame(power = fi$Power, scale_free_R2 = fi$SFT.R.sq,
                       slope = fi$slope, mean_k = fi$mean.k., median_k = fi$median.k.,
                       max_k = fi$max.k.)
  write_tab(fi_out, paste0("wgcna_softpower_", comp, ".csv"))

  pw <- pick_power(fi)
  cat(sprintf("soft power: %d (R2 = %.3f) - %s\n", pw$power, pw$r2, pw$rule))

  fig(file.path(fig_dir, paste0("08a_", comp, "_softpower")), 8, 4.2, function() {
    par(mfrow = c(1, 2), mar = c(4.2, 4.5, 2.5, 1))
    plot(fi$Power, fi$SFT.R.sq, type = "n", xlab = "soft threshold (power)",
         ylab = expression("scale-free fit, " * R^2), main = "Scale independence")
    abline(h = R2_CUT, lty = 2, col = "grey40")
    abline(h = R2_FALLBACK, lty = 3, col = "grey60")
    points(fi$Power, fi$SFT.R.sq, pch = 19, col = "steelblue")
    points(pw$power, pw$r2, pch = 19, col = "firebrick", cex = 1.6)
    text(pw$power, pw$r2, labels = pw$power, pos = 3, col = "firebrick")
    plot(fi$Power, fi$mean.k., type = "n", xlab = "soft threshold (power)",
         ylab = "mean connectivity", main = "Mean connectivity")
    points(fi$Power, fi$mean.k., pch = 19, col = "steelblue")
    points(pw$power, fi$mean.k.[fi$Power == pw$power], pch = 19, col = "firebrick", cex = 1.6)
    text(pw$power, fi$mean.k.[fi$Power == pw$power], labels = pw$power, pos = 3, col = "firebrick")
  })

  # --- network -------------------------------------------------------------
  net <- fit_network(dat, pw$power, comp)
  colors <- net$colors
  mod_levels <- sort(unique(colors))
  mod_colors <- labels2colors(colors)

  cat("modules (incl. grey):", length(mod_levels), "| genes:", length(colors), "\n")
  print(table(mod_colors))

  mes <- moduleEigengenes(dat, colors = colors, excludeGrey = TRUE)
  kme <- signedKME(dat, mes$eigengenes, corFnc = "bicor")
  # With numericLabels = TRUE the eigengenes are ME1..MEk and signedKME returns "kME_ME1" etc.
  # Normalise the columns to kME<label> so the per-gene lookup below matches on the numeric label
  # (looking them up by colour name silently produced an all-NA kME column).
  colnames(kme) <- paste0("kME", sub("^ME", "", colnames(mes$eigengenes)))

  own_kme <- rep(NA_real_, length(colors))
  own_p   <- rep(NA_real_, length(colors))
  for (i in seq_along(colors)) {
    if (colors[i] == 0) next
    col <- paste0("kME", colors[i])
    if (!col %in% colnames(kme)) next
    own_kme[i] <- kme[i, col]
    own_p[i]   <- corPvalueStudent(own_kme[i], nrow(dat))
  }

  stopifnot(sum(is.na(own_kme[colors != 0])) == 0)

  mod_tab <- data.frame(gene = colnames(dat), module = mod_colors, module_label = colors,
                        kME = own_kme, p.kME = own_p, stringsAsFactors = FALSE)
  write_tab(mod_tab, paste0("wgcna_modules_", comp, ".csv"))

  sizes <- as.data.frame(table(module = mod_colors), stringsAsFactors = FALSE)
  names(sizes)[2] <- "n_genes"
  sizes$proportion <- round(sizes$n_genes / sum(sizes$n_genes), 4)
  sizes <- sizes[order(-sizes$n_genes), ]
  write_tab(sizes, paste0("wgcna_module_sizes_", comp, ".csv"))

  me_out <- data.frame(sample = rownames(mes$eigengenes), mes$eigengenes, check.names = FALSE)
  write.csv(me_out, file.path(wk_dir, paste0("me_", comp, ".csv")), row.names = FALSE)

  saveRDS(list(net = net, mes = mes, dat = dat, traits = tr, power = pw$power,
               softpower = fi_out, module_table = mod_tab, genes_used = colnames(dat)),
          file.path(wk_dir, paste0("net_", comp, ".rds")))

  fig(file.path(fig_dir, paste0("08a_", comp, "_module_dendro")), 10, 6, function() {
    par(mar = c(1, 5, 3, 1))
    plotDendroAndColors(net$dendrograms[[1]], mod_colors[net$blockGenes[[1]]],
                        groupLabels = "module", dendroLabels = FALSE, hang = 0.03,
                        addGuide = FALSE, guideHang = 0.05,
                        main = paste0("Gene dendrogram and modules - ", cfg$label,
                                      " (power ", pw$power, ")"),
                        marAll = c(1, 6, 3, 1))
  })

  # --- module eigengene correlation heatmap --------------------------------
  me_cor <- cor(mes$eigengenes, use = "p")
  fig(file.path(fig_dir, paste0("08a_", comp, "_me_cor")), 8, 7, function() {
    par(mar = c(6, 6, 3, 1))
    labeledHeatmap(Matrix = me_cor,
                   xLabels = sub("^ME", "", colnames(me_cor)),
                   yLabels = sub("^ME", "", rownames(me_cor)),
                   ySymbols = sub("^ME", "", rownames(me_cor)),
                   colorLabels = FALSE, colors = blueWhiteRed(50),
                   textMatrix = round(me_cor, 2), setStdMargins = FALSE,
                   cex.text = 0.5, zlim = c(-1, 1),
                   main = paste0("Module eigengene correlation - ", cfg$label))
  })

  # --- sensitivity: refit without QC-flagged samples -----------------------
  keep_s <- !(rownames(dat) %in% qc_flagged$gsm[qc_flagged$set == cfg$set])
  if (sum(!keep_s) > 0) {
    net_s <- fit_network(dat[keep_s, ], pw$power, paste0(comp, "_sens"))
    fmod <- setNames(mod_colors, colnames(dat))
    smod <- setNames(labels2colors(net_s$colors), colnames(dat)[keep_s])

    # gene-partition ARI over every gene, and over genes that are non-grey in BOTH fits
    # (the all-gene value is dominated by the grey class and is not interpretable alone)
    g_all  <- intersect(names(fmod), names(smod))
    ari_all <- ari(fmod[g_all], smod[g_all])
    g_ng   <- g_all[fmod[g_all] != "grey" & smod[g_all] != "grey"]
    ari_ng <- if (length(g_ng)) ari(fmod[g_ng], smod[g_ng]) else NA_real_

    # per-module conservation: best Jaccard against any module of the sensitivity fit
    f_mods <- setdiff(unique(fmod), "grey")
    s_mods <- setdiff(unique(smod), "grey")
    stab <- do.call(rbind, lapply(f_mods, function(m) {
      A <- names(fmod)[fmod == m]
      js <- vapply(s_mods, function(n) {
        B <- names(smod)[smod == n]
        length(intersect(A, B)) / length(union(A, B))
      }, numeric(1))
      best <- s_mods[which.max(js)]
      B <- names(smod)[smod == best]
      data.frame(module = m, size_full = length(A), best_match = best,
                 size_match = length(B), jaccard = round(max(js), 4),
                 frac_recovered = round(length(intersect(A, B)) / length(A), 4),
                 stringsAsFactors = FALSE)
    }))
    stab <- stab[order(-stab$size_full), ]
    write_tab(stab, paste0("wgcna_module_stability_", comp, ".csv"))

    cat(sprintf("sensitivity (n = %d, %d flagged removed): %d modules | ARI all %.3f | ARI non-grey %.3f | median best Jaccard %.3f\n",
                sum(keep_s), sum(!keep_s), length(unique(net_s$colors)), ari_all, ari_ng,
                median(stab$jaccard)))

    sens_rows[[comp]] <- data.frame(
      compartment = comp, set = cfg$set, n_full = nrow(dat), n_sensitivity = sum(keep_s),
      n_flagged_removed = sum(!keep_s), modules_full = length(mod_levels),
      modules_sensitivity = length(unique(net_s$colors)),
      gene_partition_ARI = round(ari_all, 4),
      gene_partition_ARI_nongrey = round(ari_ng, 4),
      median_module_jaccard = round(median(stab$jaccard), 4),
      n_modules_jaccard_gt_0.5 = sum(stab$jaccard > 0.5),
      power = pw$power, stringsAsFactors = FALSE)
    rm(net_s, stab)
  }

  summary_rows[[comp]] <- data.frame(
    compartment = comp, set = cfg$set, label = cfg$label,
    n_samples = nrow(dat), n_genes_complete = sel$n_complete, n_genes_detected = sel$n_detected,
    n_genes_variable = sel$n_variable, n_genes_used = sel$n_used,
    soft_power = pw$power, scale_free_R2 = round(pw$r2, 4), power_rule = pw$rule,
    n_modules = length(mod_levels), n_modules_excl_grey = length(mod_levels) - 1,
    n_assigned = sum(colors != 0), prop_assigned = round(mean(colors != 0), 4),
    largest_module = sizes$n_genes[1], min_module = MIN_MODULE, mergeCutHeight = MERGE_CUT,
    deepSplit = DEEP_SPLIT, corType = "bicor", networkType = "signed",
    maxBlockSize = MAX_GENES, seed = SEED, stringsAsFactors = FALSE)

  rm(net, dat, mes, kme)
  gc()
}

# --- write the cross-compartment tables -------------------------------------

traits_all <- do.call(rbind, all_traits)
write_tab(traits_all, "wgcna_traits.csv")
write_tab(do.call(rbind, summary_rows), "wgcna_network_summary.csv")
if (length(sens_rows)) write_tab(do.call(rbind, sens_rows), "wgcna_outlier_sensitivity.csv")

cat("\n--- network summary ---\n")
print(do.call(rbind, summary_rows)[, c("compartment", "n_samples", "n_genes_used", "soft_power",
                                       "scale_free_R2", "n_modules_excl_grey", "prop_assigned")])
cat("\n08a done\n")

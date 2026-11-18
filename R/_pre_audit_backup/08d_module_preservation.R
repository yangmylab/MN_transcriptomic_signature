#!/usr/bin/env Rscript
# Paper 1 - Phase D, step 08d: module preservation in the independent ERCB cohorts.
#
# Protocol step 14: only modules preserved in the validation cohorts proceed to the hub-gene
# stage. Preservation is tested with WGCNA::modulePreservation, which fixes our discovery module
# assignment as the reference and asks whether the same density/connectivity structure can be
# recovered in the test cohort better than random gene sets of equal size (permutation Z scores).
#
# Tests
#   primary   : disc_glom_array -> val_glom_array     (ERCB glomerular, 0% patient overlap)
#               disc_ti_array   -> val_ti_array       (ERCB tubulointerstitial)
#   secondary : disc_glom_array -> cross_glom_rnaseq  (different platform AND different cohort -
#               expected to be weaker, reported for completeness, never used to qualify a module)
#
# API notes for WGCNA 1.74 (these cost several failed runs, see handoff P-35/P-36):
#   * referenceNetworks is a NUMERIC index into multiData - a precomputed network object is rejected
#   * our own modules are supplied through multiColor = list(<colour vector>, NA)
#   * multiData is a 2-level list: list(<set> = list(data = samples x genes))
#   * the returned statistics are the "qual" family (Zsummary.qual, Z.meanAdj.qual, ...), NOT the
#     legacy Zsummary.pres; q-values come from calculateQvalue = TRUE
#   * the statistics carry no module labels, so entries are mapped back to colours by matching the
#     module sizes recomputed on the shared gene set (the gold module is inserted in the array)
#
# Preservation is scored on Zsummary.qual > 10 (effect size, the protocol's threshold) AND
# qsummary.qual < 0.05 (significance), because the .qual Z scores run much larger than the legacy
# Zsummary.pres and a bare threshold would call almost everything preserved.
#
# Outputs
#   output/tables/wgcna_preservation_<comp>.csv        per module, primary test
#   output/tables/wgcna_preservation_secondary.csv     cross-platform test
#   output/tables/wgcna_preservation_summary.csv       per test: preserved / weak / none
#   work/wgcna/preservation_<comp>.rds                 raw modulePreservation objects
#   output/wgcna/08d_<name>_preservation.{pdf,jpeg}

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(WGCNA) })

set_dir <- file.path(PROJ, "work", "sets")
tab_dir <- file.path(OUT, "tables")
wk_dir  <- file.path(PROJ, "work", "wgcna")
fig_dir <- file.path(OUT, "wgcna")

n_threads <- max(1L, parallel::detectCores() - 1L)
allowWGCNAThreads(nThreads = n_threads)

SEED <- 20261019L
N_PERM_PRIMARY   <- 100L
N_PERM_SECONDARY <- 50L
MAX_MODULE_SIZE  <- 5000L
MAX_GOLD_SIZE    <- 1000L
Z_STRONG <- 10
Z_WEAK   <- 2
Q_CUT    <- 0.05

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

verdict_of <- function(z, q) {
  ifelse(is.na(z) | is.na(q), "not testable",
  ifelse(z > Z_STRONG & q < Q_CUT, "preserved",
  ifelse(z >= Z_WEAK & q < Q_CUT, "weak", "not preserved")))
}

# map the unlabelled statistic vector back to module colours via the shared-gene sizes
label_by_size <- function(sizes_obs, sizes_shared) {
  pool <- c(as.integer(sizes_shared), MAX_GOLD_SIZE)
  names(pool) <- c(names(sizes_shared), "gold")
  if (anyDuplicated(pool)) warning("duplicate module sizes: the label mapping is ambiguous")
  lab <- names(pool)[match(sizes_obs, pool)]
  if (anyNA(lab)) {
    cat("  observed sizes:", paste(sizes_obs, collapse = ","), "\n")
    cat("  shared sizes  :", paste(pool, collapse = ","), "\n")
    stop("cannot map preservation entries to modules")
  }
  lab
}

run_preservation <- function(comp, test_set, n_perm, label) {
  ref <- readRDS(file.path(wk_dir, paste0("net_", comp, ".rds")))
  dat <- ref$dat                       # samples x genes, exactly the network input
  cols <- setNames(ref$module_table$module[match(colnames(dat), ref$module_table$gene)],
                   colnames(dat))

  ms <- readRDS(file.path(set_dir, paste0(test_set, "_expr_symbol.rds")))
  shared <- intersect(colnames(dat), rownames(ms))
  ma <- ms[shared, , drop = FALSE]
  ma <- ma[rowSums(is.na(ma)) == 0, , drop = FALSE]
  keep_g <- intersect(colnames(dat), rownames(ma))
  dat_test <- t(ma[keep_g, , drop = FALSE])
  stopifnot(all(is.finite(dat_test)))

  cat(sprintf("\n[%s] %s (%d genes x %d samples) -> %s (%d genes x %d samples): %d shared genes, %d permutations\n",
              label, comp, ncol(dat), nrow(dat), test_set, nrow(ms), ncol(ms), length(keep_g), n_perm))

  sizes_shared <- table(cols[keep_g])

  t0 <- Sys.time()
  p <- modulePreservation(list(discovery = list(data = dat), test = list(data = dat_test)),
                          multiColor = list(cols, NA),
                          referenceNetworks = 1, testNetworks = 2,
                          networkType = "signed", corFnc = "bicor",
                          nPermutations = n_perm, randomSeed = SEED,
                          maxModuleSize = MAX_MODULE_SIZE, maxGoldModuleSize = MAX_GOLD_SIZE,
                          calculateQvalue = TRUE)
  cat("  elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "min\n")

  r  <- p$quality$Z[[1]][[2]]
  qq <- p$quality$q[[1]][[2]]
  o  <- p$quality$observed[[1]][[2]]

  labs <- label_by_size(r$moduleSize, sizes_shared)

  out <- data.frame(
    compartment = comp, test_set = test_set, comparison = label,
    module = labs, module_size_shared = as.integer(r$moduleSize),
    Zsummary = round(r$Zsummary.qual, 3),
    Z_propVarExplained = round(r$Z.propVarExplained.qual, 3),
    Z_meanSignAwareKME = round(r$Z.meanSignAwareKME.qual, 3),
    Z_meanSignAwareCorDat = round(r$Z.meanSignAwareCorDat.qual, 3),
    Z_meanAdj = round(r$Z.meanAdj.qual, 3),
    Z_meanMAR = round(r$Z.meanMAR.qual, 3),
    medianRank = round(o$medianRank.qual, 2),
    qsummary = signif(qq$qsummary.qual, 4),
    n_permutations = n_perm,
    stringsAsFactors = FALSE)
  out$verdict <- verdict_of(out$Zsummary, out$qsummary)
  out <- out[order(-out$Zsummary), ]

  cat("  preserved:", sum(out$verdict == "preserved"), "| weak:", sum(out$verdict == "weak"),
      "| not preserved:", sum(out$verdict == "not preserved"), "\n")
  print(out[, c("module", "module_size_shared", "Zsummary", "qsummary", "verdict")], row.names = FALSE)

  list(table = out, pres = p, shared = length(keep_g))
}

res_glom <- run_preservation("glom", "val_glom_array", N_PERM_PRIMARY, "primary")
write_tab(res_glom$table, "wgcna_preservation_glom.csv")
saveRDS(res_glom$pres, file.path(wk_dir, "preservation_glom.rds"))

res_ti <- run_preservation("ti", "val_ti_array", N_PERM_PRIMARY, "primary")
write_tab(res_ti$table, "wgcna_preservation_ti.csv")
saveRDS(res_ti$pres, file.path(wk_dir, "preservation_ti.rds"))

res_x <- run_preservation("glom", "cross_glom_rnaseq", N_PERM_SECONDARY, "secondary cross-platform")
write_tab(res_x$table, "wgcna_preservation_secondary.csv")

summ <- do.call(rbind, lapply(list(res_glom, res_ti, res_x), function(r) {
  t <- r$table
  data.frame(compartment = t$compartment[1], test_set = t$test_set[1],
             comparison = t$comparison[1], n_modules = nrow(t),
             shared_genes = r$shared,
             n_preserved = sum(t$verdict == "preserved"),
             n_weak = sum(t$verdict == "weak"),
             n_not_preserved = sum(t$verdict == "not preserved"),
             median_Zsummary = round(median(t$Zsummary, na.rm = TRUE), 3),
             max_Zsummary = round(max(t$Zsummary, na.rm = TRUE), 3),
             n_permutations = t$n_permutations[1], stringsAsFactors = FALSE)
}))
write_tab(summ, "wgcna_preservation_summary.csv")

for (nm in c("glom", "ti", "secondary")) {
  r <- switch(nm, glom = res_glom, ti = res_ti, secondary = res_x)
  t <- r$table
  t <- t[t$module != "gold", ]
  fig(file.path(fig_dir, paste0("08d_", nm, "_preservation")), 7.5, 5, function() {
    par(mar = c(4.5, 4.5, 3, 1))
    yl <- range(c(-2, t$Zsummary), na.rm = TRUE)
    plot(t$module_size_shared, t$Zsummary, log = "x", pch = 19, cex = 1.2,
         col = ifelse(t$verdict == "preserved", "seagreen",
               ifelse(t$verdict == "weak", "orange", "grey50")),
         xlab = "module size (genes measurable in both cohorts, log scale)",
         ylab = expression(Z[summary] ~ "(.qual)"),
         ylim = c(yl[1], max(12, yl[2] * 1.08)),
         main = paste0("Module preservation - ", t$compartment[1], " -> ", t$test_set[1]))
    abline(h = c(Z_WEAK, 5, Z_STRONG), lty = c(3, 2, 1), col = "grey40")
    text(t$module_size_shared, t$Zsummary, labels = t$module, pos = 4, cex = 0.75)
    legend("bottomright", legend = c("Z>10 & q<0.05 preserved", "2-10 weak", "<2 none"),
           pch = 19, col = c("seagreen", "orange", "grey50"), bty = "n", cex = 0.85)
  })
}

cat("\n--- preservation summary ---\n")
print(summ[, c("compartment", "test_set", "n_modules", "n_preserved", "n_weak",
               "n_not_preserved", "median_Zsummary")])
cat("\n08d done\n")

#!/usr/bin/env Rscript
# Track B - cache the harmonised matrices the panel benchmark needs (2026-11-18)
#
# The benchmark (10c) and the gate (05a/05b) must use ONE harmonisation implementation, otherwise we
# repeat P-63 (two private copies of the same statistical rule). This script builds and caches:
#
#   z_harmonised  - validation cohorts quantile-mapped onto the discovery reference, then z-scored
#                   with DISCOVERY means/sds (the valid cross-platform transfer, F62)
#   z_transductive - each cohort z-scored against ITS OWN mean/sd (the field's implicit default when
#                   people apply scale() per cohort; kept only so the benchmark can quantify it)
#
# Inputs   work/sets/{disc_glom_array,val_glom_array,cross_glom_rnaseq}_{expr_symbol.rds,samples.csv}
# Outputs  work/bench/10b_harmonised.rds          list of matrices, labels, patient ids
#          output/tables/10b_harmonisation_diagnostics.csv
# parallel: serial - one pass over three matrices, fork overhead would exceed the work

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(limma) })

set_dir <- file.path(PROJ, "work", "sets")
bench_dir <- file.path(PROJ, "work", "bench")
tab_dir <- file.path(OUT, "tables")

dir.create(bench_dir, showWarnings = FALSE, recursive = TRUE)

# --- load --------------------------------------------------------------------

load_set <- function(set) {

  ms <- readRDS(file.path(set_dir, paste0(set, "_expr_symbol.rds")))
  ann <- read.csv(file.path(set_dir, paste0(set, "_samples.csv")), stringsAsFactors = FALSE)
  stopifnot(identical(colnames(ms), ann$gsm))

  grp <- ifelse(ann$disease %in% c("Control-LD", "Control-TN"), "Control", ann$disease)
  keep <- ann$disease != "UNKNOWN" & grp != "Control"
  raw <- grp[keep]

  list(ms = ms[, keep, drop = FALSE],
       y = factor(ifelse(raw == "MN", "MN", "rest"), levels = c("rest", "MN")),
       group = raw,
       patient = ann$patient_id[keep],
       gsm = ann$gsm[keep])

}

disc <- load_set("disc_glom_array")
ercb <- load_set("val_glom_array")
rna  <- load_set("cross_glom_rnaseq")

genes <- Reduce(intersect, list(rownames(disc$ms), rownames(ercb$ms), rownames(rna$ms)))

cat("discovery:", ncol(disc$ms), "samples | ERCB:", ncol(ercb$ms), "| RNA-seq:", ncol(rna$ms), "\n")
cat("shared genes:", length(genes), "\n")

# --- harmonisation (identical transform to 05a/05b) --------------------------

std_fit <- function(M) {

  mu <- rowMeans(M)
  sdv <- apply(M, 1, sd)
  sdv[!is.finite(sdv) | sdv == 0] <- 1

  list(mu = mu, sd = sdv)

}

std <- function(M, sc) {

  Z <- sweep(sweep(M, 1, sc$mu), 1, sc$sd, "/")
  Z[!is.finite(Z)] <- 0

  Z

}

quantile_map <- function(M, ref) {

  out <- M

  for (i in seq_len(nrow(M))) {

    v <- M[i, ]
    ok <- is.finite(v)

    if (sum(ok) < 3) { out[i, ] <- NA_real_; next }

    r <- rank(v[ok], ties.method = "average")
    idx <- pmin(ncol(ref), pmax(1L, round(r / (sum(ok) + 1) * ncol(ref))))

    out[i, ] <- NA_real_
    out[i, ok] <- ref[i, idx]

  }

  out

}

M_disc <- disc$ms[genes, , drop = FALSE]
M_ercb <- ercb$ms[genes, , drop = FALSE]
M_rna  <- rna$ms[genes, , drop = FALSE]

ref_sorted <- t(apply(M_disc, 1, sort))
sc_disc <- std_fit(M_disc)

z_harmonised <- list(
  disc = std(M_disc, sc_disc),
  ercb = std(quantile_map(M_ercb, ref_sorted), sc_disc),
  rna  = std(quantile_map(M_rna,  ref_sorted), sc_disc))

# the transductive alternative: each cohort scaled against its own distribution (only for R4)
z_transductive <- list(
  disc = std(M_disc, std_fit(M_disc)),
  ercb = std(M_ercb, std_fit(M_ercb)),
  rna  = std(M_rna,  std_fit(M_rna)))

# --- diagnostics: the |z| inflation that F62 diagnosed -----------------------

diag_rows <- lapply(c("disc", "ercb", "rna"), function(k) {

  zh <- z_harmonised[[k]]
  zt <- z_transductive[[k]]

  data.frame(set = k, genes = nrow(zh), samples = ncol(zh),
             median_absz_harmonised = round(median(abs(zh)), 3),
             p999_absz_harmonised = round(quantile(abs(zh), 0.999), 2),
             median_absz_transductive = round(median(abs(zt)), 3),
             stringsAsFactors = FALSE)

})

diag <- do.call(rbind, diag_rows)

write.csv(diag, file.path(tab_dir, "10b_harmonisation_diagnostics.csv"), row.names = FALSE)

cat("\n--- |z| diagnostics (harmonised should sit near 1) ---\n")
print(diag, row.names = FALSE)

saveRDS(list(genes = genes,
             z_harmonised = z_harmonised,
             z_transductive = z_transductive,
             y = list(disc = disc$y, ercb = ercb$y, rna = rna$y),
             group = list(disc = disc$group, ercb = ercb$group, rna = rna$group),
             patient = list(disc = disc$patient, ercb = ercb$patient, rna = rna$patient),
             gsm = list(disc = disc$gsm, ercb = ercb$gsm, rna = rna$gsm)),
        file.path(bench_dir, "10b_harmonised.rds"))

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 10b_harmonised_matrices.R"),
             capture.output(print(diag, row.names = FALSE))),
           file.path(PROJ, "logs", "10b_harmonised_matrices.log"))

cat("\ncached:", file.path(bench_dir, "10b_harmonised.rds"), "\n")
cat("Done.\n")

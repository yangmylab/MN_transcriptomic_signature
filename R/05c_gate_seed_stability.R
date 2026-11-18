#!/usr/bin/env Rscript
# Paper 1 - is the harmonised separability gate stable across random-number streams?
#
# The audit of 2026-11-18 (F62) found that the *unharmonised* gate was unstable: at K = all the ERCB
# AUC ranged 0.638-0.679 across twelve folds draws, with per-diagnosis AUCs inverting. The fix
# (quantile-mapping each validation cohort onto the discovery reference before applying training
# statistics) was then validated once. This script closes that loop: it refits the same model under
# several independent RNG streams and reports the spread, so "ERCB AUC 0.861-0.922 across feature
# sizes" is a range over feature counts AND over fold partitions, not a single lucky draw.
#
# Nothing here selects a model. The reported gate stays 05a's; this measures how much of that number
# is partition noise.
#
# Inputs   work/sets/{disc_glom_array,val_glom_array}_{expr_symbol.rds,samples.csv}
# Outputs  output/tables/05c_gate_seed_stability.csv
#          output/models_gate/05c_gate_seed_stability.{pdf,jpeg}
# parallel: serial by design - 10 fits of one elastic net, fork overhead would exceed the work

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(limma); library(glmnet); library(pROC); library(ggplot2) })

set_dir <- file.path(PROJ, "work", "sets")
tab_dir <- file.path(OUT, "tables")
fig_dir <- file.path(OUT, "models_gate")

SEED0 <- 20261018L
N_DRAW <- 10L
K_GRID <- c(1000L, NA_integer_)          # NA = all genes; both are reported in 05a
ALPHA <- 0.5
K_LAB <- function(k) ifelse(is.na(k), "all", as.character(k))

mg <- function(set) {

  ms <- readRDS(file.path(set_dir, paste0(set, "_expr_symbol.rds")))
  ann <- read.csv(file.path(set_dir, paste0(set, "_samples.csv")), stringsAsFactors = FALSE)
  grp <- ifelse(ann$disease %in% c("Control-LD", "Control-TN"), "Control", ann$disease)
  grp[!(grp %in% c("MN", "FSGS", "MCD", "Control"))] <- "Other"
  keep <- ann$disease != "UNKNOWN" & grp != "Control"

  list(ms = ms[, keep, drop = FALSE], grp = grp[keep])

}

# training statistics only, after harmonising the validation cohort onto the discovery reference
# (identical transform to 05a/05b - see F62)
std_fit <- function(tr_ms) {

  mu <- rowMeans(tr_ms)
  sdv <- apply(tr_ms, 1, sd)
  sdv[!is.finite(sdv) | sdv == 0] <- 1

  list(mu = mu, sd = sdv)

}

std <- function(ms, sc) {

  z <- sweep(sweep(ms, 1, sc$mu), 1, sc$sd, "/")
  z[!is.finite(z)] <- 0

  z

}

quantile_map <- function(ms, ref) {

  out <- ms

  for (i in seq_len(nrow(ms))) {

    v <- ms[i, ]
    ok <- is.finite(v)

    if (sum(ok) < 3) { out[i, ] <- NA_real_; next }

    r <- rank(v[ok], ties.method = "average")
    idx <- pmin(ncol(ref), pmax(1L, round(r / (sum(ok) + 1) * ncol(ref))))

    out[i, ] <- NA_real_
    out[i, ok] <- ref[i, idx]

  }

  out

}

tr <- mg("disc_glom_array")
ex <- mg("val_glom_array")

genes <- intersect(rownames(tr$ms), rownames(ex$ms))

M_tr <- tr$ms[genes, , drop = FALSE]
M_ex <- quantile_map(ex$ms[genes, , drop = FALSE], t(apply(M_tr, 1, sort)))

sc <- std_fit(M_tr)
X_tr <- std(M_tr, sc)
X_ex <- std(M_ex, sc)

ybin <- factor(ifelse(tr$grp == "MN", "MN", "rest"), levels = c("rest", "MN"))
# P-22: roc() needs a two-level response; passing the multi-class group vector fails with
# "No control observation".
y_ex <- factor(ifelse(ex$grp == "MN", "MN", "rest"), levels = c("rest", "MN"))

cat("training:", ncol(X_tr), "samples |", nrow(X_tr), "genes\n")
cat("external (ERCB):", ncol(X_ex), "samples\n\n")
cat("draws:", N_DRAW, "| K:", paste(K_LAB(K_GRID), collapse = ", "), "| alpha:", ALPHA, "\n\n")

rows <- list()

for (k in K_GRID) {

  # gene ranking is deterministic given the training data, so it is fixed across draws; only the
  # cv.glmnet fold partition (and therefore lambda) varies. That is exactly the axis under test.
  d <- model.matrix(~ 0 + factor(tr$grp))
  colnames(d) <- levels(factor(tr$grp))
  fit <- eBayes(lmFit(X_tr, d), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = which(colnames(d) == "MN"), number = Inf, sort.by = "t")
  ranked <- rownames(tt)[tt$logFC != 0]
  g <- if (is.na(k)) ranked else head(ranked, k)

  Xg_tr <- t(X_tr[g, , drop = FALSE])
  Xg_ex <- t(X_ex[g, , drop = FALSE])

  for (i in seq_len(N_DRAW)) {

    seed <- SEED0 + i
    set.seed(seed)

    cv <- cv.glmnet(Xg_tr, ybin, family = "binomial", alpha = ALPHA,
                    nfolds = 5, type.measure = "auc")
    pr <- as.numeric(predict(cv$glmnet.fit, newx = Xg_ex, s = cv$lambda.min,
                             type = "response"))

    r <- roc(y_ex, pr, direction = "<", quiet = TRUE)

    rows[[length(rows) + 1]] <- data.frame(
      K = K_LAB(k), genes_in = length(g), seed = seed,
      lambda_min = signif(cv$lambda.min, 5),
      n_nonzero = sum(abs(coef(cv$glmnet.fit, s = cv$lambda.min)) > 0),
      ercb_auc = round(as.numeric(auc(r)), 4),
      ci_low = round(as.numeric(ci.auc(r, method = "delong")[1]), 4),
      ci_high = round(as.numeric(ci.auc(r, method = "delong")[3]), 4),
      stringsAsFactors = FALSE)

    cat(sprintf("  K=%-4s draw %2d  AUC %.4f  nnz %3d\n", K_LAB(k), i,
                rows[[length(rows)]]$ercb_auc, rows[[length(rows)]]$n_nonzero))

  }

}

stab <- do.call(rbind, rows)

write.csv(stab, file.path(tab_dir, "05c_gate_seed_stability.csv"), row.names = FALSE)

suite <- do.call(rbind, lapply(split(stab, stab$K), function(x) {
  data.frame(K = x$K[1], draws = nrow(x),
             auc_median = round(median(x$ercb_auc), 4),
             auc_min = min(x$ercb_auc), auc_max = max(x$ercb_auc),
             auc_range = round(max(x$ercb_auc) - min(x$ercb_auc), 4),
             auc_sd = round(sd(x$ercb_auc), 4),
             ci_low_min = min(x$ci_low), ci_high_min = min(x$ci_high),
             stringsAsFactors = FALSE)
}))

cat("\n================ SEED STABILITY ================\n")
print(suite, row.names = FALSE)

p <- ggplot(stab, aes(x = factor(K, levels = K_LAB(K_GRID)), y = ercb_auc)) +
  geom_jitter(width = 0.12, height = 0, size = 2, colour = "#2471A3", alpha = 0.85) +
  stat_summary(fun = median, geom = "crossbar", width = 0.34, colour = "#C0392B", linewidth = 0.6) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
  labs(title = "The harmonised gate is stable across fold partitions",
       subtitle = sprintf("%d independent RNG streams per feature count; red bar = median; ERCB external cohort",
                          N_DRAW),
       x = "features used (K)", y = "ERCB AUC (external)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(p)

pdf(file.path(fig_dir, "05c_gate_seed_stability.pdf"), width = 6.4, height = 4.2)
print(p)
dev.off()

jpeg(file.path(fig_dir, "05c_gate_seed_stability.jpeg"), width = 6.4, height = 4.2,
     units = "in", res = 300, quality = 95, bg = "white")
print(p)
dev.off()

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 05c_gate_seed_stability.R"),
             capture.output(print(suite, row.names = FALSE))),
           file.path(PROJ, "logs", "gate_seed_stability.log"))

cat("\nDone.\n")

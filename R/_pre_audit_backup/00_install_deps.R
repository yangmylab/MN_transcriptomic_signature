# Paper 1 - install missing R/Bioconductor dependencies into the workspace-local library.
# Run with:  Rscript R/00_install_deps.R   (do NOT override R_LIBS_USER)
local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })
cat("Installing into :", WS_LIB, "\n")
cat("libPaths        :", paste(.libPaths(), collapse = " | "), "\n")
cat("Bioc version    :", as.character(BiocManager::version()), "\n")

options(repos = c(CRAN = "https://mirrors.aliyun.com/CRAN/"))

cran_pkgs <- c("WGCNA", "glmnet", "randomForest", "caret", "Boruta", "xgboost",
               "metafor", "rms", "ggpubr", "corrplot", "patchwork", "matrixStats",
               "reshape2", "data.table", "preprocessCore")
bioc_pkgs <- c("oligo", "ConsensusClusterPlus", "org.Hs.eg.db",
               "hugene21sttranscriptcluster.db", "hgu133plus2.db",
               "hgu133a.db", "hgu133b.db", "ComplexHeatmap", "msigdbr", "DOSE",
               "arrayQualityMetrics", "impute")

install_one <- function(p, bioc = FALSE) {
  if (requireNamespace(p, quietly = TRUE)) { cat("SKIP (present):", p, "\n"); return(invisible(TRUE)) }
  cat("==== installing", p, if (bioc) "(Bioc)" else "(CRAN)", "====\n")
  ok <- tryCatch({
    if (bioc) BiocManager::install(p, lib = WS_LIB, ask = FALSE, update = FALSE)
    else install.packages(p, lib = WS_LIB)
    requireNamespace(p, quietly = TRUE)
  }, error = function(e) { cat("ERROR", p, ":", conditionMessage(e), "\n"); FALSE })
  cat(if (isTRUE(ok)) "OK   " else "FAIL ", p, "\n")
  isTRUE(ok)
}

res <- c()
for (p in cran_pkgs) res[p] <- install_one(p, bioc = FALSE)
for (p in bioc_pkgs) res[p] <- install_one(p, bioc = TRUE)

want <- c(cran_pkgs, bioc_pkgs)
missing <- names(res)[!res]
cat("\n================ SUMMARY ================\n")
cat("installed OK :", paste(names(res)[res], collapse = ", "), "\n")
cat("FAILED       :", if (length(missing)) paste(missing, collapse = ", ") else "NONE", "\n")

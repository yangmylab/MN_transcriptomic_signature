#!/usr/bin/env Rscript
# Paper 1 - scRNA arm: install the packages the integration plan needs.
# Mirrors: Tsinghua CRAN + Bioconductor (works from China). Installs into .Rlib via R/_lib.R.

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(Ncpus = max(1, parallel::detectCores() - 1), timeout = 1800)

# core: everything the Paper-1 integration arm uses (07b-07f)
need <- c("scDblFinder", "miloR", "muscat", "decoupleR", "SoupX", "nnls", "quadprog")

# optional: cell-cell communication belongs to Paper 2; liana is not published for Bioc 3.23
optional <- c("liana")

have <- need[vapply(need, requireNamespace, logical(1), quietly = TRUE)]
miss <- setdiff(need, have)

cat("already installed:", if (length(have)) paste(have, collapse = ", ") else "none", "\n")
cat("to install       :", if (length(miss)) paste(miss, collapse = ", ") else "none", "\n")

if (length(miss)) {
  BiocManager::install(miss, update = FALSE, ask = FALSE)
}

miss_opt <- optional[!vapply(optional, requireNamespace, logical(1), quietly = TRUE)]

if (length(miss_opt)) {
  message("optional (not installed): ", paste(miss_opt, collapse = ", "))
}

ok <- need[vapply(need, requireNamespace, logical(1), quietly = TRUE)]

cat("\n================ INSTALL RESULT ================\n")

for (p in need) {
  cat(sprintf("%-14s %-6s %s\n", p, requireNamespace(p, quietly = TRUE),
              tryCatch(as.character(packageVersion(p)), error = function(e) "-")))
}

cat("\navailable:", length(ok), "of", length(need), "\n")

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 07z_install_scrna_pkgs.R"),
             capture.output(print(data.frame(package = need,
                                              installed = vapply(need, requireNamespace,
                                                                 logical(1), quietly = TRUE))))),
           file.path(PROJ, "logs", "scrna_install.log"))

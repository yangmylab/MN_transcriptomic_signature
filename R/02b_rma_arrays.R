#!/usr/bin/env Rscript
# Paper 1 - Phase A step 2b: RMA-normalize every array cohort of the APPROVED design, at gene level
# on the Brainarray Entrez-Gene definitions the original authors used.
#
#   GPL19983 (HuGene-2.1-ST) : affy cannot read this array type, so gene-level RMA is rebuilt from
#                              oligo's probe-level values using the Brainarray CDF probe indices
#                              (validated against the authors' own matrix, Pearson r = 0.986).
#                              The box has only 14 GB RAM and a 244-array cohort needs ~12 GB in one
#                              gulp, so samples are processed in CHUNKS and quantile-normalised to a
#                              reference distribution built in a first pass.
#   GPL22945/24120/19109/19184 (HG-U133) : affy + the matching Brainarray CDF, straight RMA.
#
# GSE108109 and GSE108112 are deliberately not processed: they are 100% redundant releases that the
# approved design drops.

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(affy); library(Biobase); library(oligo); library(preprocessCore); library(affyio) })

rma_dir   <- file.path(PROJ, "work", "rma")
cel_root  <- file.path(PROJ, "work", "cel")
cache_root <- file.path(PROJ, "work", "rma_cache")

for (d in c(rma_dir, cel_root, cache_root)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

NCPU   <- min(6, max(1, parallel::detectCores() - 1))
CHUNK  <- 40    # ~330 MB of probe matrix per chunk; keeps peak well under the 14 GB budget

# ---- helpers ------------------------------------------------------------------------------------

gunzip_bin <- function(src, dst) {
  ci <- gzfile(src, "rb"); co <- file(dst, "wb")
  repeat { ch <- readBin(ci, "raw", n = 1048576); if (!length(ch)) break; writeBin(ch, co) }
  close(ci); close(co)
}

decompress_cels <- function(acc) {
  d <- file.path(cel_root, acc)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  src <- list.files(file.path(DATA, acc, "raw"), pattern = "\\.CEL\\.gz$", full.names = TRUE)
  dst <- file.path(d, sub("\\.gz$", "", basename(src)))
  for (i in which(!file.exists(dst))) gunzip_bin(src[i], dst[i])
  dst
}

gsm_platform <- function(path) {
  l <- readLines(gzfile(path), warn = FALSE)
  pull <- function(k) {
    x <- l[grepl(paste0("^!", k), l)]
    if (!length(x)) return(NULL)
    gsub('^"|"$', "", strsplit(sub("^[^\t]*\t", "", x), "\t")[[1]])
  }
  g <- pull("Sample_geo_accession"); p <- pull("Sample_platform_id")
  if (is.null(g) || is.null(p) || length(g) != length(p)) return(NULL)
  data.frame(gsm = g, gpl = p, stringsAsFactors = FALSE)
}

# quantile-normalise one column onto a fixed reference distribution
qn_to_ref <- function(col, ref) {
  o <- order(col)
  v <- col
  v[o] <- ref
  v
}

# ---- gene-level RMA for HuGene-2.1-ST on Brainarray gene definitions (memory-bounded) ------------

hugene_gene_rma <- function(cels, acc, cdf_pkg = "hugene21sthsentrezgcdf") {
  if (!requireNamespace(cdf_pkg, quietly = TRUE)) stop("CDF package not installed: ", cdf_pkg)
  env <- get(cdf_pkg, envir = asNamespace(cdf_pkg))
  keys <- ls(env)

  cache <- file.path(cache_root, acc)
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)

  n <- length(cels)
  chunks <- split(seq_len(n), ceiling(seq_len(n) / CHUNK))

  # --- pass 1: background-correct + log2 per chunk, cache it, accumulate the quantile reference ---

  ref_sum <- NULL
  ref_n <- 0
  row_names <- NULL

  for (k in seq_along(chunks)) {
    f <- file.path(cache, sprintf("chunk_%03d.rds", k))
    if (!file.exists(f)) {
      cat("    chunk", k, "/", length(chunks), ": bg-correcting", length(chunks[[k]]), "arrays\n")
      fs <- oligo::read.celfiles(cels[chunks[[k]]], verbose = FALSE)
      bc <- oligo::backgroundCorrect(fs, method = "rma")
      rm(fs); invisible(gc(verbose = FALSE))
      m <- oligo::pm(bc, target = "probes")
      rm(bc); invisible(gc(verbose = FALSE))
      if (max(m, na.rm = TRUE) > 100) m <- log2(m)
      saveRDS(m, f)
      rm(m); invisible(gc(verbose = FALSE))
    }
    m <- readRDS(f)
    if (is.null(row_names)) row_names <- as.integer(rownames(m))
    s <- apply(m, 2, sort)
    ref_sum <- if (is.null(ref_sum)) rowSums(s) else ref_sum + rowSums(s)
    ref_n <- ref_n + ncol(m)
    rm(m, s); invisible(gc(verbose = FALSE))
  }

  ref <- ref_sum / ref_n
  rm(ref_sum); invisible(gc(verbose = FALSE))

  stopifnot(!anyNA(row_names))
  cat("    quantile reference built from", ref_n, "arrays over", length(ref), "probes\n")

  # --- pre-resolve each gene's probe rows (identical in every chunk) -------------------------------

  # map every CEL cell index to its matrix row ONCE. match() against a 1.02M-element vector rebuilds
  # a hash table on every call, which cost 436 s for 25,582 genes and dominated the entire run.
  pos_of <- integer(max(row_names))
  pos_of[row_names] <- seq_along(row_names)

  rows <- lapply(keys, function(g) {
    idx <- as.integer(env[[g]][, "pm"]); idx <- idx[!is.na(idx)]
    pos <- pos_of[idx]; pos[pos > 0]
  })

  # --- pass 2: normalise each chunk to the reference and median-polish per gene -------------------

  out <- matrix(NA_real_, length(keys), n, dimnames = list(sub("_at$", "", keys), basename(cels)))

  # median polish per gene (the RMA summarisation step). preprocessCore::rcModelPLM returns probe
  # and sample effects CONCATENATED (length nrow+ncol), and the medians form below is the variant
  # that was validated against the authors' own matrix (Pearson r = 0.986).
  est_one <- function(i, m) {
    pos <- rows[[i]]
    if (length(pos) < 3) return(rep(NA_real_, ncol(m)))
    mp <- stats::medpolish(m[pos, , drop = FALSE], na.rm = TRUE, trace.iter = FALSE)
    mp$overall + mp$col
  }

  for (k in seq_along(chunks)) {
    cat("    chunk", k, "/", length(chunks), ": summarising to genes\n")
    m <- readRDS(file.path(cache, sprintf("chunk_%03d.rds", k)))
    for (j in seq_len(ncol(m))) m[, j] <- qn_to_ref(m[, j], ref)
    est <- parallel::mclapply(seq_along(keys), est_one, m = m, mc.cores = NCPU)
    out[, chunks[[k]]] <- do.call(rbind, est)
    rm(m, est); invisible(gc(verbose = FALSE))
  }

  out
}

# ---- affy RMA for the U133 cohorts --------------------------------------------------------------

affy_gene_rma <- function(cels, cdf) {
  ab <- ReadAffy(filenames = cels, cdfname = cdf, verbose = FALSE)
  e <- affy::rma(ab, verbose = FALSE)
  m <- exprs(e)
  rownames(m) <- sub("_at$", "", rownames(m))
  m
}

save_mat <- function(m, acc, tag) {
  f <- file.path(rma_dir, paste0(acc, "_", tag, "_gene_rma.txt.gz"))
  gz <- gzfile(f, "w")
  write.table(data.frame(gene = rownames(m), m, check.names = FALSE), gz, sep = "\t",
              quote = FALSE, row.names = FALSE)
  close(gz)
  saveRDS(m, file.path(rma_dir, paste0(acc, "_", tag, "_gene_rma.rds")))
  cat("  saved", basename(f), ":", nrow(m), "x", ncol(m), "\n")
}

# ---- run: HuGene-2.1-ST cohorts ------------------------------------------------------------------

cat("=== HuGene-2.1-ST cohorts (gene-level RMA on Brainarray definitions) ===\n")

for (acc in c("GSE133288", "GSE200818", "GSE200828")) {
  out_f <- file.path(rma_dir, paste0(acc, "_GPL19983_gene_rma.rds"))
  if (file.exists(out_f)) { cat(acc, ": already done\n"); next }
  cat(acc, ": decompressing CEL files\n")
  cel <- decompress_cels(acc)
  t0 <- Sys.time()
  m <- hugene_gene_rma(cel, acc)
  colnames(m) <- sub("_.*$", "", basename(cel))
  cat(" ", acc, "elapsed", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
  save_mat(m, acc, "GPL19983")
  rm(m); invisible(gc(verbose = FALSE))
}

# ---- run: U133 cohorts (split by GPL, affy RMA) --------------------------------------------------

cat("\n=== HG-U133 cohorts (affy + Brainarray CDF) ===\n")

u133 <- list(
  GSE104948 = c(PLUS2 = "hgu133plus2hsentrezgcdf", U133A = "hgu133ahsentrezgcdf"),
  GSE104954 = c(PLUS2 = "hgu133plus2hsentrezgcdf", U133A = "hgu133ahsentrezgcdf"),
  GSE99339  = c(PLUS2 = "hgu133plus2hsentrezg2cdf", U133A = "hgu133ahsentrezg2cdf"),
  GSE99325  = c(PLUS2 = "hgu133plus2hsentrezg2cdf", U133A = "hgu133ahsentrezg2cdf")
)

# the series matrices of a two-platform release list OVERLAPPING sample ids, so the platform is read
# from each CEL's own header instead (HG-U133A and HG-U133 Plus 2 have different dimensions)
cel_platform <- function(f) {
  h <- affyio::read.celfile.header(f, info = "basic")
  if (grepl("PLUS", toupper(h$cdfName))) "PLUS2" else "U133A"
}

for (acc in names(u133)) {
  out_f <- file.path(rma_dir, paste0(acc, "_gene_rma_merged.rds"))
  if (file.exists(out_f)) { cat(acc, ": already done\n"); next }
  cel <- decompress_cels(acc)
  gsm <- sub("_.*$", "", basename(cel))
  plat <- vapply(cel, cel_platform, "")
  cat(acc, ":", paste(names(table(plat)), table(plat), sep = "=", collapse = " "), "\n")
  mats <- list()
  for (p in unique(plat)) {
    idx <- which(plat == p)
    mats[[p]] <- affy_gene_rma(cel[idx], u133[[acc]][[p]])
    colnames(mats[[p]]) <- gsm[idx]
    cat("  ", acc, p, "->", nrow(mats[[p]]), "genes x", ncol(mats[[p]]), "samples\n")
  }
  # the two array types of a release carry the same samples but different gene sets, so align on the
  # UNION of genes and fill NA where a gene is absent from one platform
  genes <- Reduce(union, lapply(mats, rownames))
  smp <- unlist(lapply(mats, colnames))
  merged <- matrix(NA_real_, length(genes), length(smp), dimnames = list(genes, smp))
  for (k in seq_along(mats)) merged[rownames(mats[[k]]), colnames(mats[[k]])] <- mats[[k]]
  saveRDS(merged, out_f)
  cat("  merged:", nrow(merged), "x", ncol(merged), "\n")
  rm(mats, merged); invisible(gc(verbose = FALSE))
}

cat("\n=========== DONE ===========\n")

print(list.files(rma_dir))

#!/usr/bin/env Rscript
# Paper 1 - Phase A step 0c: repair the corrupt GSE197307 Ensembl IDs and put the three RNA-seq
# count matrices into a uniform "processed" layout (one clean, tab-separated, gzipped matrix each).

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

log_dir <- file.path(PROJ, "logs")

dir.create(log_dir, showWarnings = FALSE)

# ---- GSE197307 (glomeruli, NEPTUNE) ------------------------------------------------------------
# As published, the supplementary file has a mangled first column: the header reads
# "ENSEMBL_ene_ID" and EVERY Ensembl gene ID lost its "G" (ENS00000000003 vs ENSG00000000003).
# The corruption is uniform (every row matches ^ENS[0-9]+$), so the repair is deterministic.

p197_in <- file.path(DATA, "GSE197307", "GSE197307_NEPTUNE_Glom_RNA_counts_2022.txt.gz")

p197_out <- file.path(DATA, "GSE197307", "processed", "GSE197307_counts_ensg.txt.gz")

dir.create(dirname(p197_out), showWarnings = FALSE, recursive = TRUE)

x197 <- read.delim(gzfile(p197_in), check.names = FALSE, stringsAsFactors = FALSE)

cat("GSE197307 as read:", nrow(x197), "genes x", ncol(x197) - 1, "samples\n")

cat("  raw header[1]:", colnames(x197)[1], "\n")

ids197 <- x197[[1]]

stopifnot(all(grepl("^ENS[0-9]+$", ids197)))

x197[[1]] <- sub("^ENS", "ENSG", ids197)

colnames(x197)[1] <- "ensembl_gene_id"

stopifnot(!anyDuplicated(x197[[1]]), all(grepl("^ENSG[0-9]+$", x197[[1]])))

gz <- gzfile(p197_out, "w")

write.table(x197, gz, sep = "\t", quote = FALSE, row.names = FALSE)

close(gz)

cat("  wrote", p197_out, "(", nrow(x197), "x", ncol(x197), ")\n")

# ---- GSE182380 (tubulointerstitium, NEPTUNE) and GSE216841 (glomeruli, LCM-FFPE) ---------------
# Both are already well-formed: 182380 is keyed by HGNC symbol, 216841 by Ensembl gene ID with the
# symbol and description as extra columns. Re-emit them so every RNA-seq matrix sits in one place
# with one naming convention.

write_gz <- function(d, out_gz) {
  gz <- gzfile(out_gz, "w")
  write.table(d, gz, sep = "\t", quote = FALSE, row.names = FALSE)
  close(gz)
}

p182_out <- file.path(DATA, "GSE182380", "processed", "GSE182380_counts_symbol.txt.gz")

p216_out <- file.path(DATA, "GSE216841", "processed", "GSE216841_counts_ensg.txt.gz")

dir.create(dirname(p182_out), showWarnings = FALSE, recursive = TRUE)

dir.create(dirname(p216_out), showWarnings = FALSE, recursive = TRUE)

d182 <- read.delim(gzfile(file.path(DATA, "GSE182380", "GSE182380_NEPTUNE_TI_adult_RNA_counts_2021.txt.gz")), check.names = FALSE, stringsAsFactors = FALSE)

colnames(d182)[1] <- "hgnc_symbol"

write_gz(d182, p182_out)

cat(" GSE182380:", nrow(d182), "genes x", ncol(d182) - 1, "samples -> ", p182_out, "\n")

d216 <- read.delim(gzfile(file.path(DATA, "GSE216841", "GSE216841_MNd_ncounts_annot.txt.gz")), check.names = FALSE, stringsAsFactors = FALSE)

# column 2 merely repeats column 1: the file ships both "rowname" and an identical
# "ensembl_gene_id", so drop the duplicate and keep the first as the canonical key.

d216 <- d216[, -2]

colnames(d216)[1] <- "ensembl_gene_id"

write_gz(d216, p216_out)

cat(" GSE216841:", nrow(d216), "genes x", ncol(d216) - 3, "samples -> ", p216_out, "\n")

# ---- Ensembl -> symbol map, harvested from GSE216841 (no new download needed) -------------------

map <- unique(d216[, c("ensembl_gene_id", "hgnc_GRCh38p12")])

map <- map[nzchar(map$hgnc_GRCh38p12), ]

colnames(map) <- c("ensembl_gene_id", "hgnc_symbol")

map_path <- file.path(PROJ, "metadata", "ensg2symbol.tsv")

dir.create(dirname(map_path), showWarnings = FALSE, recursive = TRUE)

write.table(map, map_path, sep = "\t", quote = FALSE, row.names = FALSE)

cov197 <- mean(x197$ensembl_gene_id %in% map$ensembl_gene_id)

cat("ensg2symbol map:", nrow(map), "entries;", sprintf("%.1f%%", 100 * cov197), "of GSE197307 genes covered\n")

# ---- Run log -----------------------------------------------------------------------------------

writeLines(c(
  paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 00c_repair_counts.R"),
  paste0("  GSE197307: " , nrow(x197), " genes x ", ncol(x197) - 1, " samples -> ", p197_out),
  paste0("  GSE182380: ", nrow(d182), " genes x ", ncol(d182) - 1, " samples -> ", p182_out),
  paste0("  GSE216841: ", nrow(d216), " genes x ", ncol(d216) - 1, " samples -> ", p216_out),
  paste0("  ensg2symbol: ", nrow(map), " entries -> ", map_path),
  paste0("  GSE197307 symbol coverage: ", sprintf("%.1f%%", 100 * cov197))
), file.path(log_dir, "repair_counts.log"))

cat("Done.\n")

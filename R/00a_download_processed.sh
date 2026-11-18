#!/usr/bin/env bash
# Paper 1 - Track A: processed-data download, PARALLEL (8 streams -> aggregates the host's
# capped per-connection bandwidth; measured 8-conn ~235 KB/s vs ~25 KB/s single).
# Series matrices + NEPTUNE processed matrices + RNA-seq counts (~190 MB total).
set -u
ROOT="${PAPER1_PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
DEST="$(dirname "$ROOT")/Datasets"; LOGD="$ROOT/logs"
BASE=https://ftp.ncbi.nlm.nih.gov/geo/series
mkdir -p "$DEST" "$LOGD"
LOG="$LOGD/download_processed.log"
LIST="$ROOT/logs/trackA_urls.txt"
export MNFETCH_LOG="$LOG"

p(){ echo "$1" | sed -E 's/[0-9]{3}$/nnn/'; }

: > "$LIST"
add(){ echo "$BASE/$(p "$2")/$2/$1 $DEST/$2" >> "$LIST"; }

for g in GSE133288 GSE200828 GSE200818 GSE108109 GSE108112 GSE73953 GSE197307 GSE182380 GSE216841; do
  add "matrix/${g}_series_matrix.txt.gz" "$g"
done
for g in GSE104948 GSE104954; do
  add "matrix/${g}-GPL22945_series_matrix.txt.gz" "$g"
  add "matrix/${g}-GPL24120_series_matrix.txt.gz" "$g"
done
for g in GSE99339 GSE99325; do
  add "matrix/${g}-GPL19109_series_matrix.txt.gz" "$g"
  add "matrix/${g}-GPL19184_series_matrix.txt.gz" "$g"
done
add "suppl/GSE200828_processed_data.txt.gz" GSE200828
add "suppl/GSE200818_processed_data.txt.gz" GSE200818
add "suppl/GSE197307_NEPTUNE_Glom_RNA_counts_2022.txt.gz" GSE197307
add "suppl/GSE182380_NEPTUNE_TI_adult_RNA_counts_2021.txt.gz" GSE182380
add "suppl/GSE216841_MNd_ncounts_annot.txt.gz" GSE216841

echo "[$(date '+%F %T')] TRACK A start: $(wc -l < "$LIST") files" | tee -a "$LOG"
xargs -P 8 -n 2 -a "$LIST" bash "$ROOT/R/_fetch.sh" 2>&1 | tee -a "$LOG"
echo "[$(date '+%F %T')] TRACK A DONE" | tee -a "$LOG"
du -sh "$DEST" | tee -a "$LOG"

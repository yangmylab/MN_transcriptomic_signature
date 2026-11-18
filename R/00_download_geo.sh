#!/usr/bin/env bash
# Paper 1 - Phase A step 1: download GEO datasets (microarray RAW.tar + RNA-seq counts + series matrices)
# Datasets live in the SHARED MN collection <MN>/Datasets/<GSE>/ (moved 2026-10-19; Paper_2
# and Paper_3 reuse several accessions). This mirrors DATA in R/_lib.R.
set -u
ROOT="${PAPER1_PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
DEST="$(dirname "$ROOT")/Datasets"
LOGD="$ROOT/logs"
BASE=https://ftp.ncbi.nlm.nih.gov/geo/series
mkdir -p "$DEST" "$LOGD"
LOG="$LOGD/download_geo.log"

say(){ echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

# GSE -> prefix dir (GSE133 -> GSE133nnn)
pref(){ echo "$1" | sed -E 's/[0-9]{3}$/nnn/'; }

fetch(){ # $1 = url, $2 = outdir
  local url="$1" out="$2" f
  f="$out/$(basename "$url")"
  if [[ -s "$f" ]]; then say "SKIP (exists) $(basename "$f")"; return 0; fi
  mkdir -p "$out"
  say "GET  $url"
  wget -c -q --tries=6 --timeout=60 --waitretry=10 -O "$f.part" "$url" \
    && mv "$f.part" "$f" && say "OK   $(basename "$f") $(du -h "$f" | cut -f1)" \
    || { say "FAIL $url"; return 1; }
}

MICROARRAY=(GSE133288 GSE200828 GSE200818 GSE108109 GSE108112 GSE104948 GSE104954 GSE99339 GSE99325 GSE73953)

for g in "${MICROARRAY[@]}"; do
  p=$(pref "$g")
  fetch "$BASE/$p/$g/suppl/${g}_RAW.tar" "$DEST/$g"
  fetch "$BASE/$p/$g/matrix/${g}_series_matrix.txt.gz" "$DEST/$g"
  fetch "$BASE/$p/$g/suppl/filelist.txt" "$DEST/$g"
done

# additional processed matrices supplied by NEPTUNE
fetch "$BASE/GSE200nnn/GSE200828/suppl/GSE200828_processed_data.txt.gz" "$DEST/GSE200828"
fetch "$BASE/GSE200nnn/GSE200818/suppl/GSE200818_processed_data.txt.gz" "$DEST/GSE200818"

# RNA-seq count matrices
fetch "$BASE/GSE197nnn/GSE197307/suppl/GSE197307_NEPTUNE_Glom_RNA_counts_2022.txt.gz" "$DEST/GSE197307"
fetch "$BASE/GSE182nnn/GSE182380/suppl/GSE182380_NEPTUNE_TI_adult_RNA_counts_2021.txt.gz" "$DEST/GSE182380"
fetch "$BASE/GSE216nnn/GSE216841/suppl/GSE216841_MNd_ncounts_annot.txt.gz" "$DEST/GSE216841"
fetch "$BASE/GSE197nnn/GSE197307/matrix/GSE197307_series_matrix.txt.gz" "$DEST/GSE197307"
fetch "$BASE/GSE182nnn/GSE182380/matrix/GSE182380_series_matrix.txt.gz" "$DEST/GSE182380"
fetch "$BASE/GSE216nnn/GSE216841/matrix/GSE216841_series_matrix.txt.gz" "$DEST/GSE216841"

say "DONE. Total size:"; du -sh "$DEST" | tee -a "$LOG"; find "$DEST" -type f | wc -l | xargs echo "files:" | tee -a "$LOG"

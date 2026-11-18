#!/usr/bin/env bash
# Paper 1 - Track B: raw CEL tarball download (~7.5 GB). Long-running background job given the
# ~60-100 KB/s uplink to NCBI from this host (~20 h). Series matrices (Track A) unblock the
# analysis meanwhile; these raw files allow an optional raw-level RMA recomputation later.
set -u
ROOT="${PAPER1_PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
DEST="$(dirname "$ROOT")/Datasets"; LOGD="$ROOT/logs"
BASE=https://ftp.ncbi.nlm.nih.gov/geo/series
mkdir -p "$DEST" "$LOGD"
LOG="$LOGD/download_raw.log"
say(){ echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
pref(){ echo "$1" | sed -E 's/[0-9]{3}$/nnn/'; }

fetch(){
  local url="$1" out="$2" f
  f="$out/$(basename "$url")"
  if [[ -s "$f" ]]; then say "SKIP $f"; return 0; fi
  mkdir -p "$out"
  wget -c -q --tries=20 --timeout=120 --waitretry=30 --read-timeout=120 -O "$f.part" "$url" \
    && mv "$f.part" "$f" && say "OK   $f ($(du -h "$f"|cut -f1))" || say "FAIL $url"
}

# discovery first, then validation, so partial completion yields usable raw data soonest
ORDER=(GSE133288 GSE200828 GSE200818 GSE108109 GSE108112 GSE104948 GSE104954 GSE99339 GSE99325 GSE73953)
for g in "${ORDER[@]}"; do
  fetch "$BASE/$(pref $g)/$g/suppl/${g}_RAW.tar" "$DEST/$g"
  fetch "$BASE/$(pref $g)/$g/suppl/filelist.txt" "$DEST/$g"
done
say "TRACK B DONE"; du -sh "$DEST" | tee -a "$LOG"

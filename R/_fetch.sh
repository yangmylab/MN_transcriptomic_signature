#!/usr/bin/env bash
# Single-file resumable fetch used by the parallel download drivers.
# usage: _fetch.sh <url> <outdir>
set -u
url="${1:-}"; out="${2:-}"
if [[ -z "$url" || -z "$out" ]]; then echo "USAGE: _fetch.sh <url> <outdir>"; exit 2; fi
f="$out/$(basename "$url")"
mkdir -p "$out"
if [[ -s "$f" ]]; then echo "SKIP $f"; exit 0; fi
if curl -sS --fail --location --retry 10 --retry-delay 5 --retry-all-errors \
        -C - -m 7200 --connect-timeout 30 -o "$f.part" "$url"; then
  mv "$f.part" "$f"
  echo "OK   $f ($(du -h "$f" | cut -f1))"
else
  echo "FAIL $url"; exit 1
fi

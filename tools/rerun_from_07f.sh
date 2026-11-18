#!/usr/bin/env bash
# Continue the 2026-10-20 audit re-run from 07f (05a/05b/06a already completed).
set -u
# Resolve the project root from this script's own location: the volume is remounted under a
# different user dir between sessions (/run/media/{meng,yang}/MY_Data), so a hard-coded cd breaks.
ROOT="${PAPER1_PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$ROOT" || exit 1
LOG=logs/audit_rerun.log

step() {
  local script=$1; shift
  echo "=== $(date '+%F %T') START $script $*" | tee -a "$LOG"
  Rscript "R/$script" "$@" >> "$LOG" 2>&1
  local rc=$?
  echo "=== $(date '+%F %T') END $script exit=$rc" | tee -a "$LOG"
  return $rc
}
fail() { echo "=== ABORT: $1 failed" | tee -a "$LOG"; exit 1; }

step 07f_deconvolution.R         || fail 07f
step 08a_wgcna_network.R         || fail 08a
step 08b_module_trait.R          || fail 08b
step 08c_wgcna_stability.R       || fail 08c
step 08d_module_preservation.R   || fail 08d
mv -f permutedStats-actualModules_partial.RData work/wgcna/ 2>/dev/null
step 08e_hub_genes.R             || fail 08e
step 08f_module_enrichment.R --force || fail 08f
step 09a_panel_selection.R       || fail 09a
step 09b_panel_validation.R      || fail 09b
step 09c_panel_robustness.R      || fail 09c
step 09d_panel_mechanism.R       || fail 09d
step 09e_panel_sensitivity.R     || fail 09e
step 09f_pla2r_strata.R          || fail 09f

echo "=== $(date '+%F %T') ALL STEPS DONE" | tee -a "$LOG"

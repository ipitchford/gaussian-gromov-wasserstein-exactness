#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
export MPLCONFIGDIR="${TMPDIR:-/tmp}/gaussian-gw-matplotlib"
# Every checker raises SystemExit on failure (also under python -O); with
# set -e the build stops at the first failed check. Outputs keep the shipped
# filenames in this directory (and results/ for verify_gw.py).
python3 verify_gw.py --out results > results/run_log.json
python3 independent_symbolic_audit.py > independent_symbolic_audit.json
python3 verify_combined_flow.py --out . > combined_flow_run.log
python3 proof_digest_checks.py --out . > proof_digest_checks_run.log
python3 gw_gaussian_bracket.py > gw_gaussian_bracket_demo.log
python3 extension_checks.py --out . > extension_checks_run.log
python3 negative_controls.py --out . > negative_controls_run.log
bash build_pdf.sh

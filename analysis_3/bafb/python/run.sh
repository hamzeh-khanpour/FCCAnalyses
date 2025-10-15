#!/usr/bin/env bash
# Glue: Stage-1 (signed-cosÎ¸, JC3) -> Stage-2 (Poisson fit) -> Plots
# Exclusive Z->bb sample, no backgrounds.

set -euo pipefail

# -------------------------- helpers --------------------------
has_flag() { local f="$1"; shift; for a in "$@"; do [[ "$a" == "$f" || "$a" == $f=* ]] && return 0; done; return 1; }
flag_value_next() { local f="$1"; shift; local i; for ((i=1;i<=$#;i++)); do
  if [[ "${!i}" == "$f" ]]; then local j=$((i+1)); (( j <= $# )) && { eval "echo \${$j}"; return 0; }; fi; done; echo ""; }
print_cmd() { printf '   %q ' "$@"; echo; }

# ---------------------- env / defaults -----------------------
export PYTHONPATH="$(pwd)/analysis_3/bafb/python:${PYTHONPATH:-}"

# If no args, run with --test
if [[ $# -eq 0 ]]; then set -- --test; fi

# Normalize legacy flag
ARGS=( "$@" )
for i in "${!ARGS[@]}"; do
  [[ "${ARGS[$i]}" == "--files-list" ]] && ARGS[$i]="--input-file-list"
done

# Threads: prefer CLI; else env BABFB_THREADS; avoid dup
NTHREADS_ENV=()
if ! has_flag "--n-threads" "${ARGS[@]}" && [[ -n "${BABFB_THREADS:-}" ]]; then
  NTHREADS_ENV=( --n-threads "${BABFB_THREADS}" )
fi
echo ">> [run] Threads requested: ${NTHREADS_ENV[*]:-(from CLI or managed)}"

# -------------------------- Stage-1 --------------------------
# Honor user --output; else default under outputs/
OUT1="outputs/bafb/stage1/bafb_stage1.root"
user_out=false
for i in "${!ARGS[@]}"; do
  case "${ARGS[$i]}" in
    --output) val="$(flag_value_next --output "${ARGS[@]}")"; [[ -n "$val" ]] && OUT1="$val" && user_out=true ;;
    --output=*) OUT1="${ARGS[$i]#--output=}"; user_out=true ;;
  esac
done
mkdir -p "$(dirname "$OUT1")"

# Build Stage-1 command (only recognized fccanalysis flags)
STAGE1_CMD=( fccanalysis run analysis_3/bafb/python/analysis_stage1.py "${ARGS[@]}" )
$user_out || STAGE1_CMD+=( --output "$OUT1" )
# Append env threads only if user didn't pass one
(( ${#NTHREADS_ENV[@]} > 0 )) && STAGE1_CMD+=( "${NTHREADS_ENV[@]}" )

echo ">> [run] Stage-1 output: $OUT1"
echo ">> [run] Running Stage-1:"
print_cmd "${STAGE1_CMD[@]}"
"${STAGE1_CMD[@]}"

# Verify Stage-1 output
if [[ ! -f "$OUT1" ]]; then
  echo "!! [run] Stage-1 output not found: $OUT1" >&2
  exit 2
fi

# -------------------------- Stage-2 --------------------------
mkdir -p outputs/bafb/stage2/results outputs/bafb/stage2/plots

STAGE2_CMD=( fccanalysis run analysis_3/bafb/python/analysis_stage2.py --input "$OUT1" )
# (Optional) If you trust pass-through, you may add analysis flags here later.

# Append env threads only if user didn't pass one
if ! has_flag "--n-threads" "${ARGS[@]}" && (( ${#NTHREADS_ENV[@]} > 0 )); then
  STAGE2_CMD+=( "${NTHREADS_ENV[@]}" )
fi

echo ">> [run] Running Stage-2 on: $OUT1"
print_cmd "${STAGE2_CMD[@]}"
"${STAGE2_CMD[@]}"

# --------------------------- Plots ---------------------------
ROOT_OUT="outputs/bafb/stage2/afbb_stage2.root"
JSON_OUT="outputs/bafb/stage2/results/AFBb_fit.json"

if [[ -f "$ROOT_OUT" && -f "$JSON_OUT" ]]; then
  PLOT_CMD=( python3 analysis_3/bafb/python/analysis_plots.py
             --root "$ROOT_OUT" --json "$JSON_OUT" --out outputs/bafb/stage2/plots )
  echo ">> [run] Generating plots"
  print_cmd "${PLOT_CMD[@]}"
  "${PLOT_CMD[@]}"
else
  echo "!! [run] Skipping plots; missing:"
  [[ -f "$ROOT_OUT" ]] || echo "   - $ROOT_OUT"
  [[ -f "$JSON_OUT" ]] || echo "   - $JSON_OUT"
  exit 3
fi

# -------------------------- Summary --------------------------
echo ">> [run] Done. Outputs:"
echo "- Stage1: ${OUT1}"
echo "- Stage2: outputs/bafb/stage2/afbb_stage2.root"
echo "- Results: outputs/bafb/stage2/results/AFBb_fit.json (and .txt)"
echo "- Plots  : outputs/bafb/stage2/plots/AFBb_signed_fit.png/pdf"

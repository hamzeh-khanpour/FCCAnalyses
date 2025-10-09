#!/usr/bin/env bash
# Glue script: runs stage-1 -> stage-2 -> plots.
# - Ensures analysis modules import cleanly
# - Defaults to --test when no args
# - Translates legacy --files-list -> --input-file-list
# - Forwards --n-threads to BOTH stages (env fallback: BABFB_THREADS)
# - Respects a user-provided --output for Stage-1, otherwise uses a canonical path
# - Echoes commands for reproducibility and does minimal sanity checks

set -euo pipefail

# Make our Python modules importable
export PYTHONPATH="$(pwd)/analysis_2/bafb/python:${PYTHONPATH:-}"

# Default to --test if no args
if [ $# -eq 0 ]; then
  set -- --test
fi

# Normalize args (legacy flag)
ARGS=("$@")
for i in "${!ARGS[@]}"; do
  [[ "${ARGS[$i]}" == "--files-list" ]] && ARGS[$i]="--input-file-list"
done

# Extract/compose --n-threads flag (CLI has priority, env fallback BABFB_THREADS)
NTHREADS_FLAG=()
found_threads=false
for i in "${!ARGS[@]}"; do
  case "${ARGS[$i]}" in
    --n-threads)
      if (( i+1 < ${#ARGS[@]} )); then
        NTHREADS_FLAG=(--n-threads "${ARGS[$((i+1))]}")
        found_threads=true
      fi
      ;;
    --n-threads=*)
      NTHREADS_FLAG=("${ARGS[$i]}")
      found_threads=true
      ;;
  esac
done
if ! $found_threads && [[ -n "${BABFB_THREADS:-}" ]]; then
  NTHREADS_FLAG=(--n-threads "${BABFB_THREADS}")
fi

# Determine Stage-1 output path:
# - if user already provided --output/--output=..., honor it
# - else write to the canonical path
OUT1="outputs/bafb/stage1/bafb_stage1.root"
user_out=false
for i in "${!ARGS[@]}"; do
  case "${ARGS[$i]}" in
    --output)
      if (( i+1 < ${#ARGS[@]} )); then
        OUT1="${ARGS[$((i+1))]}"
        user_out=true
      fi
      ;;
    --output=*)
      OUT1="${ARGS[$i]#--output=}"
      user_out=true
      ;;
  esac
done
mkdir -p "$(dirname "$OUT1")"

echo ">> [run] Threads requested: ${NTHREADS_FLAG[*]:-(managed by framework)}"
echo ">> [run] Stage-1 output: $OUT1"

# Build Stage-1 command
STAGE1_CMD=(fccanalysis run analysis_2/bafb/python/analysis_stage1.py "${ARGS[@]}")
if ! $user_out; then
  STAGE1_CMD+=(--output "$OUT1")
fi
# Always add threads flag if present (even when user provided --output)
if (( ${#NTHREADS_FLAG[@]} > 0 )); then
  STAGE1_CMD+=("${NTHREADS_FLAG[@]}")
fi

echo ">> [run] Running stage-1:"
printf '   %q ' "${STAGE1_CMD[@]}"; echo
"${STAGE1_CMD[@]}"

# Minimal sanity check on Stage-1 output
if [[ ! -f "$OUT1" ]]; then
  echo "!! [run] Stage-1 output not found: $OUT1" >&2
  exit 2
fi

echo ">> [run] Running stage-2 on: $OUT1"
STAGE2_CMD=(fccanalysis run analysis_2/bafb/python/analysis_stage2.py --input "$OUT1")
if (( ${#NTHREADS_FLAG[@]} > 0 )); then
  STAGE2_CMD+=("${NTHREADS_FLAG[@]}")
fi
printf '   %q ' "${STAGE2_CMD[@]}"; echo
"${STAGE2_CMD[@]}"

echo ">> [run] Making plots (standalone script)"
PLOT_CMD=(python3 analysis_2/bafb/python/analysis_plots.py)
printf '   %q ' "${PLOT_CMD[@]}"; echo
"${PLOT_CMD[@]}"

echo ">> [run] Done. Outputs:"
echo "   - Stage1: $OUT1"
echo "   - Stage2: outputs/bafb/stage2/bafb_outputs.root and summary.json"
echo "   - Plots : outputs/bafb/plots/*.png and *.pdf"

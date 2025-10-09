#!/usr/bin/env bash
# Glue script: runs stage-1 -> stage-2 -> plots.
# - Ensures our analysis modules import cleanly
# - Defaults to --test when no args
# - Translates legacy --files-list -> --input-file-list
# - Writes stage-1 output to a known path and feeds it to stage-2
set -euo pipefail

# Ensure our analysis modules are importable
export PYTHONPATH="$(pwd)/analysis_1/bafb/python:${PYTHONPATH:-}"

# Default to --test if no args provided
if [ $# -eq 0 ]; then
  set -- --test
fi

# Normalize args (legacy flag)
ARGS=("$@")
for i in "${!ARGS[@]}"; do
  [[ "${ARGS[$i]}" == "--files-list" ]] && ARGS[$i]="--input-file-list"
done

# Extract --n-threads (if provided) so we can forward it to both stages
NTHREADS_FLAG=()
for i in "${!ARGS[@]}"; do
  case "${ARGS[$i]}" in
    --n-threads)
      if (( i+1 < ${#ARGS[@]} )); then
        NTHREADS_FLAG=(--n-threads "${ARGS[$((i+1))]}")
      fi
      ;;
    --n-threads=*)
      NTHREADS_FLAG=("${ARGS[$i]}")
      ;;
  esac
done

# Stage-1 output path
OUT1="outputs/bafb/stage1/bafb_stage1.root"
mkdir -p "$(dirname "$OUT1")"

echo ">> [run] Stage-1 output: $OUT1"
echo ">> [run] Running stage-1 with args: ${ARGS[*]}"
fccanalysis run analysis_1/bafb/python/analysis_stage1.py "${ARGS[@]}" --output "$OUT1"

echo ">> [run] Running stage-2 on: $OUT1"
fccanalysis run analysis_1/bafb/python/analysis_stage2.py "${NTHREADS_FLAG[@]}" --input "$OUT1"

echo ">> [run] Making plots (standalone script)"
python3 analysis_1/bafb/python/analysis_plots.py

echo ">> [run] Done. Outputs:"
echo "   - Stage1: $OUT1"
echo "   - Stage2: outputs/bafb/stage2/bafb_outputs.root and summary.json"
echo "   - Plots : outputs/bafb/plots/*.png"

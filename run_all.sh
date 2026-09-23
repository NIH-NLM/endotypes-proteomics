#!/usr/bin/env bash
# The machine path. The same notebooks are the human path in JupyterLab.
#
#   ./run_all.sh              steps 00-07 and 10-16, the pipeline
#   ./run_all.sh --optional   also runs 08 and 09, the projections
#
# Notebooks execute in place, so the committed file carries its own output.
#
# Every step reads cohorts/ (committed) and writes data/run_artifacts/
# (gitignored, regenerable). Step 00 is the exception: it is the one step that
# writes cohorts/, and it will not redraw a split that already exists.
set -euo pipefail
cd "$(dirname "$0")/ipynb"

STEPS=(00_prepare_data 01_soft_threshold 02_modules 03_eigenproteins \
       04_module_traits 05_endotypes 06_heatmap 07_federation 10_federated_modules \
       11_modules_per_cohort 12_panels_per_cohort 13_cluster_both_axes \
       14_heatmaps 15_project_and_federate 16_cohort_diagnostics)
[[ "${1:-}" == "--optional" ]] && STEPS+=(08_project_healthy 09_project_timepoints)

for nb in "${STEPS[@]}"; do
  printf '%-24s ' "$nb"
  if jupyter nbconvert --to notebook --execute --inplace \
       --ExecutePreprocessor.timeout=7200 "$nb.ipynb" >/dev/null 2>"/tmp/$nb.err"; then
    echo "OK"
  else
    echo "FAILED"; tail -20 "/tmp/$nb.err"; exit 1
  fi
done
echo "all steps completed"

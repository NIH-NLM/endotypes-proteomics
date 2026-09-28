#!/usr/bin/env bash
# The machine path. The same notebooks are the human path in JupyterLab.
#
#   ./run_all.sh              every step, 00-20
#
# Notebooks execute in place, so the committed file carries its own output.
#
# Every step reads cohorts/ (committed) and writes data/run_artifacts/
# (gitignored, regenerable). Step 00 is the exception: it is the one step that
# writes cohorts/, and it will not redraw a split that already exists.
set -euo pipefail
cd "$(dirname "$0")/ipynb"

STEPS=(00_prepare_data 01_soft_threshold 01b_interferon_panel 02_modules \
       02b_interferon 03_eigenproteins \
       04_module_traits 05_endotypes 06_heatmap 07_federation \
       08_project_healthy 09_project_timepoints 10_federated_modules \
       11_modules_per_cohort 12_panels_per_cohort 13_cluster_both_axes \
       14_heatmaps 15_kmeans_arm 16_project_and_federate \
       17_cohort_diagnostics 18_federate_per_module 19_federation_benefit \
       20_module_preservation)
# --optional is accepted and ignored. Steps 08 and 09 used to sit behind it and
# were therefore skipped by a plain run, which left their committed notebooks
# with no output at all. They take 5 and 4 seconds; the flag was not worth it.
[[ "${1:-}" == "--optional" ]] && true

# 02b is the slow one, about 10 minutes: it runs one ReactomePA over-representation
# test per module, 52 of them, and enrichPathway rebuilds its Reactome mapping on
# every call. It is not cached, because a cached enrichment would outlive a change
# to the fit or to reactome.db without saying so.

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

#!/bin/bash
#SBATCH --job-name=emdq_data
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%A_%a.log
# Data-matched analysis (jetana_data) of the kept event files, one file per
# array element (about 15 s each), e.g.
#   sbatch --account=<account> --array=1-$(wc -l < tasks/bestfit.txt)%50 \
#          --export=ALL,TASKS=tasks/bestfit.txt slurm/data_array.sbatch
# PIONMASS=1 in --export gives every charged particle the pion mass, as in
# ALICE's track-based jets, and writes spec_data_pimass/ (used for ALICE).
# IMAGE/BIND as in array.sbatch.
cd "$SLURM_SUBMIT_DIR" || exit 1
T=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${TASKS:?set TASKS=tasks/<group>.txt}")
[ -z "$T" ] && exit 0
O=spec_data; [ "$PIONMASS" = 1 ] && O=spec_data_pimass
mkdir -p $O
if [ -n "$IMAGE" ]; then
  apptainer exec --env JETANA_PIONMASS=${PIONMASS:-0} ${BIND:+-B "$BIND"} "$IMAGE" ./jetana_data $O/$T.spec events/$T.dat.gz
else
  JETANA_PIONMASS=${PIONMASS:-0} ./jetana_data $O/$T.spec events/$T.dat.gz
fi

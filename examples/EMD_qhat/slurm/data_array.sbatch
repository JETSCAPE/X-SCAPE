#!/bin/bash
#SBATCH --job-name=emdq_data
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%A_%a.log
# Data-matched analysis (jetana_data) of the kept event files, one file per
# array element (about 15 s each), e.g.
#   sbatch --account=<account> --array=1-$(wc -l < tasks/bestfit.txt)%150 \
#          --export=ALL,TASKS=tasks/bestfit.txt slurm/data_array.sbatch
# IMAGE/BIND as in array.sbatch.
cd "$SLURM_SUBMIT_DIR" || exit 1
T=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${TASKS:?set TASKS=tasks/<group>.txt}")
[ -z "$T" ] && exit 0
mkdir -p spec_data
if [ -n "$IMAGE" ]; then
  apptainer exec ${BIND:+-B "$BIND"} "$IMAGE" ./jetana_data spec_data/$T.spec events/$T.dat.gz
else
  ./jetana_data spec_data/$T.spec events/$T.dat.gz
fi

#!/bin/bash
#SBATCH --job-name=emdq
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.log
# One hydro task per array element. Submit from this folder, e.g.
#   sbatch --account=<account> --array=1-$(wc -l < tasks/routeB.txt)%80 \
#          --export=ALL,TASKS=tasks/routeB.txt slurm/array.sbatch
# LHC tasks take about 1.5 h, RHIC tasks about 0.5 h (2 threads).
# Add KEEP_EVENTS=1 to --export to keep the final-state hadrons (events/).
# Optional container: IMAGE=<jetscape image .sif> and BIND=<host paths to
# mount, comma-separated, covering this folder and the X-SCAPE build>.
# Do not rebuild the X-SCAPE build used here while tasks run: a running
# runJetscape that loses its shared libraries dies with SIGBUS.
cd "$SLURM_SUBMIT_DIR" || exit 1
T=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${TASKS:?set TASKS=tasks/<group>.txt}")
[ -z "$T" ] && exit 0
if [ -n "$IMAGE" ]; then
  apptainer exec ${BIND:+-B "$BIND"} "$IMAGE" bash run_task.sh "$T" 2
else
  bash run_task.sh "$T" 2
fi

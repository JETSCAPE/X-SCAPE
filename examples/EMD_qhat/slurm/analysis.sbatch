#!/bin/bash
#SBATCH --job-name=emdq_ana
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.log
# Run one analysis script on a compute node, e.g.
#   sbatch --account=<account> --export=ALL,SCRIPT=option1_analysis.py slurm/analysis.sbatch
# IMAGE/BIND as in array.sbatch (the scripts need numpy and matplotlib).
cd "$SLURM_SUBMIT_DIR" || exit 1
if [ -n "$IMAGE" ]; then
  apptainer exec ${BIND:+-B "$BIND"} "$IMAGE" python3 -W ignore ${SCRIPT:?set SCRIPT}
else
  python3 -W ignore ${SCRIPT:?set SCRIPT}
fi

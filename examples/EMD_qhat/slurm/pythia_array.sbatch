#!/bin/bash
#SBATCH --job-name=emdq_pythia
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%A_%a.log
# Standalone PYTHIA 8 pp reference (pythia_pp.cc), 25k events per element
# (about 5 min), then jetana_data with and without the pion mass, e.g.
#   sbatch --account=<account> --array=1-16%50 slurm/pythia_array.sbatch
# Lines of slurm/pythia_tasks.txt: <name> <nEvents> <seed> <mpi 1|0>.
# IMAGE/BIND as in array.sbatch.
cd "$SLURM_SUBMIT_DIR" || exit 1
read n nev seed mpi <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" slurm/pythia_tasks.txt)"
[ -z "$n" ] && exit 0
mkdir -p events spec_data spec_data_pimass
run() { if [ -n "$IMAGE" ]; then apptainer exec ${BIND:+-B "$BIND"} "$IMAGE" "$@"; else "$@"; fi; }
run ./pythia_pp events/$n.dat.gz $nev $seed $mpi > logs/$n.pythia.log 2>&1 || exit 1
run ./jetana_data spec_data/$n.spec events/$n.dat.gz
run env JETANA_PIONMASS=1 ./jetana_data spec_data_pimass/$n.spec events/$n.dat.gz

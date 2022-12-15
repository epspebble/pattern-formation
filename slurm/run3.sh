#!/bin/bash
#SBATCH --account=def-simontse
#SBATCH --cpus-per-task=1        
#SBATCH --mem-per-cpu=32G
#SBATCH --time=00:12:00           # time (DD-HH:MM)
#SBATCH --array=1-45
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simon.tse@twu.ca

module load julia/1.8.1

FN="~/project/pattern-formation/slurm/run3.par"
CMD="sed -n '${SLURM_ARRAY_TASK_ID}p' ${FN}"

cd ~/project/GrayScott/nc
julia ~/project/pattern-formation/gs_gen.jl 50000 $(eval "$CMD")

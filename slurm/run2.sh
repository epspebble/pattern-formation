#!/bin/bash
#SBATCH --account=def-simontse
#SBATCH --cpus-per-task=1        
#SBATCH --mem-per-cpu=64G
#SBATCH --time=1-00:00           # time (DD-HH:MM)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simon.tse@twu.ca

module load julia/1.8.1
cd ~/project/GrayScott/nc
julia ~/project/pattern-formation/gs_gen.jl 100000 0.062 0.0609

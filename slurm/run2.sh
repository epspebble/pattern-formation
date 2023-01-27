#!/bin/bash
#SBATCH --account=def-simontse
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-00:20
#SBATCH --mail-user=simon.tse@twu.ca
#SBATCH --mail-type=ALL
#SBATCH --profile=all
#SBATCH -e %x.e%A-%a
#SBATCH -o %x.o%A-%a
#SBATCH --array=1-231
##SBATCH --mail-type=FAIL
##SBATCH --mail-type=BEGIN,END

module load matlab
cd ~/project/pattern-formation/zebrafish
matlab -batch "run2($SLURM_ARRAY_TASK_ID)"

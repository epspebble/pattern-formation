#!/bin/bash
#SBATCH --account=def-simontse
#SBATCH --time=0-12:00
#SBATCH --mem=1G
#SBATCH --mail-user=simon.tse@twu.ca
##SBATCH --mail-type=FAIL
##SBATCH --mail-type=BEGIN,END
#SBATCH --mail-type=ALL
#SBATCH -e %x.e%A-%a
#SBATCH -o %x.o%A-%a
#SBATCH --profile=all

RUN_NUM=3
RUN_DIR=run$RUN_NUM
PROG=jt_run$RUN_NUM

cd ~/project/cytoneme
cp ./run_template/juxtacrine2d_v2.f ./$RUN_DIR
cp ./run_template/params_bkp.inc ./$RUN_DIR

cd ./$RUN_DIR
sed s/nx=50,ny=50/nx=75,ny=75/ <params_bkp.inc >params.inc
gfortran juxtacrine2d_v2.f -w -o $PROG 

./$PROG

#for ((i=$1;i<$2;++i)); 
#do
#    let i0=$i*200;
#    let i1=$i0+200;
#    ~/scripts/batch_prepare.py 144x $i0 $i1;
#done


# Job ID: 48709971
# Cluster: cedar
# User/Group: simontse/simontse
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 01:00:37
# CPU Efficiency: 99.70% of 01:00:48 core-walltime
# Job Wall-clock time: 01:00:48
# Memory Utilized: 727.81 MB
# Memory Efficiency: 71.08% of 1.00 GB

#!/bin/csh
#PBS -l nodes=1:ppn=1
#PBS -q workq
#PBS -l walltime=48:00:00
#PBS -N runWM
cd $PBS_O_WORKDIR
module load python/2.7.6
python WMrateget.py 34

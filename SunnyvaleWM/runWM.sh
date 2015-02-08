#!/bin/csh
#PBS -l nodes=3:ppn=8
#PBS -q workq
#PBS -l walltime=08:00:00
#PBS -N runWM
cd $PBS_O_WORKDIR
module load python/2.7.6
python WMrateget1.py & python WMrateget2.py & python WMrateget3.py & python WMrateget4.py & python WMrateget5.py & python WMrateget6.py & python WMrateget7.py & python WMrateget8.py & python WMrateget9.py & python WMrateget10.py & python WMrateget11.py & python WMrateget12.py & python WMrateget13.py & python WMrateget14.py & python WMrateget15.py & python WMrateget16.py & python WMrateget17.py & python WMrateget18.py & python WMrateget19.py & python WMrateget20.py & python WMrateget21.py & python WMrateget22.py & python WMrateget23.py 

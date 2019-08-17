#!/bin/bash 
  
#PBS -P Project_Name_of_Job  
#PBS -j oe 
#PBS -N Job_Name 

#PBS -q parallel12 
#PBS -l select=1:ncpus=12:mpiprocs=12:mem=40GB 

###--- Setting for parallel8, parallel12, parallel24, parallel20 ---- 
###--- To activate the setting, change ###PBS to #PBS for the line. ---- 
###PBS -q parallel8 
###PBS -l select=1:ncpus=8:mpiprocs=8:mem=40GB 

###PBS -q parallel12 
###PBS -l select=1:ncpus=12:mpiprocs=12:mem=40GB 

###PBS -q parallel24 
###PBS -l select=1:ncpus=24:mpiprocs=24:mem=160GB 

###PBS -q parallel20 
###PBS -l select=2:ncpus=20:mpiprocs=20:mem=160GB -l place=pack 

###--- Set different walltime/running time (max=720 hours) ---- 
###--- To activate the setting, change ###PBS to #PBS for the line. ---- 
###PBS -l walltime=720:00:00 

cd $PBS_O_WORKDIR;   ## this line is needed, do not delete and change.
np=$( cat  ${PBS_NODEFILE} |wc -l );  ### get number of CPUs, do not change 

##--- Put your exec/application commands below ---  
##--- For example: 
source /etc/profile.d/rec_modules.sh 
module load xe_2015;  module list 

mpirun -f ${PBS_NODEFILE} ./a.out 

##--- END HERE --- 

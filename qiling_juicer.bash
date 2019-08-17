#!/bin/bash 
  
#PBS -P qilling_hic  
#PBS -j oe 
#PBS -N juicer 

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
bash /hpctmp/e0056363/juicer/scripts/juicer.sh \
-y /root/qiling/juicer/restriction_sites/mm10_noScaffold_Arima.txt \
-z /root/qiling/juicer/references/mm10_noScaffold.fasta \
-p /hpctmp/e0056363/juicer/references/mm10.chrom.sizes_noScaffold &> juicer_run.log

##--- END HERE --- 

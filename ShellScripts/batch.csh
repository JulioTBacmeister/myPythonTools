#!/bin/csh
#PBS -N timesh

### Charging account
#PBS -A P93300642 
### Request one chunk of resources with 1 CPU and 10 GB of memory
#PBS -l select=1:ncpus=1:mem=10GB
### Allow job to run up to 30 minutes
#PBS -l walltime=11:00:00
### Route the job to the casper queue
#PBS -q casper
### Join output and error streams into single file
#PBS -j oe



# This job's working directory
#
echo `date`
if ( $?PBS_O_WORKDIR == 1 ) then
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
endif


module load conda

conda activate npl

python data_for_sam_ml.py

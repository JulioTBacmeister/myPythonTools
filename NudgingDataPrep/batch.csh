#!/bin/csh
#PBS -N ERA5proc

### Charging account
#PBS -A P93300642 
### Request one chunk of resources with N CPU and M GB of memory
#PBS -l select=1:ncpus=16:mem=256GB
### 
#PBS -l walltime=2:00:00
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

conda activate npl-2022b

./GenRegrid.py --year=2021 --month=6 --day=1 --hour=99

#!/bin/csh
#PBS -N prectrx

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

module load conda

conda activate npl

cd /glade/work/juliob/myPythonTools/TropicalCyclones

./create_PrecipExceedance.py --CAM=1 --sstA='sst2' --year0=2070 --year1=2099

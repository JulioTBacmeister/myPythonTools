#!/bin/csh
#PBS -N gwdrive

### Charging account
#PBS -A P93300642 
### Request one chunk of resources with N CPU and M GB of memory
#PBS -l select=1:ncpus=1:mem=128GB
### 
#PBS -l walltime=2:00:00
### Route the job to the casper queue
#PBS -q casper
### Join output and error streams into single file
#PBS -j oe



# This job's working directory
#


module load conda

conda activate npl-2023b
#./gw_driver.py --year=2011 --month=1 --SourceMethod=uniform
#./gw_driver.py --year=2011 --month=1 --SourceMethod=vort500
./gw_driver.py --year=2010 --month=10 --SourceMethod=SavedCLUBBmomflux
#./gw_driver.py --year=2010 --month=10 --SourceMethod=vort500

#!/bin/csh
#PBS -N gwanana

### Charging account
#PBS -A P93300642 
### Request one chunk of resources with N CPU and M GB of memory
#PBS -l select=1:ncpus=1:mem=64GB
### 
#PBS -l walltime=12:00:00
### Route the job to the casper queue
#PBS -q casper
### Join output and error streams into single file
#PBS -j oe



# This job's working directory
#


module load conda

conda activate npl-2023b

./zonal_avg_wrt.py


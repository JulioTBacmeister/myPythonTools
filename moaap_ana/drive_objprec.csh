#!/bin/csh
#PBS -N moaap

### Charging account
#PBS -A P93300042 
### Request one chunk of resources with N CPU and M GB of memory
#PBS -l select=1:ncpus=1:mem=128GB
### 
#PBS -l walltime=06:00:00
### Route the job to the casper queue
#PBS -q casper
### Join output and error streams into single file
#PBS -j oe


./object_precip_proc.py --year 2005

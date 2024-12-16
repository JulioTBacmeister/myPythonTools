#!/bin/csh
#PBS -N gwanana

### Charging account
#PBS -A P93300042 
### Request one chunk of resources with N CPU and M GB of memory
#PBS -l select=1:ncpus=4:mem=64GB
### 
#PBS -l walltime=6:00:00
### Route the job to the casper queue
#PBS -q casper
### Join output and error streams into single file
#PBS -j oe



# This job's working directory
#


module load conda

conda activate npl-2024b


#./regrid_gwana_wrt.py
#./regrid_gwana_latSlice_wrt.py
#./phs2_xy_gwana.py

./regrid_vorticity_wrt.py
./regrid_PBLq_wrt.py
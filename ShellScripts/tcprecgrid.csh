#!/bin/csh
#PBS -N tcprgrid

### Charging account
#PBS -A P93300642 
### Request one chunk of resources with 1 CPU and 10 GB of memory
#PBS -l select=1:ncpus=1:mem=10GB
### Allow job to run up to 11 hours
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

./create_TRMM_TCprecip_grid.py 
./create_TCprecip_grid.py --ens=1 --year0=1979 --year1=2012
./create_TCprecip_grid.py --ens=2 --year0=1979 --year1=2012
./create_TCprecip_grid.py --ens=3 --year0=1979 --year1=2012
./create_TCprecip_grid.py --sstA='sst1' --year0=2070 --year1=2099
./create_TCprecip_grid.py --sstA='sst2' --year0=2070 --year1=2099
./create_TCprecip_grid.py --sstA='sst7' --year0=2070 --year1=2099
./create_TCprecip_grid.py --sstA='sst6' --year0=2070 --year1=2099


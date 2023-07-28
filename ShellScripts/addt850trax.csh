#!/bin/csh
#PBS -N adt1850trk

### Charging account
#PBS -A P93300642 
### Request one chunk of resources with 1 CPU and 10 GB of memory
#PBS -l select=1:ncpus=1:mem=64GB
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

./add_t850_trax.py --sstA='sst1' --year0=2070 --year1=2099
./add_t850_trax.py --sstA='sst2' --year0=2070 --year1=2099
./add_t850_trax.py --ens=1 --year0=1979 --year1=2012
./add_t850_trax.py --ens=2 --year0=1979 --year1=2012
./add_t850_trax.py --ens=3 --year0=1979 --year1=2012
./add_t850_trax.py --sstA='sst1' --ens=2 --year0=2070 --year1=2099
./add_t850_trax.py --sstA='sst1' --ens=3 --year0=2070 --year1=2099



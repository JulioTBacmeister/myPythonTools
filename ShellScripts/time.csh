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

# Submit like this:            
#   qsub -N poo -j oe -A P93300642 -q casper -l select=1:ncpus=1:mem=10GB,walltime=08:00:00 -- /glade/work/juliob/SAMwrf_grids/regrid.csh 2015 6



# This job's working directory
#
echo `date`
if ( $?PBS_O_WORKDIR == 1 ) then
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
endif


module load conda

conda activate npl

#./drv_time_management.py --year=${year} --month=${month} 
#./drv_time_management.py --year=2010 --month=9
#./drv_time_management.py --year=2010 --month=10
#./drv_time_management.py --year=2010 --month=11

./regrid_vertical.py --year=2010 --month=9
./regrid_vertical.py --year=2010 --month=10
./regrid_vertical.py --year=2010 --month=11


./drv_time_management.py --year=2010 --month=12
./regrid_vertical.py --year=2010 --month=12

./drv_time_management.py --year=2011 --month=1
./regrid_vertical.py --year=2011 --month=1

./drv_time_management.py --year=2011 --month=2
./regrid_vertical.py --year=2011 --month=2

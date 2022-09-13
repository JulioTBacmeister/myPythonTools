#!/bin/csh
#PBS -N zonav
#PBS -r n
#PBS -e zonav.err
#PBS -o zonav.log
#PBS -m ae
#PBS -q long
#
# Number of nodes, number of processors
#
# nodes = physical host
# ppn   = processors per node (i.e., number of cores)
#
#PBS -l nodes=1:ppn=48

# Submit like this:            
#   qsub -F "2010 4" zonav.sh

if ( "$#argv" != 2) then
  echo "Wrong number of arguments specified:"
  echo "  -arg 1 year"
  echo "  -arg 2 month"
  exit
endif

set n = 1
set year = "$argv[$n]"
set n = 2
set month = "$argv[$n]"

#
# This job's working directory
#
echo `date`
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

conda activate adf_v0.07

./drv_hf_zonal_mean.py --year=${year} --month=${month}

exit 0



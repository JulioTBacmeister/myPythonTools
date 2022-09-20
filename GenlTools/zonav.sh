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
#   qsub -F "2010 6 c6_3_59.f09_L58.CTL01" zonav.sh

if ( "$#argv" != 2) then
  echo "Wrong number of arguments specified:"
  echo "  -arg 1 year"
  echo "  -arg 2 month"
endif

set lonfill = 0
set n = 1
set year = "$argv[$n]"
set n = 2
set month = "$argv[$n]"
set n = 3
set case = "$argv[$n]"


if ( "$#argv" == 4) then
   set lonfill = 1
endif
#
# This job's working directory
#
echo `date`
if ( $?PBS_O_WORKDIR == 1 ) then
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
endif

conda activate adf_v0.07

if ( $lonfill == 1 ) then
echo " Will fill Longitudes for 3D output "
./drv_hf_zonal_mean.py --year=${year} --month=${month} --case=${case} -F
else
echo " 2D Y-Z output "
./drv_hf_zonal_mean.py --year=${year} --month=${month} --case=${case}
endif


#./drv_hf_zonal_yz.py --year=${year} --month=${month}

exit 0



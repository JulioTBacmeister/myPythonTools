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
#   qsub -F "2010 6 c6_3_59.f09_L58.CTL01" addphis.csh

set n = 1
set year = "$argv[$n]"
set n = 2
set month = "$argv[$n]"
set n = 3
set case = "$argv[$n]"


#
# This job's working directory
#
echo `date`
if ( $?PBS_O_WORKDIR == 1 ) then
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
endif

conda activate adf_v0.07

./drv_hf_add_phis.py --year=${year} --month=${month} --case=${case}

exit 0



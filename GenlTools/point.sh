#!/bin/csh
#PBS -N pointsh
#PBS -r n
#PBS -e pointsh.err
#PBS -o pointsh.log
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
#   qsub -F point.sh

# This job's working directory
#
echo `date`
if ( $?PBS_O_WORKDIR == 1 ) then
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
endif

conda activate adf_v0.07

./drv_hf_zonal_mean.py -y 2010 -m 2 -P 262.5,36.6 -X c6_3_59.f09_L58.CTL02
./drv_hf_zonal_mean.py -y 2010 -m 3 -P 262.5,36.6 -X c6_3_59.f09_L58.CTL02
./drv_hf_zonal_mean.py -y 2010 -m 4 -P 262.5,36.6 -X c6_3_59.f09_L58.CTL02
./drv_hf_zonal_mean.py -y 2010 -m 5 -P 262.5,36.6 -X c6_3_59.f09_L58.CTL02
./drv_hf_zonal_mean.py -y 2010 -m 6 -P 262.5,36.6 -X c6_3_59.f09_L58.CTL02
./drv_hf_zonal_mean.py -y 2010 -m 7 -P 262.5,36.6 -X c6_3_59.f09_L58.CTL02
./drv_hf_zonal_mean.py -y 2010 -m 8 -P 262.5,36.6 -X c6_3_59.f09_L58.CTL02
./drv_hf_zonal_mean.py -y 2010 -m 9 -P 262.5,36.6 -X c6_3_59.f09_L58.CTL02
./drv_hf_zonal_mean.py -y 2010 -m 10 -P 262.5,36.6 -X c6_3_59.f09_L58.CTL02
./drv_hf_zonal_mean.py -y 2010 -m 11 -P 262.5,36.6 -X c6_3_59.f09_L58.CTL02
./drv_hf_zonal_mean.py -y 2010 -m 12 -P 262.5,36.6 -X c6_3_59.f09_L58.CTL02

#./drv_hf_zonal_yz.py --year=${year} --month=${month}

exit 0



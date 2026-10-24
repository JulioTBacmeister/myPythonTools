#!/bin/tcsh
#PBS -N concat_monthly_to_bimonthly
#PBS -A CESM0024 
#PBS -q casper
#PBS -l select=1:ncpus=1:mem=300GB
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -o logs/

# Go to submission directory
cd $PBS_O_WORKDIR

# Make sure logs directory exists
mkdir -p logs

# Initialize modules and load NCO
if ( -f /etc/profile.d/modules.csh ) then
    source /etc/profile.d/modules.csh
endif
module purge
module load nco

echo "Using ncrcat: `which ncrcat`"
echo "`date`: Starting concatenation"

# Paths
set indir  = /glade/derecho/scratch/hannay/MOAPP/monthly_files
set outdir = /glade/derecho/scratch/hannay/MOAPP/bimonthly_files

mkdir -p $outdir

# Files
set file1 = $indir/CAM7_L58_HR_hourly.199001.nc
set file2 = $indir/CAM7_L58_HR_hourly.199002.nc
set outfile = $outdir/CAM7_L58_HR_hourly.199001.nc

echo "Input 1: $file1"
echo "Input 2: $file2"
echo "Output : $outfile"

# Run concatenation
##ncrcat -h -O $file1 $file2 $outfile
set status = $status



echo "`date`: ncrcat finished with status $status"
exit $status
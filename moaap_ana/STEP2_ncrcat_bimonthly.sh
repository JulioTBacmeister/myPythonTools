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

# Loop over years
foreach year (`seq 1990 2010`)

    # Loop over months 01,03,05,...,11
    foreach month (01 03 05 07 09 11)

        # Compute the second month of the pair
        @ m2 = $month + 1
        if ($m2 < 10) then
            set month2 = 0$m2
        else
            set month2 = $m2
        endif

        # Input files
        set file1   = $indir/CAM7_L58_HR_hourly.${year}${month}.nc
        set file2   = $indir/CAM7_L58_HR_hourly.${year}${month2}.nc

        # Output file
        set outfile = $outdir/CAM7_L58_HR_hourly.${year}${month}.nc

        # Check existence
        if (! -e $file1) then
            echo "MISSING: $file1 — skipping"
            continue
        endif
        if (! -e $file2) then
            echo "MISSING: $file2 — skipping"
            continue
        endif

        # Skip if output already exists
        if (-e $outfile) then
            echo "Exists → skipping $outfile"
            continue
        endif

        echo "Concatenating:"
        echo "   $file1"
        echo " + $file2"
        echo " → $outfile"

        ncrcat -h -O $file1 $file2 $outfile
        set status = $status
        echo "Status = $status"
        echo "-------------------------"

    end
end
#!/bin/tcsh
#PBS -N concat_hourly_to_monthly
#PBS -A CESM0024      
#PBS -q casper                 
#PBS -l select=1:ncpus=1:mem=100GB
#PBS -l walltime=12:00:00
#PBS -o logs/
#PBS -j oe                   


# Go to the directory where you submit the job
cd $PBS_O_WORKDIR

# Load NCO
module load nco

# Directories
set daily_dir   = /glade/derecho/scratch/hannay/MOAPP/daily_files
set monthly_dir = /glade/derecho/scratch/hannay/MOAPP/monthly_files

mkdir -p $monthly_dir

echo "`date` : Starting monthly concatenation on casper"
echo "  daily_dir   = $daily_dir"
echo "  monthly_dir = $monthly_dir"
echo ""

# Loop over years 1990–2010
@ year = 1990
while ( $year <= 2010 )

  # Loop over months (zero-padded)
  foreach month (01 02 03 04 05 06 07 08 09 10 11 12)

    set pattern = ${daily_dir}/CAM7_L58_HR_hourly.${year}-${month}-??.nc
    set outfile = ${monthly_dir}/CAM7_L58_HR_hourly.${year}${month}.nc

    # If output already exists, skip
    if ( -e $outfile ) then
      echo "`date` : Skipping ${year}-${month} (output exists: $outfile)"
      continue
    endif

    # Expand pattern into a list of files
    set files = ( $pattern )

    # If no files found, skip
    if ( $#files == 0 ) then
      echo "`date` : No daily files for ${year}-${month}, skipping"
      continue
    endif

    echo ""
    echo "`date` : Concatenating ${year}-${month}"
    echo "  Input files example: $files[1]"
    echo "  Output file        : $outfile"
    echo ""

    # Concatenate along record dim (time)
    ncrcat -h -O $files $outfile

    if ( $status != 0 ) then
      echo "`date` : ERROR during ncrcat for ${year}-${month}, status = $status"
      echo "           Check files and output; stopping job."
      exit 1
    else
      echo "`date` : Finished ${year}-${month}"
    endif

  end

  @ year++
end

echo ""
echo "`date` : All done."

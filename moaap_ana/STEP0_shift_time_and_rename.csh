#!/bin/csh

# renaming files because the time units between 001 and 002 is different.
# - change time variable: time = time + 4015
# - and change units time:units = "days since 2001-01-01 00:00:00"  
#                to  time:units = "days since 1990-01-01 00:00:00" ;

# to make it safe, I make a copy of the 002 files into: 
# glade/campaign/cesm/community/amwg/cam_hr/f.e30_cam6_4_116.FHISTC_LTso.ne120_L58.c6aero.001/atm/hist/002 
# and make my changes over there. 

# There is some overlap for the year 2001 between Julio and Cecile's files because 
# - Julio ran until mid 2001
# - Julio's files from December 2000 can have timesteps into 2001
# so I created a directory with Julio original files: 
# /glade/campaign/cesm/community/amwg/cam_hr/f.e30_cam6_4_116.FHISTC_LTso.ne120_L58.c6aero.001/atm/hist/julio_orig_files

module load nco   
cd /glade/campaign/cesm/community/amwg/cam_hr/f.e30_cam6_4_116.FHISTC_LTso.ne120_L58.c6aero.001/atm/hist/002

# Loop through all c6aero.002 files
foreach f002 ( f.e30_cam6_4_116.FHISTC_LTso.ne120_L58.c6aero.002.cam.h*a*.nc )

    # Create new filename by replacing c6aero.002 with c6aero.001
    set f001 = `echo $f002 | sed 's/c6aero\.002/c6aero\.001/'`

    # Skip file if the output already exists
    if ( -e $f001 ) then
        echo "Output file $f001 already exists."
        continue
    endif
    
    echo "Processing: $f002 to $f001"
    ncap2 -O \
    -s 'time = time + 4015; time@units = "days since 1990-01-01 00:00:00"; time_bounds = time_bounds + 4015; time_bounds@units = "days since 1990-01-01 00:00:00"' \
    $f002 $f001

 end
 
foreach f002 ( f.e30_cam6_4_116.FHISTC_LTso.ne120_L58.c6aero.002.cam.h*i*.nc )

    # Create new filename by replacing c6aero.002 with c6aero.001
    set f001 = `echo $f002 | sed 's/c6aero\.002/c6aero\.001/'`

    # Skip file if the output already exists
    if ( -e $f001 ) then
        echo "Output file $f001 already exists."
        continue
    endif
    
    echo "Processing: $f002 to $f001"
    ncap2 -O \
    -s 'time = time + 4015; time@units = "days since 1990-01-01 00:00:00"' \
    $f002 $f001

 end
 

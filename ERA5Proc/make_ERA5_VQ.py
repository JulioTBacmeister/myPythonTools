################################################
# New style 
################################################
import sys
import os

workdir_ = '/glade/work/juliob'
if ( workdir_ not in sys.path ):
    sys.path.append(workdir_)
    print( f" a path to {workdir_} added in {__name__} ")

import importlib
import argparse as arg
import time
import numpy as np
import xarray as xr

# import modules in other directories
from PyRegridding.Regridder import GenRegrid as GnR
from PyRegridding.Regridder import Initialize as Ini
from PyRegridding.Regridder import WriteDST as Wrt
from PyRegridding.Regridder import ReadInSrc as Rd 

from myPythonTools.Utils import utils as uti



import WriteVQ as wVQ

# The usual
from datetime import date

# Some other useful packages 
import copy
import time
import cftime

importlib.reload( Rd )

def main(year, month ):
    import calendar
    
    tic_total = time.perf_counter()
    days_in_month = calendar.monthrange(year,month)[1]

    print( f"About to process hourly VQ {year:04d}-{month:02d}",flush=True )

    rc=Ini.prep(Dst = 'fv1x1', DstVgrid='L93',  Src='ERA5', WOsrf=False , RegridMethod=None , IC_for_pg=False )

    ndays = uti.days_in_month(year, month, check_for_leap_year=True)
    
    for d in np.arange( start=1, stop=ndays+1):
        for h in np.arange( start=0,stop=24,step=6):
        
            rc=Rd.get_ERA5_just_VQ( year=year, month=month, day=d, hour0=h)
    
            rc=wVQ.write_netcdf(year=year, month=month)
    
    #rc=wVQ.write_monthly_netcdf(year=year, month=month)
    
    
    code = 1
    toc_total = time.perf_counter()

    pTime = f"Total processing time was  {toc_total - tic_total:0.4f} seconds"
    print(pTime,flush=True)


def daily(year, month ):
    import calendar
    
    tic_total = time.perf_counter()
    print( f"About to make daily averages for {year:04d}-{month:02d}",flush=True )
    ndays = uti.days_in_month(year, month, check_for_leap_year=True)
    
    for d in np.arange( start=1, stop=ndays+1):
        rc = wVQ.write_daily_netcdf(year=year,month=month,day=d,return_dataset=False)
        
    code = 1
    toc_total = time.perf_counter()

    pTime = f"Total processing time was  {toc_total - tic_total:0.4f} seconds"
    print(pTime,flush=True)

def monthly(year, month ):
    import calendar
    
    tic_total = time.perf_counter()
    print( f"About to make monthly average for {year:04d}-{month:02d}",flush=True )

    rc = wVQ.write_monthly_netcdf(year=year,month=month,return_dataset=False)
    
    code = 1
    toc_total = time.perf_counter()

    pTime = f"Total processing time was  {toc_total - tic_total:0.4f} seconds"
    print(pTime,flush=True)



if __name__ == "__main__":
    
    my_parser = arg.ArgumentParser()
    my_parser.add_argument("--month",    type=int, default=1)
    my_parser.add_argument("--year",     type=int, default=2010)
    args = my_parser.parse_args()
    main(args.year, args.month )


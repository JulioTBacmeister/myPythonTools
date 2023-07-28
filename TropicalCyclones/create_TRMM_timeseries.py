#!/usr/bin/env python

# Import packages 
import sys
import argparse as arg

import time
import xarray as xr
import numpy as np
import pandas as pd

def main(year):
    
    tic_total = time.perf_counter()

    print( f"About to process Year {year} ")
    yearA = str(year).zfill(4)
    
    files = '/glade/p/cgd/amp/juliob/TRMM/HDF/' + yearA + '/3B42.' + yearA + '*.nc'
    ofile = '/glade/p/cgd/amp/juliob/TRMM/HDF/' + yearA + '/3B42_3hrly.' + yearA + '.nc'

    dS=xr.open_mfdataset( files )
    dS.to_netcdf( ofile )
    
    toc_total = time.perf_counter()
    pTime = f"Total processing time was  {toc_total - tic_total:0.4f} seconds"
    print(pTime)

if __name__ == "__main__":
    # argument: indir -> get all nc files in this directory
    # argument: map -> the offlinemap file already prepared
    # argument: outdir -> directory where remapped files should go
    # my_parser = arg.ArgumentParser()
    # my_parser.add_argument("--month", type=int)
    # my_parser.add_argument("--year", type=int)
    # args = my_parser.parse_args()
    
    my_parser = arg.ArgumentParser()
    my_parser.add_argument("--year",     type=int, default=2010)
    args = my_parser.parse_args()
    main( args.year )

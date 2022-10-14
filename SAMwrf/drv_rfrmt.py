#!/usr/bin/env python

import xarray as xr
import os
import subprocess as sp
import sys
import numpy as np
import argparse as arg
from pathlib import Path


def main(year,month):

    input_dir = "/glade/p/cesm/amwg_dev/juliob/SAMwrf/ne30x16/2015/nc4/"
    nc3_dir = "/glade/p/cesm/amwg_dev/juliob/SAMwrf/ne30x16/2015/nc3/"
    cdf5_dir = "/glade/p/cesm/amwg_dev/juliob/SAMwrf/ne30x16/2015/cdf5/"
    work_dir = "/glade/work/juliob/SAMwrf_grids/"

    ifile_root = "f.e22r.SAMwrf01.f09.L70.NODEEP_2015_01.cam.h1."
    ofile_root = "f.e22r.SAMwrf01.f09.L70.NODEEP_2015_01.cam.h1."




    days_in_month =[31 , 28, 31, 30, 31, 30, 31, 31, 30, 31,30, 31 ]
    
    iy=year #2016
    imm=month # 6
    idd=1
    iss=3600

    print(" Days in month=",days_in_month[imm-1])


    ndays = days_in_month[imm-1]
    for idd in np.arange(2,ndays+1):
        for ihh in np.arange(0,24):
        
            yy=str(iy).zfill(4)
            mm=str(imm).zfill(2)
            dd=str(idd).zfill(2)
            ss=str(ihh*3600).zfill(5)
            tstamp=yy+'-'+mm+'-'+dd+'-'+ss+'.nc'
            
            ifile=input_dir+ifile_root+tstamp
            nc3file=nc3_dir+ifile_root+tstamp
 
            print("hour =",ihh,ifile)
            a=xr.open_dataset(ifile)
            #b=a
            a.to_netcdf( nc3file, format="NETCDF3_CLASSIC" )


if __name__ == "__main__":
    # argument: indir -> get all nc files in this directory
    # argument: map -> the offlinemap file already prepared
    # argument: outdir -> directory where remapped files should go
    my_parser = arg.ArgumentParser()
    my_parser.add_argument("--month", type=int)
    my_parser.add_argument("--year", type=int)
    args = my_parser.parse_args()
    main(args.year, args.month )

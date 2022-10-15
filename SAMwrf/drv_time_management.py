#!/usr/bin/env python

import xarray as xr
import os
import subprocess as sp
import sys
import numpy as np
import argparse as arg
from pathlib import Path

def get_cam_time(tim):
    ymd = get_ymd(tim)
    sec = get_seconds_str(tim)
    return ["-".join([ymd[i].item(),s.item()]) for i,s in enumerate(sec)]

def get_ymd(tim):
    return tim.dt.strftime("%Y-%m-%d")

def get_seconds_str(tim):
    seconds = tim.dt.hour*3600 + tim.dt.minute*60 + tim.dt.second
    seconds_str = [f'{s.item():05d}' for s in seconds]
    return xr.DataArray(seconds_str, dims='time', coords={'time':tim})

def main(year,month):

    input_dir  = "/glade/p/cesm/amwg_dev/juliob/SAMwrf/ne30x16/i/"
    output_dir = "/glade/p/cesm/amwg_dev/juliob/SAMwrf/ne30x16/o/"
    work_dir = "/glade/work/juliob/SAMwrf_grids/"

    ifile_root = "f.e22r.SAMwrf01.ne30.L32.NODEEP_2010_01.cam.h1."
    ofile_root = "f.e22r.SAMwrf01.ne30.L32.NODEEP_2010_01.cam.h1."




    days_in_month =[31 , 28, 31, 30, 31, 30, 31, 31, 30, 31,30, 31 ]
    
    iy=year #2016
    imm=month # 6
    idd=1
    iss=3600

    print(" Days in month=",days_in_month[imm-1])


    ndays = days_in_month[imm-1]
    for idd in np.arange(1,ndays+1):
        
        yy=str(iy).zfill(4)
        mm=str(imm).zfill(2)
        dd=str(idd).zfill(2)
        ss=str(iss).zfill(5)
        tstamp=yy+'-'+mm+'-'+dd+'-'+ss+'.nc'
        
        ifile=input_dir+ifile_root+tstamp
 

        a=xr.open_dataset(ifile)

        tees = get_cam_time( a.time )
        time = a.time

        for i in np.arange(24):
            g=a.sel( time= time[i] )
            ofile=output_dir+ofile_root + tees[i] + '.nc'
            print(ofile)
            g.to_netcdf( ofile , format="NETCDF3_CLASSIC" )


if __name__ == "__main__":
    # argument: indir -> get all nc files in this directory
    # argument: map -> the offlinemap file already prepared
    # argument: outdir -> directory where remapped files should go
    my_parser = arg.ArgumentParser()
    my_parser.add_argument("--month", type=int)
    my_parser.add_argument("--year", type=int)
    args = my_parser.parse_args()
    main(args.year, args.month )

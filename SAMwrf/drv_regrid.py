#!/usr/bin/env python


import os
import subprocess as sp
import sys
import numpy as np
import argparse as arg
from pathlib import Path



def main(year,month):

    input_dir =  "/glade/campaign/cgd/projects/NCGD0051/ENSO_2010/L32/f.e22r.SAMwrf01.ne30x16.L32.NODEEP_2010_01/atm/hist/"
    output_dir = "/glade/p/cesm/amwg_dev/juliob/SAMwrf/ne30x16/i/"
    work_dir = "/glade/work/juliob/SAMwrf_grids/"

    ifile_root = "f.e22r.SAMwrf01.ne30x16.L32.NODEEP_2010_01.cam.h1."
    #ofile_root = "f.e22r.SAMwrf01.f09.L32.NODEEP_2010_01.cam.h1."
    ofile_root = "f.e22r.SAMwrf01.ne30.L32.NODEEP_2010_01.cam.h1."




    days_in_month =[31 , 28, 31, 30, 31, 30, 31, 31, 30, 31,30, 31 ]
    
    iy=year #2016
    imm=month # 6
    idd=1
    iss=3600

    print(" Days in month=",days_in_month[imm-1])


    ndays = days_in_month[imm-1]
    for idd in np.arange(1,ndays+1):
    #for idd in np.arange(1,5):
        
        yy=str(iy).zfill(4)
        mm=str(imm).zfill(2)
        dd=str(idd).zfill(2)
        ss=str(iss).zfill(5)
        tstamp=yy+'-'+mm+'-'+dd+'-'+ss+'.nc'
        
        ifile=input_dir+ifile_root+tstamp
        ofile=output_dir+ofile_root+tstamp

        #wgtsfile = work_dir+ "SAMwrf_ne30x16_TO_f09-cnsrv.nc"
        wgtsfile = work_dir+ "SAMwrf_ne30x16_TO_ne30pg3-cnsrv.nc"
        
        # create command
        cmd1='module load nco/4.7.9'
        cmd2='ncremap -4 -m '+ wgtsfile + ' '+ ifile + ' '+ ofile
        #cmd2='touch ' + ofile
        xcmd=cmd1+'; '+cmd2
        #print(cmd2)
        sp.run("source /etc/profile.d/modules.sh ; module load nco/4.7.9 ; "+ cmd2  , shell=True)
        #sp.run( cmd1 , shell=True)
        #sp.run(['ls' , '-l'])

if __name__ == "__main__":
    # argument: indir -> get all nc files in this directory
    # argument: map -> the offlinemap file already prepared
    # argument: outdir -> directory where remapped files should go
    my_parser = arg.ArgumentParser()
    my_parser.add_argument("--month", type=int)
    my_parser.add_argument("--year", type=int)
    args = my_parser.parse_args()
    main(args.year, args.month )

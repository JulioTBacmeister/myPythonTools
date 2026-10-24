#!/usr/bin/env python
################################################

import datetime
import pandas as pd
import subprocess
import matplotlib.path as mplPath
import sys
import calendar
from calendar import monthrange

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import re
import argparse


import ncutils as nc
import importlib

importlib.reload(nc)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Object analysis.")
    parser.add_argument("--year", type=int, default=-9999,  help="Year to do ")
    return parser.parse_args()


def drive():
    
    #year=2004
    #objn='MRG_Objects'

    # Parse command-line arguments
    args = parse_arguments()

    year = args.year
    
    f0=f'/glade/derecho/scratch/hannay/MOAPP/moaap_output/monthly/{year:04}*_CAM7_L58_HR_ObjectMasks__dt-1h_MOAAP-masks.nc'
    fL=f'/glade/derecho/scratch/hannay/MOAPP/moaap_output/monthly_ne30/{year:04}*_CAM7_L58_LR_ObjectMasks__dt-1h_MOAAP-masks.nc'
    fobs=f'/glade/derecho/scratch/fudan/pub/MOAAP_ERA5/{year:04}*_ERA5_ObjectMasks__dt-1h_MOAAP-masks.nc'
    
    
    
    
    Xc=xr.open_mfdataset( f0 ,data_vars='different', coords='different')
    Xo=xr.open_mfdataset( fobs ,data_vars='different', coords='different')
    XL=xr.open_mfdataset( fL ,data_vars='different', coords='different')




    do_vars=[
            'MCS_Tb_Objects', 
            'MRG_Objects', 
            'IGW_Objects', 
            'Kelvin_Objects', 
            'EIG0_Objects', 
            'ER_Objects',
            'FR_Objects', 
            'CY_Objects', 
            'ACY_Objects', 
            'MS_Objects', 
            'IVT_Objects', 
            'AR_Objects', 
            'CY_z500_Objects', 
            'ACY_z500_Objects', 
            'COL_Objects', 
            'JET_Objects', 
            'BT_Objects', 
            'TC_Objects', 
            ]

    latL=XL.lat.values
    lonL=XL.lon.values
    
    lato=Xo.lat.values
    lono=Xo.lon.values
    
    latc=Xc.lat.values
    lonc=Xc.lon.values
    
    
    ntL,nyL,nxL=XL.sizes['time'],XL.sizes['yc'],XL.sizes['xc']
    ndL=ntL//24

    nto,nyo,nxo=Xo.sizes['time'],Xo.sizes['yc'],Xo.sizes['xc']
    ndo=nto//24
    
    ntc,nyc,nxc=Xc.sizes['time'],Xc.sizes['yc'],Xc.sizes['xc']
    ndc=ntc//24
    
    for objn in do_vars:
    
        print( f" Doing {objn}" )

        objs_d_L= np.zeros( ( ndL, nyL, nxL ) )
        objprec_d_L= np.zeros( ( ndL, nyL, nxL ) )
        prec_d_L= np.zeros( ( ndL, nyL, nxL ) )
        
        for d in np.arange(ndc):
            # read in object masks and precip in one day chunks
            objs_L  = XL[objn][d*24:d*24+24,:,:].values
            prec    = XL['PR'][d*24:d*24+24,:,:].values
            # flatten obj mask to indicate only 1 or 0 (present / not present)
            boole   = (objs_L > 0)
            obj01 = boole.astype(int)
            # perform time averages to obtain daily means
            objs_d_L[d,:,:] = np.average( obj01, axis=0)    
            objprec_d_L[d,:,:] = np.average( obj01*prec, axis=0)    
            prec_d_L[d,:,:] = np.average( prec, axis=0)   
        
            print(f'{d}',end=',', flush=True)
    
        Xout = nc.wrt_objprec(year=year, exp='CAM7_LR', lat=latL, lon=lonL, objName=objn, objmsk=objs_d_L, objprec=objprec_d_L, prec=prec_d_L)
        print(f'finished CAM_LR {year}')
        
        objs_d_o= np.zeros( ( ndo, nyo, nxo ) )
        objprec_d_o= np.zeros( ( ndo, nyo, nxo ) )
        prec_d_o= np.zeros( ( ndo, nyo, nxo ) )
                
        for d in np.arange(ndo):
            # read in object masks and precip in one day chunks
            objs_o  = Xo[objn][d*24:d*24+24,:,:].values
            prec    = Xo['PR'][d*24:d*24+24,:,:].values
            # flatten obj mask to indicate only 1 or 0 (present / not present)
            boole   = (objs_o > 0)
            obj01   = boole.astype(int)
            # perform time averages to obtain daily means
            objs_d_o[d,:,:] = np.average( obj01 , axis=0)   # object 'presence' 
            objprec_d_o[d,:,:] = np.average( obj01 *prec, axis=0)  #precip from these objs 
            prec_d_o[d,:,:] = np.average( prec, axis=0)  # total precip
        
            print(f'{d}',end=',', flush=True)
        
        Xout = nc.wrt_objprec(year=year, exp='ERA5_Obs', lat=lato, lon=lono, objName=objn, objmsk=objs_d_o, objprec=objprec_d_o, prec=prec_d_o)
        print(f'finished ERA5 {year}')
        
        objs_d_c= np.zeros( ( ndc, nyc, nxc ) )
        objprec_d_c= np.zeros( ( ndc, nyc, nxc ) )
        prec_d_c= np.zeros( ( ndc, nyc, nxc ) )
        
        for d in np.arange(ndc):
            # read in object masks and precip in one day chunks
            objs_c  = Xc[objn][d*24:d*24+24,:,:].values
            prec    = Xc['PR'][d*24:d*24+24,:,:].values
            # flatten obj mask to indicate only 1 or 0 (present / not present)
            boole   = (objs_c > 0)
            obj01 = boole.astype(int)
            # perform time averages to obtain daily means
            objs_d_c[d,:,:] = np.average( obj01, axis=0)    
            objprec_d_c[d,:,:] = np.average( obj01*prec, axis=0)    
            prec_d_c[d,:,:] = np.average( prec, axis=0)   
        
            print(f'{d}',end=',', flush=True)
    
        Xout = nc.wrt_objprec(year=year, exp='CAM7_HR', lat=latc, lon=lonc, objName=objn, objmsk=objs_d_c, objprec=objprec_d_c, prec=prec_d_c)
        print(f'finished CAM_HR {year}')
        print(f"next ...", flush=True)


if __name__ == "__main__":
    drive()



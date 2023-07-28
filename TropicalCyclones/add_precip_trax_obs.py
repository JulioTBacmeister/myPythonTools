#!/usr/bin/env python
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature


import trax_util as trx
import TRMM_util as trmm
import ibtracs as IBT

import importlib

import time as TimeUtils

import sys
# import modules in other directories
sys.path.append('/glade/work/juliob/PyRegridding/Regridder/')
import argparse as arg

import esmfRegrid as erg

importlib.reload( trx )

def Parse_ymdh( ymdh ):
    year=ymdh//1000000
    mdh=ymdh-year*1000000
    month=mdh//10000
    dh=mdh-month*10000
    day=dh//100
    hour=dh-day*100
    return year,month,day,hour

def main(year0,year1):

    Ret=6378.1

    ##################
    # Get track file
    ###################
    power_wind=1.
    trk=IBT.readtrx()
    ##################
    # Get BaseName
    ###################

    BaseName = 'TRMM_IBTrACS'
    TrakFile = 'IBTrACS'
    
    print( f"BaseName is  {BaseName} " )
    print( f"TrakFile is  {TrakFile} " )


    nstorms,ntraxtime = np.shape( trk.lat )
    print( trk.hour[0,0:2])
    prectrk=np.zeros((nstorms,ntraxtime) )


    lonR=trk.lon.reshape( nstorms*ntraxtime )
    latR=trk.lat.reshape( nstorms*ntraxtime )
    yearR=trk.year.reshape( nstorms*ntraxtime )
    monthR=trk.month.reshape( nstorms*ntraxtime )
    dayR=trk.day.reshape( nstorms*ntraxtime )
    hourR=trk.hour.reshape( nstorms*ntraxtime )

    prectrkR = np.zeros( (nstorms*ntraxtime ))

    tic_overall = TimeUtils.perf_counter()

    basename = BaseName + 'cam.h4.PRECT.'
    # Directories
    drc='/glade/campaign/cgd/ccr/nanr/HPSS/tokeep/TIMESLICE/atm/proc/tseries/hourly3/PRECT/'

    for y in np.arange(start=year0,stop=year1+1):
        yA=str(y).zfill(4)
        fname   =  "/glade/p/cgd/amp/juliob/TRMM/3B42/timeseries/3B42_3hrly."+yA+".nc"
        dQ=xr.open_dataset( fname )
        nt,ny,nx = np.shape( dQ.precip )
        print(f"Doing {fname}")

        lonQd=dQ.lon.values
        latQd=dQ.lat.values
    

        for t in np.arange(nt):

            ymdh=dQ.yyyymmddhh[t].values.item()
            Tyear,Tmonth,Tday,Thour = Parse_ymdh( ymdh )

            prcQd_yx = dQ.precip[t,:,:].fillna(0).values
            lonQd0,prcQd_yx = trmm.shiftlon( lonQd, prcQd_yx )

            oop=np.where((yearR==Tyear) &
                     (monthR==Tmonth) &
                     (dayR==Tday) &
                     (hourR==Thour) )

            for u in oop[0]:
                if((yearR[u] == Tyear) and
                    (monthR[u] == Tmonth) and 
                    (dayR[u] == Tday) and
                    (hourR[u] == Thour) ):
                    prectrkR[u] = trx.Prec500( lonR[u],latR[u],lonQd0,latQd,prcQd_yx)

            if ((t % 100)==0):
                print( nt, t, Tyear,Tmonth,Tday,Thour )


    toc_overall = TimeUtils.perf_counter()
    pTime = f"Overall time creating prec-500   {toc_overall - tic_overall:0.4f} seconds"
    print(pTime)

    prectrk = np.reshape( prectrkR, (nstorms,ntraxtime))
    print(np.shape(prectrk))

    yAyA = str(year0).zfill(4)+'-'+str(year1).zfill(4)
    
    ## Should be done 
    ## np.savez('arrays.npz', arr1=arr1, arr2=arr2, arr3=arr3)
    npSaveFile = '/glade/work/juliob/NumPySaves/'+BaseName+'-prectrax-'+yAyA+'.npz'
    np.savez( npSaveFile , prectrk=prectrk  )


if __name__ == "__main__":
    my_parser = arg.ArgumentParser()
    my_parser.add_argument("--year0",     type=int, default=1979)
    my_parser.add_argument("--year1",     type=int, default=1979)
    my_parser.add_argument("--ens",       type=int, default=1)
    my_parser.add_argument("--sstA",      type=str, default='sst1' )
    args = my_parser.parse_args()
    main( year0=args.year0, year1=args.year1 )

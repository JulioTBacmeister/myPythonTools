#!/usr/bin/env python
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature


import trax_util as trx

import importlib

import time as TimeUtils

import sys
# import modules in other directories
sys.path.append('/glade/work/juliob/PyRegridding/Regridder/')
import argparse as arg


import esmfRegrid as erg


importlib.reload( trx )

def main(sstA,ens,year0,year1):

    Ret=6378.1

    scrip_dir = '/glade/p/cesmdata/cseg/inputdata/share/scripgrids/'
    dst_scrip = scrip_dir + 'fv0.23x0.31_141008.nc'
    src_scrip = scrip_dir + 'ne120np4_pentagons_100310.nc'

    dst1d_scrip = scrip_dir + 'fv0.9x1.25_141008.nc'

    #ne120 ==> latlon 0.25 degree
    regrd, srcf, dstf = erg.Regrid(srcScrip = src_scrip , 
                                    srcType  = 'mesh'  ,
                                    dstScrip = dst_scrip ,
                                    dstType  = 'grid'   )


    Scr=xr.open_dataset( dst_scrip )
    latQd_yx = Scr.grid_center_lat.values.reshape( (768, 1152 ) )
    lonQd_yx = Scr.grid_center_lon.values.reshape( (768, 1152 ) )
    lonQd = lonQd_yx[0,:]
    latQd = latQd_yx[:,0]

    ##################
    # Get track file
    ###################
    power_wind=(1./6.)**(0.11)

    if (year0>2023):
        TrakFile = trx.rcp85fname(sst=sstA) 
    else:
        TrakFile = trx.pdfname(ens=ens) 

    trk=trx.readtrx( TrakFile  , power_wind=power_wind )

    if (sstA == 'sst4'):
        trk=trx.Heal_rc4( trk )
        print(f" 'Fixed' years in SST4 tracks. Might wanna check " )
    if (sstA == 'sst6'):
        trk=trx.Heal_rc6( trk )
        print(f" 'Fixed' years in SST6 tracks. Might wanna check " )
    
    ##################
    # Get BaseName
    ###################
    if (year0>2023):
        BaseName = trx.rcp85fname(sst=sstA,justBaseName=True) 
    else:
        BaseName = trx.pdfname(ens=ens,justBaseName=True) 

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
        
    sys.stdout.flush()

    basename = BaseName + 'cam.h4.PRECT.'
    # Directories
    drc='/glade/campaign/cgd/ccr/nanr/HPSS/tokeep/TIMESLICE/atm/proc/tseries/hourly3/PRECT/'

    for y in np.arange(start=year0,stop=year1+1):
        yA=str(y).zfill(4)
        fname=drc+basename+yA+'010100Z-'+yA+'123121Z.nc'
        dQ=xr.open_dataset( fname )
        nt = dQ.dims['time']

        for t in np.arange(nt):
            #time=dQ.time_bnds[t,0].values.item()
            time=dQ.time[t].values.item()
            #print(time.year,time.month,time.day,time.hour)

            prcQd = dQ.PRECT[t,:].values

            prcQd_yx = erg.HorzRG(aSrc = prcQd , 
                        regrd = regrd , 
                        srcField=srcf , 
                        dstField=dstf , 
                        srcGridkey='c' ,
                        dstGridkey='yx' ) 

            oop=np.where((yearR==time.year) &
                         (monthR==time.month) &
                         (dayR==time.day) &
                         (hourR==time.hour) )

            for u in oop[0]:
                if((yearR[u] == time.year) and
                    (monthR[u] == time.month) and 
                    (dayR[u] == time.day) and
                    (hourR[u] == time.hour) ):
                    prectrkR[u] = trx.Prec500( lonR[u],latR[u],lonQd,latQd,prcQd_yx)

            if ((t % 100)==0):
                print( nt, t, time.year,time.month,time.day,time.hour )

        sys.stdout.flush()


    toc_overall = TimeUtils.perf_counter()
    pTime = f"Overall time creating prec-500   {toc_overall - tic_overall:0.4f} seconds"
    print(pTime)

    prectrk = np.reshape( prectrkR, (nstorms,ntraxtime))
    print(np.shape(prectrk))

    yAyA = str(year0).zfill(4)+'-'+str(year1).zfill(4)
    
    ## Should be done 
    ## np.savez('arrays.npz', arr1=arr1, arr2=arr2, arr3=arr3)
    npSaveFile = '/glade/work/juliob/NumPySaves/PrecTrax/'+BaseName+'-prectrax-'+yAyA+'.npz'
    np.savez( npSaveFile , prectrk=prectrk  )
    
    
if __name__ == "__main__":

    my_parser = arg.ArgumentParser()
    my_parser.add_argument("--year0",     type=int, default=1979)
    my_parser.add_argument("--year1",     type=int, default=1979)
    my_parser.add_argument("--ens",       type=int, default=1)
    my_parser.add_argument("--sstA",      type=str, default='sst1' )
    args = my_parser.parse_args()
    main( sstA=args.sstA,ens=args.ens,year0=args.year0,year1=args.year1 )

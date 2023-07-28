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

def write_netcdf(time,time_bnds,lat,lon,data,fname):
    #
    # time is xarray.DataArray.data
    # time_bnds xarray.DataArray
    #------------------------------------
    print( "No dims specified in write_netdcf. Inheriting frpm dataarrays " )
    coords = dict( 
        lon  = ( ["lon"],lon ),
        lat  = ( ["lat"],lat ),
        time = ( ["time"], time ), 
            )
    dS = xr.Dataset( coords=coords  )
    dS["time_bnds"]=time_bnds
    Dar = xr.DataArray( data=data, dims=('time','lat','lon',),
                        attrs=dict( description='TC-precipitation',units='ms-1',) ,) 
    dS["TCPRECT"]=Dar
    
    dS.to_netcdf( fname )

    
def main(sstA,ens,year0,year1):


    print(f" Doing it again ... again ... ")
    
    scrip_dir = '/glade/p/cesmdata/cseg/inputdata/share/scripgrids/'
    dst_scrip = scrip_dir + 'fv0.23x0.31_141008.nc'
    src_scrip = scrip_dir + 'ne120np4_pentagons_100310.nc'
    dst1d_scrip = scrip_dir + 'fv0.9x1.25_141008.nc' #Not used at the moment

    print( f"About to create ESMF regridding objects ")
    #ne120 ==> latlon 0.25 degree
    regrd, srcf, dstf = erg.Regrid(srcScrip = src_scrip , 
                                    srcType  = 'mesh'  ,
                                    dstScrip = dst_scrip ,
                                    dstType  = 'grid'   )

    DstFieldData = dstf.data
    nx,ny = np.shape( DstFieldData ) # Note ESMF objects are transposed 
    print(f"Size of precip grid {nx,ny}")

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
    windR=trk.hour.reshape( nstorms*ntraxtime )

    tic_overall = TimeUtils.perf_counter()

    basename = BaseName + 'cam.h4.PRECT.'
    basename_o = BaseName + 'cam.h4.TCPRECT_xTS.'
    # input and output Directories
    drc='/glade/campaign/cgd/ccr/nanr/HPSS/tokeep/TIMESLICE/atm/proc/tseries/hourly3/PRECT/'
    drc_o='/glade/campaign/cgd/amp/juliob/TC-cesm1/precip/'

    for y in np.arange(start=year0,stop=year1+1):
        yA=str(y).zfill(4)
        fname=drc+basename+yA+'010100Z-'+yA+'123121Z.nc'
        fname_o=drc_o+basename_o+yA+'010100Z-'+yA+'123121Z.nc'
        dQ=xr.open_dataset( fname )
        print( f"Opened file {fname}" )
        nt = dQ.dims['time']
        prec_TC = np.zeros( (nt , ny, nx ) )

        for t in np.arange(nt):
            #for t in np.arange(110):
            time=dQ.time[t].values.item()

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
                    (hourR[u] == time.hour) and 
                    (windR[u] >= 18.0) ):
                    precTC,js,jn,iw,ie = trx.Prec500grid( lonR[u],latR[u],lonQd,latQd,prcQd_yx)
                    #print( np.shape(precTC),np.shape(prec_TC[t,js:jn,iw:ie]),js,jn,iw,ie )
                    prec_TC[t,js:jn,iw:ie] = prec_TC[t,js:jn,iw:ie] + precTC


            if ((t % 100)==0):
                print( nt, t, time.year,time.month,time.day,time.hour )

        write_netcdf(dQ.time.data,dQ.time_bnds,latQd,lonQd,prec_TC,fname_o)
        print( f"Wrote file {fname_o}" )

    toc_overall = TimeUtils.perf_counter()
    pTime = f"Overall time creating prec-500   {toc_overall - tic_overall:0.4f} seconds"
    
if __name__ == "__main__":

    my_parser = arg.ArgumentParser()
    my_parser.add_argument("--year0",     type=int, default=1979)
    my_parser.add_argument("--year1",     type=int, default=1979)
    my_parser.add_argument("--ens",       type=int, default=1)
    my_parser.add_argument("--sstA",      type=str, default='sst1' )
    args = my_parser.parse_args()
    main( sstA=args.sstA,ens=args.ens,year0=args.year0,year1=args.year1 )

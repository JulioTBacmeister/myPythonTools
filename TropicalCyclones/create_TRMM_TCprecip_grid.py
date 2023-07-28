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


import esmfRegrid as erg

from datetime import datetime




importlib.reload( trmm )

def write_netcdf(yyyymmddhh,lat,lon,data,fname):
    #dims   = ["lon","lat","time"]
    coords = dict( 
        lon  = ( ["lon"],lon ),
        lat  = ( ["lat"],lat ),
        time = ( ["time"], yyyymmddhh ), 
            )
    dS = xr.Dataset( coords=coords  )
    #dS["time_bnds"]=time_bnds
    Dar = xr.DataArray( data=data, dims=('time','lat','lon',),
                        attrs=dict( description='TC-precipitation in TRMM',units='mm hr-1',)) 
    dS["TCPRECT"]=Dar
    
    dS.attrs["Description"]= ("TRMM 3B42 merged with IBTrACS.since1980.v04r00.nc. This file contains 3hrly precipitation in a ~500km radius from IBTrACS storm centers."
                              +" Longitudes have been shifted to 0-360.  TRMM NaNs are replaced with Zero." )

    
    # Get the current date
    current_date = datetime.now()
    # Format it as "yyyy-mm-dd"
    formatted_date = current_date.strftime('%Y-%m-%d')
    dS.attrs["CreationDate"]=formatted_date
    
    dS.to_netcdf( fname )
    
    
def Parse_ymdh( ymdh ):
    year=ymdh//1000000
    mdh=ymdh-year*1000000
    month=mdh//10000
    dh=mdh-month*10000
    day=dh//100
    hour=dh-day*100
    return year,month,day,hour

    




year0=1998
year1=2019


##################
# Get track file
###################
power_wind=1. 


trk=IBT.readtrx()


nstorms,ntraxtime = np.shape( trk.lat )
print( trk.hour[0,0:2])
prectrk=np.zeros((nstorms,ntraxtime) )


lonR=trk.lon.reshape( nstorms*ntraxtime )
latR=trk.lat.reshape( nstorms*ntraxtime )
yearR=trk.year.reshape( nstorms*ntraxtime )
monthR=trk.month.reshape( nstorms*ntraxtime )
dayR=trk.day.reshape( nstorms*ntraxtime )
hourR=trk.hour.reshape( nstorms*ntraxtime )
windR=trk.wind.reshape( nstorms*ntraxtime )


tic_overall = TimeUtils.perf_counter()

for y in np.arange(start=year0,stop=year1+1):
    yA=str(y).zfill(4)
    fname   =  "/glade/p/cgd/amp/juliob/TRMM/3B42/timeseries/3B42_3hrly."+yA+".nc"
    fname_o =  "/glade/p/cgd/amp/juliob/TRMM/3B42_TCs/timeseries/3B42_3hrly_TCPRECT_xTS."+yA+".nc"
    dQ=xr.open_dataset( fname )
    print( f"Opened {fname}" )
    nt,ny,nx = np.shape( dQ.precip )
    prec_TC = np.zeros( (nt , ny, nx ) )
    
    lonQd=dQ.lon.values
    latQd=dQ.lat.values
    
    #lonQd0=np.roll(lonQd,shift=720,axis=0)
    #oo=np.where(lonQd0<0)
    #lonQd0[oo[0]]=lonQd0[oo[0]]+360.
    
    for t in np.arange(nt):
        #for t in np.arange(110):
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
                (hourR[u] == Thour) and 
                (windR[u] >= 18.0) ):
                precTC,js,jn,iw,ie = trx.Prec500grid( lonR[u],latR[u],lonQd0,latQd,prcQd_yx)
                #print( np.shape(precTC),np.shape(prec_TC[t,js:jn,iw:ie]),js,jn,iw,ie )
                prec_TC[t,js:jn,iw:ie] = prec_TC[t,js:jn,iw:ie] + precTC
            
                
        if ((t % 100)==0):
            print( nt, t, Tyear,Tmonth,Tday,Thour )
            
    write_netcdf(dQ.yyyymmddhh.data,latQd,lonQd0,prec_TC,fname_o)
    print( f"Wrote {fname_o}" )

toc_overall = TimeUtils.perf_counter()
pTime = f"Overall time creating prec-500   {toc_overall - tic_overall:0.4f} seconds"

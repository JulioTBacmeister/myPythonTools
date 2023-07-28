import xarray as xr
import numpy as np

import importlib
import cftime

def TS(year0=2006,year1=2100):
    
    
    drc='/glade/campaign/cesm/collections/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/TS/'
    fils='b.e11.BRCP85C5CNBDRD.f09_g16.0**.cam.h0.TS.208101-210012.nc'
    dse=xr.open_mfdataset( drc+fils ,concat_dim='ensemble',combine='nested')
    lats=dse.lat.values
    lons=dse.lon.values


    drcLF='/glade/campaign/cesm/collections/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/LANDFRAC/'
    filLF='b.e11.BRCP85C5CNBDRD.f09_g16.001.cam.h0.LANDFRAC.208101-210012.nc'
    dsLF=xr.open_dataset( drcLF+filLF )
    landf=dsLF.LANDFRAC.values[0,:,:]

    
def FullTS():    
    
    ts,landf,lat,lon = bespokeTS()
    drc='/glade/collections/cdg/data/cmip5/output1/NSF-DOE-NCAR/CESM1-CAM5/rcp85/mon/atmos/Amon/r3i1p1/files/ts_20140129/'
    file=drc+'ts_Amon_CESM1-CAM5_rcp85_r3i1p1_200601-210012.nc'
    dS1 = xr.open_dataset( file )
    ts1 = dS1.ts.values[-31*12:,:,:]
    nt,ny,nx = np.shape( ts1 )
    ts1=ts1.reshape(1,nt,ny,nx)
    tsens = (ts.mean( axis=0 )).reshape(1,nt,ny,nx)
    #print( np.shape(ts1), np.shape(ts), np.shape(tsens) )
    ts_x = np.concatenate( (ts1, ts, tsens), axis=0 )
    print( "SHape of returned TS: ",np.shape(ts_x) )

    return ts_x,landf,lat,lon

def PresentDayTS():    
    
    drc='/glade/campaign/cesm/collections/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/TS/'

    fils='b.e11.B20TRC5CNBDRD.f09_g16.0**.cam.h0.TS.192001-200512.nc'
    dse0=xr.open_mfdataset( drc+fils ,concat_dim='ensemble',combine='nested')
    lat=dse0.lat.values
    lon=dse0.lon.values

    print( "This is an exmaple of how you can begin to"+ 
        " use the god-awful impenetrable shit-show that "+
        " are cftime" )
    poop = dse0.time_bnds[0,0,1].values.item()
    print(poop.year,poop.month,poop.day,poop.hour)

    drcLF='/glade/campaign/cesm/collections/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/LANDFRAC/'
    filLF='b.e11.BRCP85C5CNBDRD.f09_g16.001.cam.h0.LANDFRAC.208101-210012.nc'
    dsLF=xr.open_dataset( drcLF+filLF )
    landf =dsLF.LANDFRAC.values[0,:,:]

    ts0 = dse0.TS[:,-30*12:,:,:].values
 
    #ts_x = np.concatenate( (ts0, ts1), axis=1 )
    ts_x = ts0

    return ts_x,landf,lat,lon


def bespokeTS():
    
    # This is a quick and dirty function to return the LENS1 TS fields for 2070-2100
    
    # Here is where LENS1 lives
    drc='/glade/campaign/cesm/collections/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/TS/'

    
    # Get the last 20 years of 21st C
    fils='b.e11.BRCP85C5CNBDRD.f09_g16.0**.cam.h0.TS.208101-210012.nc'
    dse1=xr.open_mfdataset( drc+fils ,concat_dim='ensemble',combine='nested')
    lat=dse1.lat.values
    lon=dse1.lon.values

    # Get 2006-2080
    fils='b.e11.BRCP85C5CNBDRD.f09_g16.0**.cam.h0.TS.200601-208012.nc'
    dse0=xr.open_mfdataset( drc+fils ,concat_dim='ensemble',combine='nested')

    print( "This is an exmaple of how you can begin to"+ 
        " use the god-awful impenetrable shit-show that "+
        " are cftime" )
    poop = dse0.time_bnds[0,0,1].values.item()
    print(poop.year,poop.month,poop.day,poop.hour)

    # Get LANDFRAC for LENS
    drcLF='/glade/campaign/cesm/collections/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/LANDFRAC/'
    filLF='b.e11.BRCP85C5CNBDRD.f09_g16.001.cam.h0.LANDFRAC.208101-210012.nc'
    dsLF=xr.open_dataset( drcLF+filLF )
    landf =dsLF.LANDFRAC.values[0,:,:]

    # Pick out last 11 years of 2006-2080
    # and pre-pend to 2081-2100
    ts0 = dse0.TS[:,900-132:,:,:].values
    ts1 = dse1.TS.values

    ts_x = np.concatenate( (ts0, ts1), axis=1 )

    return ts_x,landf,lat,lon
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
    
    # This is a quick and dirty function to return the LENS2 TS fields for 2070-2100
    
    # Here is where LENS2 TS lives
    drc=lens2dir()+'TS/'
    
    """
    #EXAMPLE
    b.e21.BSSP370cmip6.f09_g17.LE2-1281.005.cam.h0.TS.205501-206412.nc
    b.e21.BSSP370cmip6.f09_g17.LE2-1281.005.cam.h0.TS.206501-207412.nc
    b.e21.BSSP370cmip6.f09_g17.LE2-1281.005.cam.h0.TS.207501-208412.nc
    b.e21.BSSP370cmip6.f09_g17.LE2-1281.005.cam.h0.TS.208501-209412.nc
    b.e21.BSSP370cmip6.f09_g17.LE2-1281.005.cam.h0.TS.209501-210012.nc
    """
    
    date_tags=[ '206501-207412', '207501-208412', '208501-209412', '209501-210012' ]
    
    iRd=0
    for tag in date_tags:
        fils='b.e21.BSSP370cmip6.f09_g17.LE2-*.*.cam.h0.TS.'+tag+'.nc'
        print( drc+fils )
        dseN=xr.open_mfdataset( drc+fils ,concat_dim='ensemble',combine='nested')
        tsN = dseN.TS.values
        if (iRd==0):
            ts_x=tsN
        else:
            ts_x =  np.concatenate( (ts_x, tsN), axis=1 )
        iRd = iRd+1
    
    lat=dseN.lat.values
    lon=dseN.lon.values

    drcLF=lens2dir()+'LANDFRAC/'
    filLF='b.e21.BSSP370cmip6.f09_g17.LE2-1281.005.cam.h0.LANDFRAC.205501-206412.nc'
    print( f"Getting {drcLF+filLF}" )
    dsLF=xr.open_dataset( drcLF+filLF )
    landf =dsLF.LANDFRAC.values[0,:,:]

    ts_x = ts_x[:,5*12:,:,:]

    return ts_x,landf,lat,lon

def longtermTS():
    
    # This is a quick and dirty function to return the LENS2 TS fields for 2070-2100
    
    # Here is where LENS2 TS lives
    drc=lens2dir()+'TS/'
    
    """
    #EXAMPLE
    b.e21.BSSP370cmip6.f09_g17.LE2-1281.005.cam.h0.TS.205501-206412.nc
    b.e21.BSSP370cmip6.f09_g17.LE2-1281.005.cam.h0.TS.206501-207412.nc
    b.e21.BSSP370cmip6.f09_g17.LE2-1281.005.cam.h0.TS.207501-208412.nc
    b.e21.BSSP370cmip6.f09_g17.LE2-1281.005.cam.h0.TS.208501-209412.nc
    b.e21.BSSP370cmip6.f09_g17.LE2-1281.005.cam.h0.TS.209501-210012.nc
    """

    """
    date_tags=[ '200501-201412', 
                '201501-203412', 
                '202501-203412', 
                '203501-204412', 
                '204501-205412', 
                '205501-206412', 
                '206501-207412', 
               '207501-208412', 
               '208501-209412', 
               '209501-210012' ]
    """
 
    date_tags=[ '208501-209412', 
               '209501-210012' ]
    
    iRd=0
    for tag in date_tags:
        fils='b.e21.BSSP370cmip6.f09_g17.LE2-*.*.cam.h0.TS.'+tag+'.nc'
        print( drc+fils )
        dseN=xr.open_mfdataset( drc+fils ,concat_dim='ensemble',combine='nested')
        tsN = dseN.TS.values
        if (iRd==0):
            ts_x=tsN
        else:
            ts_x =  np.concatenate( (ts_x, tsN), axis=1 )
        iRd = iRd+1
    
    lat=dseN.lat.values
    lon=dseN.lon.values

    drcLF=lens2dir()+'LANDFRAC/'
    filLF='b.e21.BSSP370cmip6.f09_g17.LE2-1281.005.cam.h0.LANDFRAC.205501-206412.nc'
    print( f"Getting {drcLF+filLF}" )
    dsLF=xr.open_dataset( drcLF+filLF )
    landf =dsLF.LANDFRAC.values[0,:,:]

    ts_x = ts_x[:,5*12:,:,:]

    return ts_x,landf,lat,lon

def decadalTS( contains ):
    
    # This is a quick and dirty function to return the LENS2 TS fields for 2070-2100
    
    # Here is where LENS2 TS lives
    drc=lens2dir()+'TS/'
    
    """
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BHISTcmip6.f09_g17.LE2-1301.001.cam.h0.TS.186001-186912.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BHISTcmip6.f09_g17.LE2-1301.001.cam.h0.TS.187001-187912.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BHISTcmip6.f09_g17.LE2-1301.001.cam.h0.TS.188001-188912.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BHISTcmip6.f09_g17.LE2-1301.001.cam.h0.TS.189001-189912.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BHISTcmip6.f09_g17.LE2-1301.001.cam.h0.TS.190001-190912.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BHISTcmip6.f09_g17.LE2-1301.001.cam.h0.TS.191001-191912.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BHISTcmip6.f09_g17.LE2-1301.001.cam.h0.TS.192001-192912.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BHISTcmip6.f09_g17.LE2-1301.001.cam.h0.TS.193001-193912.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BHISTcmip6.f09_g17.LE2-1301.001.cam.h0.TS.194001-194912.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BHISTcmip6.f09_g17.LE2-1301.001.cam.h0.TS.195001-195912.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BHISTcmip6.f09_g17.LE2-1301.001.cam.h0.TS.196001-196912.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BHISTcmip6.f09_g17.LE2-1301.001.cam.h0.TS.197001-197912.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BHISTcmip6.f09_g17.LE2-1301.001.cam.h0.TS.198001-198912.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BHISTcmip6.f09_g17.LE2-1301.001.cam.h0.TS.199001-199912.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BHISTcmip6.f09_g17.LE2-1301.001.cam.h0.TS.200001-200912.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BHISTcmip6.f09_g17.LE2-1301.001.cam.h0.TS.201001-201412.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BSSP370cmip6.f09_g17.LE2-1301.001.cam.h0.TS.201501-202412.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BSSP370cmip6.f09_g17.LE2-1301.001.cam.h0.TS.202501-203412.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BSSP370cmip6.f09_g17.LE2-1301.001.cam.h0.TS.203501-204412.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BSSP370cmip6.f09_g17.LE2-1301.001.cam.h0.TS.204501-205412.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BSSP370cmip6.f09_g17.LE2-1301.001.cam.h0.TS.205501-206412.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BSSP370cmip6.f09_g17.LE2-1301.001.cam.h0.TS.206501-207412.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BSSP370cmip6.f09_g17.LE2-1301.001.cam.h0.TS.207501-208412.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BSSP370cmip6.f09_g17.LE2-1301.001.cam.h0.TS.208501-209412.nc
    /glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/TS/b.e21.BSSP370cmip6.f09_g17.LE2-1301.001.cam.h0.TS.209501-210012.nc
    """

    """
    date_tags=[ '200501-201412', 
                '201501-203412', 
                '202501-203412', 
                '203501-204412', 
                '204501-205412', 
                '205501-206412', 
                '206501-207412', 
               '207501-208412', 
               '208501-209412', 
               '209501-210012' ]
    """
 
    if ( contains >= 1860 and contains <= 2009 ):
        date_tags=[ f'{contains}01-{contains+9}12',] 
    elif ( contains >= 2010 and contains <= 2014 ):
        date_tags=[ f'{contains}01-{contains+4}12',] 
    elif ( contains >= 2015 and contains <= 2094 ):
        date_tags=[ f'{contains}01-{contains+9}12',] 
    elif ( contains >= 2095 and contains <= 2100 ):
        date_tags=[ f'{contains}01-{contains+5}12',]
    else:
        print(f'No files contain {contains}')
        return -999,-999

    if (contains <=2014):
        rootName = 'b.e21.BHISTcmip6.f09_g17.LE2-'
    else:
        rootName = 'b.e21.BSSP370cmip6.f09_g17.LE2-'

    iRd=0
    for tag in date_tags:
        #fils='b.e21.BSSP370cmip6.f09_g17.LE2-*.*.cam.h0.TS.'+tag+'.nc'
        fils=f'{rootName}*.*.cam.h0.TS.{tag}.nc'
        print( drc+fils )
        
        dseN=xr.open_mfdataset( drc+fils ,concat_dim='ensemble',combine='nested')
        tsN = dseN.TS.values
        if (iRd==0):
            ts_x=tsN
        else:
            ts_x =  np.concatenate( (ts_x, tsN), axis=1 )
        iRd = iRd+1

    lat=dseN.lat.values
    lon=dseN.lon.values

    drcLF=lens2dir()+'LANDFRAC/'
    filLF='b.e21.BSSP370cmip6.f09_g17.LE2-1281.005.cam.h0.LANDFRAC.205501-206412.nc'
    print( f"Getting {drcLF+filLF}" )
    dsLF=xr.open_dataset( drcLF+filLF )
    landf =dsLF.LANDFRAC.values[0,:,:]

    return ts_x,landf,lat,lon



######################################
def lens2dir():
    # Here is where LENS2 lives
    drc='/glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/month_1/'
    
    return drc


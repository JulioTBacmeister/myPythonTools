#!/usr/bin/env python

import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import scipy.interpolate as si
import argparse as arg
from pathlib import Path
import subprocess as sp

"""
            This module is supposed to interpolate results from regridded
            SAMwrf ne30x16 ==> ne30 data in the vertical to a default ne30 vertical grid. 
            The issue is that the regridded data is on "rougher" topography than 
            the default ne30 model.  Assume you know of nthing more about the 
            default ne30 grid than phis.
"""


def regrid(year,month,day,hour,reftopopath,ifname,idir):

    y4=str(year).zfill(4)
    m2=str(month).zfill(2)
    d2=str(day).zfill(2)
    s5=str( hour*3600 ).zfill(5)
    xtnc = y4+'-'+m2+'-'+d2+'-'+s5+'.nc'




    ofname = ifname+'_vrg'
    odir   = idir+'vrg/'
    sp.run( 'mkdir -p '+odir, shell=True )

    fM= reftopopath  

    fR=idir + ifname + '.cam.h1.'+xtnc    

    fMo=odir + ofname + '.cam.h1.'+xtnc 

    print(fR)

    #plt.ion()

    aM=xr.open_mfdataset( fM )  # def ne30
    aR=xr.open_mfdataset( fR )  # rough from ne30x16

    print(list(aR.variables))

    dimsM=aM.dims
    dimsR=aR.dims
    assert(dimsM['ncol'] == dimsR['ncol']), "ncol is different in the datasets. Can't proceed." 

    ncols=dimsR['ncol']

    lon=aR['lon']
    lat=aR['lat']

    hyai=aR['hyai'].values
    hyam=aR['hyam'].values
    hybi=aR['hybi'].values
    hybm=aR['hybm'].values

    #get surface geopotential
    phMx=aM['PHIS']
    phR=aR['PHIS']
    phM=phMx[0,:]

    phM=phM.values
    phR=phR.values

    dph=phR - phM

    #get Temperature on "rough" grid.
    TeR=aR['T'].values

    #get Specific humidity on "rough" grid.
    qR=aR['Q'].values
    
    #get U,V winds on "rough" grid.
    uR=aR['U'].values
    vR=aR['V'].values

    #get surface pressure
    psMx=aM['PS']
    psR=aR['PS']
    psMq=psMx[0,:]

    psMq=psMq.values
    psR=psR.values
    dps=psR-psMq

    fig, (ax1, ax2) = plt.subplots(nrows=2)
    vlon=lon.values
    vlat=lat.values

    #vdps=dps.values
    #ax1.tricontour(vlon,vlat, vdps, levels=14, linewidths=0.5, colors='k')
    #ax2.tricontourf(lon, lat, vdps , levels=30, cmap="RdBu_r")
    #plt.show()

    Tes=TeR[31,:]
    psMo=psR * np.exp ( (phR-phM) / (287.*Tes ) )
            
    dims=psR.shape

     #imin=dps.argmin()  # Picks out a point where rough topo is substantially higher than default topo

    TeMo = TeR
    uMo = uR
    vMo = vR
    qMo = qR

    for icol in np.arange(ncols):
        ptrg=hyam*100000. + hybm*psMo[icol]
        psrc=hyam*100000. + hybm*psR[icol]
        if ( ( phR[icol]>10.*9.8 or phM[icol]>10.*9.8)  and vlon[icol]>250. and vlon[icol]<350. and vlat[icol]>-60 and vlat[icol]<17.5 ):
            #Temperature
            f=si.interp1d(psrc,TeR[:,icol], fill_value='extrapolate')
            TeMo[:,icol]=f(ptrg)
            #Humidity
            f=si.interp1d(psrc,qR[:,icol], fill_value='extrapolate'  )
            qMo[:,icol]=f(ptrg)

        if ( phR[icol]>phM[icol] and vlon[icol]>250. and vlon[icol]<350. and vlat[icol]>-60 and vlat[icol]<17.5):
            #Zonal wind
            f=si.interp1d(psrc,uR[:,icol], fill_value= 0.,bounds_error=False  )
            uMo[:,icol]=f(ptrg)
            #Meridional wind
            f=si.interp1d(psrc,vR[:,icol], fill_value= 0.,bounds_error=False  )
            vMo[:,icol]=f(ptrg)

        if( (icol % 10000) == 0 ):
            print("Regridded through icol=",icol)


    aMo=xr.Dataset( coords=aR.coords )

    # add PHIS
    xphMo= xr.DataArray( data=phM, dims=aR['PHIS'].dims, coords=aR['PHIS'].coords, attrs=dict( description='Surface Geopotential Height',units='m+2 s-2',) ,) 
    aMo['PHIS']=xphMo
    # add PS
    xpsMo= xr.DataArray( data=psMo, dims=aR['PS'].dims, coords=aR['PS'].coords, attrs=dict( description='Surface Pressure',units='Pa',) ,) 
    aMo['PS']=xpsMo
    # add T
    xTeMo= xr.DataArray( data=TeMo, dims=aR['T'].dims, coords=aR['T'].coords, attrs=dict( description='Temperature',units='K',) ,) 
    aMo['T']=xTeMo
    # add Q
    xqMo= xr.DataArray( data=qMo, dims=aR['Q'].dims, coords=aR['Q'].coords, attrs=dict( description='Specific Humidity',units='kg kg-1',) ,) 
    aMo['T']=xTeMo
    # add U
    xuMo= xr.DataArray( data=uMo, dims=aR['U'].dims, coords=aR['U'].coords, attrs=dict( description='Zonal wind',units='m s-1',) ,) 
    aMo['U']=xuMo
    xvMo= xr.DataArray( data=vMo, dims=aR['V'].dims, coords=aR['V'].coords, attrs=dict( description='Meridional wind',units='m s-1',) ,) 
    aMo['V']=xvMo

    aMo.to_netcdf(fMo)

def main(year=2010,month=6):
    
    days_in_month =[31 , 28, 31, 30, 31, 30, 31, 31, 30, 31,30, 31 ]
    imm=month

    reftopopath ='/project/amp/juliob/CAM/c6_3_59.ne30pg3_L32_SAMwrf.ndg01/c6_3_59.ne30pg3_L32_SAMwrf.ndg01.cam.h0.2010-06.nc'
    ifname = 'f.e22r.SAMwrf01.ne30.L32.NODEEP_2010_01'
    idir   = '/project/amp/juliob/CAM/SAMwrf.ne30_L32/'

    ndays = days_in_month[imm-1]
    for idd in np.arange(1,ndays+1):
        for ihh in np.arange(24):
            regrid(year=year,month=month,day=idd,hour=ihh,idir=idir,ifname=ifname,reftopopath=reftopopath )
            #print( year,month,idd,ihh)


#def main(year,month):

    #month(year,month)

if __name__ == "__main__":
    # argument: indir -> get all nc files in this directory
    # argument: map -> the offlinemap file already prepared
    # argument: outdir -> directory where remapped files should go
    my_parser = arg.ArgumentParser()
    my_parser.add_argument("--month", type=int)
    my_parser.add_argument("--year", type=int)
    args = my_parser.parse_args()
    main(args.year, args.month )

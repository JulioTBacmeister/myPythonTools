#!/usr/bin/env python
# Import packages 
import sys
import os 

import xarray as xr
import numpy as np
import pandas as pd
import scipy

import time as TimeUtil

from PyRegridding.Regridder.GlobalVarClass import Gv
from PyRegridding.Utils import MyConstants as Co



def write_netcdf( year=None , month=None  ):

    grav=Co.grav()
    pdTime_ERA = Gv.pdTime_ERA
    TimeStamp = Gv.SrcTimeStamp
    
    ps_ERA = Gv.ps_ERA 
    q_ERA = Gv.q_ERA  
    v_ERA = Gv.v_ERA  

    delp_ERA = Gv.delp_ERA 

    vq_ERA = q_ERA * v_ERA

    dpvq_ERA = delp_ERA * q_ERA * v_ERA

    lon_ERA = Gv.lon_ERA 
    lat_ERA = Gv.lat_ERA 

    #-----------------------------------------------
    # Get hybrid eta-coordinate coefficients for ERA5
    #-----------------------------------------------
    a_model = Gv.amid_ERA 
    b_model = Gv.bmid_ERA 
    a_half = Gv.aint_ERA  
    b_half = Gv.bint_ERA 
    print( "shape of a_model ", np.shape( Gv.amid_ERA ) ) 

    
    SuperDir = "/glade/campaign/cgd/amp/juliob/ERA5"

    ntime = np.shape(pdTime_ERA)[0]
    print(ntime)
    #Bdiro="/glade/derecho/scratch/juliob/ERA5/" + Gv.MyDst    
    #Bdiro=f"{SuperDir}/Fluxes/VQ/"
    Bdiro=f"{SuperDir}/Fluxes/VQ/hourly/{str(year).zfill(4)}-{str(month).zfill(2)}/"
    #######
    os.makedirs( Bdiro , exist_ok=True )

    #version='test_netcdf4_default'
    #Bfilo="/glade/derecho/scratch/juliob/ERA5/" + Gv.MyDst + "/" + Gv.MySrc +"_x_"+ Gv.MyDst + "_"+ Gv.MyDstVgrid + "_" + version 
    Bfilo= f"{Bdiro}/e5.oper.an.ml.VQ"

    ilev = a_half  + b_half  * 100_000. #* 100_000.
    lev  = a_model + b_model * 100_000. #* 100_000.

    
    dims   = ["lon","lat","time","lev","ilev"]
    coords = dict( 
        longitude  = ( ["longitude"],lon_ERA ),
        latitude   = ( ["latitude"],lat_ERA ),
        lev  = ( ["lev"],lev),
        ilev = ( ["ilev"],ilev),
        time = ( ["time"],pdTime_ERA ) , #pd.to_datetime( pdTime_ERA[itim] ) ),
    )

    Wds = xr.Dataset( coords=coords  )
    #Wds["Time"] = pdTime_ERA[itim] )
    #Wds["P_00"] = 100_000.

    Dar = xr.DataArray( data=a_half, dims=('ilev',),
                        attrs=dict( description='interface hybrid eta coordinate A-coeff ',units='1',) ,) 
    Wds['a_half'] = Dar

    Dar = xr.DataArray( data=b_half, dims=('ilev',),
                        attrs=dict( description='interface hybrid eta coordinate B-coeff ',units='1',) ,) 
    Wds['b_half'] = Dar

    Dar = xr.DataArray( data=a_model , dims=('lev',),
                        attrs=dict( description='mid-level hybrid eta coordinate A-coeff ',units='1',) ,) 
    Wds['a_model'] = Dar

    Dar = xr.DataArray( data=b_model , dims=('lev',),
                        attrs=dict( description='mid-level hybrid eta coordinate B-coeff ',units='1',) ,) 
    Wds['b_model'] = Dar


    Dar = xr.DataArray( data=ps_ERA , 
                        dims=('time','latitude','longitude',),
                        attrs=dict( description='Surface Pressure',units='Pa',) ,) 
    Wds['SP'] = Dar

    Dar = xr.DataArray( data=vq_ERA , 
                        dims=('time','lev','latitude','longitude',),
                        attrs=dict( description='Meridional water flux',units='m s-1 (kg kg-1)',) ,) 
    Wds['VQ'] = Dar

    Dar = xr.DataArray( data=delp_ERA , 
                        dims=('time','lev', 'latitude','longitude',),
                        attrs=dict( description='Pressure thickness',units='Pa',) ,) 
    Wds['DELP'] = Dar

    Dar = xr.DataArray( data=dpvq_ERA , 
                        dims=('time','lev', 'latitude','longitude',),
                        attrs=dict( description='Pressure thickness',units='Pa',) ,) 
    Wds['dPVQ'] = Dar

    """
    Dar = xr.DataArray( data=pint_ERA , 
                        dims=('time','ilev', 'latitude','longitude',),
                        attrs=dict( description='Ifc. Pressure',units='Pa',) ,) 
    Wds['PINT'] = Dar
    """

    filo= f"{Bfilo}.{TimeStamp}.nc"
    Wds.to_netcdf( filo ) #,format="NETCDF3_CLASSIC" )
    print( f' .. Wrote {filo}' , flush=True)

    code = 1
    return code


def write_daily_netcdf(year=None, month=None, day=None, return_dataset=False ):

    SuperDir = f"/glade/campaign/cgd/amp/juliob/ERA5" 
    #e5.oper.an.ml.VQ.1990101006_1990101011.nc
    Bdiri=f"{SuperDir}/Fluxes/VQ/hourly/{str(year).zfill(4)}-{str(month).zfill(2)}/"
    Bfili= f"{Bdiri}/e5.oper.an.ml.VQ.{str(year).zfill(4)}{str(month).zfill(2)}{str(day).zfill(2)}*"
    Bdiro=f"{SuperDir}/Fluxes/VQ/daily/{str(year).zfill(4)}-{str(month).zfill(2)}/"
    Bfilo= f"{Bdiro}/e5.oper.an.ml.VQ.{str(year).zfill(4)}-{str(month).zfill(2)}-{str(day).zfill(2)}.nc"
    #######
    os.makedirs( Bdiro , exist_ok=True )


    X = xr.open_mfdataset( Bfili, data_vars='different', coords='different' , parallel=False  ) # No DASK Bullshit
    
    
    Xa = X.mean(dim='time').compute()
    Xa['DateStamp'] = f'{str(year).zfill(4)}-{str(month).zfill(2)}-{str(day).zfill(2)}'
    
    da = xr.DataArray( X.time.values, dims=['times'] , name='times_in_avg' )
    Xa['times_in_avg' ]=da
    
    if (return_dataset==False):
        Xa.to_netcdf( Bfilo )
        print(f' DATA Written to {Bfilo}' , flush=True )
        rcode=1
        return rcode
    else:
        print(f' Would be written to {Bfilo}' , flush=True )
        return Xa,X

def write_monthly_netcdf(year=None, month=None, return_dataset=False ):

    SuperDir = f"/glade/campaign/cgd/amp/juliob/ERA5" 
    #e5.oper.an.ml.VQ.1990101006_1990101011.nc
    Bdiri=f"{SuperDir}/Fluxes/VQ/daily/{str(year).zfill(4)}-{str(month).zfill(2)}/"
    Bfili= f"{Bdiri}/e5.oper.an.ml.VQ.{str(year).zfill(4)}-{str(month).zfill(2)}-*.nc"
    Bdiro=f"{SuperDir}/Fluxes/VQ/monthly/{str(year).zfill(4)}/"
    Bfilo= f"{Bdiro}/e5.oper.an.ml.VQ.{str(year).zfill(4)}-{str(month).zfill(2)}.nc"
    #######
    os.makedirs( Bdiro , exist_ok=True )


    X=xr.open_mfdataset( Bfili, combine='nested',
                            concat_dim='time',
                            data_vars='different',
                            coords='different')
    
    X_subset = X.drop_vars('DateStamp')
    X_subset = X_subset.drop_vars('times_in_avg')
    
    Xa = X_subset.mean(dim='time').compute()

    Xa['DateStamp']=X['DateStamp']
    Xa['times_in_avg']=X['times_in_avg']
    
    if (return_dataset==False):
        Xa.to_netcdf( Bfilo )
        print(f' DATA Written to {Bfilo}', flush=True )
        rcode=1
        return rcode
    else:
        print(f' Would be written to {Bfilo}', flush=True )
        return Xa,X


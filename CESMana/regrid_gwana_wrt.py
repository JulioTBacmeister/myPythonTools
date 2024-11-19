#!/usr/bin/env python

################################################
# New style 
################################################
import sys
import os

workdir_ = '/glade/work/juliob'
if ( workdir_ not in sys.path ):
    sys.path.append(workdir_)
    print( f" a path to {workdir_} added in {__name__} {utils_path} ")

from myPythonTools.Utils import utils as uti
from myPythonTools.Utils import numerical_utils as nuti
from myPythonTools.Utils import AveragingUtils as Av
from myPythonTools.Utils import validation_data as Val
from myPythonTools.Utils import PlotUtil as Pu
from myPythonTools.Plotting import LatLonMaps as LL 
from myPythonTools.CASutils import filter_utils as fu 

from PyRegridding.Utils import GridUtils as GrU
from PyRegridding.Utils import MakePressures as MkP
from PyRegridding.Drivers import RegridField as RgF


# The usual
from datetime import date
import numpy as np
import xarray as xr

# Some other useful packages 
import copy
import time
import cftime
import yaml
#from box import Box




#def procField( FLD = None ):
#    global A, RgOb_ne240_x_ne16, RgOb_ne16_x_ne240  # Declare as global to make accessible in procField()
def procField(A, RgOb_ne240_x_ne16, RgOb_ne16_x_ne240, FLD=None):

    ###
    tootrange = f"{(A.X.time[0].values.item() ).strftime('%Y-%m-%d')}_{(A.X.time[-1].values.item() ).strftime('%Y-%m-%d')}"
    fname=f'{A.exp}_GWana_{FLD}_{tootrange}.nc'
    Bdiro=f'/glade/derecho/scratch/juliob/archive/{A.exp}/atm/GWana/'
    #######
    
    fname = f'{Bdiro}{fname}'
    print(f" 'bout to proc {FLD} into: \n{fname}")
    
    
    os.makedirs( Bdiro , exist_ok=True )
    fld = FLD.lower()
    fldO = f'{fld}O'
    fldOx2 = f'{fld}Ox2'
    fldOx2xO = f'{fld}Ox2xO'
    
    Topo=xr.open_dataset(A.X.topography_file)
    topoO=Topo.PHIS.values / 9.8
    
    ##########################################
    # Regrid 'fld' and calculate residual
    # from 2x2 (ne16)
    ##########################################
    ogO=A.X[FLD].values
    long_name = A.X[FLD].long_name
    units = A.X[FLD].units

    print( f'FLD={FLD} {long_name} {units}', flush=True )
    
    ogOx2=RgF.Horz(xfld_Src=ogO , Src='ne240pg3', Dst='ne16pg3' , RegridObj_In=  RgOb_ne240_x_ne16  ) 
    
    ogOx2xO=RgF.Horz(xfld_Src=ogOx2 , Src='ne16pg3' , Dst='ne240pg3', RegridObj_In= RgOb_ne16_x_ne240  ) 
    
    latO,lonO = GrU.latlon( grid='ne240pg3' )
    lat2,lon2 = GrU.latlon( grid='ne16pg3' )
    lev=A.X.lev.values
    
    
    dims   = ["ncol","ncol2","time","lev"]
    coords = dict( 
        lon  = ( ["ncol"],lonO),
        lat  = ( ["ncol"],latO ),
        lon2  = ( ["ncol2"],lon2),
        lat2  = ( ["ncol2"],lat2 ),
        lev  = ( ["lev"],lev),
        time = ( ["time"], A.X.time.values  ), #pd.to_datetime( pdTime_ERA[itim] ) ),
    )
    
    Xw = xr.Dataset( coords=coords  )
    Dar = xr.DataArray( data=topoO, 
                        dims=('ncol',),
                        attrs=dict( long_name='Topography',units='m',) ,) 
    Xw['topoO'] = Dar

    if (A.X[FLD].dims == ('time', 'lev','ncol')):
        Dar = xr.DataArray( data=ogO, 
                            dims=('time','lev','ncol',),
                            attrs=dict( long_name=long_name,units=units ,) ,) 
        Xw[fldO] = Dar
        
        Dar = xr.DataArray( data=ogOx2, 
                            dims=('time','lev','ncol2',),
                            attrs=dict( long_name=f'{long_name} coarsegrained' ,units=units ,) ,) 
        Xw[fldOx2] = Dar
        
        Dar = xr.DataArray( data=ogOx2xO, 
                            dims=('time','lev','ncol',),
                            attrs=dict( long_name=f'{long_name} coarsegrained-prolonged' ,units=units ,) ,) 
        Xw[fldOx2xO] = Dar
        
    elif (A.X[FLD].dims == ('time','ncol')):
        Dar = xr.DataArray( data=ogO, 
                            dims=('time','ncol',),
                            attrs=dict( long_name=long_name,units=units ,) ,) 
        Xw[fldO] = Dar
        
        Dar = xr.DataArray( data=ogOx2, 
                            dims=('time','ncol2',),
                            attrs=dict( long_name=f'{long_name} coarsegrained' ,units=units ,) ,) 
        Xw[fldOx2] = Dar
        
        Dar = xr.DataArray( data=ogOx2xO, 
                            dims=('time','ncol',),
                            attrs=dict( long_name=f'{long_name} coarsegrained-prolonged' ,units=units ,) ,) 
        Xw[fldOx2xO] = Dar
    
    Xw.to_netcdf( fname )
        

#def driver():
#    global A, RgOb_ne240_x_ne16, RgOb_ne16_x_ne240  # Declare as global to make accessible in procField()
def driver():

    with open('configure_gwana.yaml', 'r') as file:
        cfg = yaml.safe_load(file)
    
    print( f"In regrid_gwana_wrt", flush=True )
    ######################################################
    
    x='A'
    exp, subd, Src, Hkey, Dst, useri = cfg[x]['name'] , cfg[x]['subdir'] , cfg[x]['SrcGrid'] , cfg[x]['Hkey'] , cfg[x]['DstGrid'] , cfg[x]['user'] 
    ymdPatLs = cfg[x]['ymdPat']
    hsPat = cfg[x]['hsPat']
    print( exp, subd, Src, Hkey, Dst, useri , flush=True )
    print( ymdPatLs , flush=True )

    created_RegridObjs = False
    
    for ymdPat in ymdPatLs:
        print( f"\n \n \n ################## \n Gonna do - {ymdPat}" , flush=True )
    
        A = uti.MakeDict4Exp( exp=exp , user=useri, subd=subd , 
                             hsPat=hsPat , ymdPat=ymdPat,verbose=True, open_dataset=True )

        target_plevs = [ 932.0, 856., 499., 227.0,  58., 2.8]
        plevs = A.X.lev.values
        lev_sel = uti.find_nearest_plev_indices(plev=plevs , target_levels=target_plevs )
        A.X = A.X.isel( lev=lev_sel )
        print( f'Selected out levels {A.X.lev.values}' , flush=True )
    

        if (created_RegridObjs == False ):
            ##########################
            # Create regridding objects
            #############################
            RgOb_ne240_x_ne16 = RgF.Horz( Src='ne240pg3', Dst='ne16pg3' ) 
            
            RgOb_ne16_x_ne240 = RgF.Horz(  Src='ne16pg3' , Dst='ne240pg3' , RegridMethod='BILINEAR' ) 
            
            RgOb_ne240_x_llOxO = RgF.Horz( Src='ne240pg3', Dst='latlonOxO' ) 
            created_RegridObjs = True 
        else:
            print( f"Using existing Regrid Objects ", flush=True )
    
        # Call procField with required arguments
        if ('TS' in A.X ):
            procField(A, RgOb_ne240_x_ne16, RgOb_ne16_x_ne240, FLD='TS')
        procField(A, RgOb_ne240_x_ne16, RgOb_ne16_x_ne240, FLD='OMEGA')
        procField(A, RgOb_ne240_x_ne16, RgOb_ne16_x_ne240, FLD='U')
        procField(A, RgOb_ne240_x_ne16, RgOb_ne16_x_ne240, FLD='V')

if __name__ == "__main__":
    ################################################

    driver()
    

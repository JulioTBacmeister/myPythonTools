#!/usr/bin/env python
################################################
# New style 
# ###############################################
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
def writeFile(A, zeta_x1, zeta2_x1):

    FLD='zeta'
    ###
    tootrange = f"{(A.X.time[0].values.item() ).strftime('%Y-%m-%d')}_{(A.X.time[-1].values.item() ).strftime('%Y-%m-%d')}"
    fname=f'{A.exp}_GWana_{FLD}_{tootrange}.nc'
    Bdiro=f'/glade/derecho/scratch/juliob/archive/{A.exp}/atm/GWana/'
    #######
    
    fname = f'{Bdiro}{fname}'
    print(f" 'bout to proc {FLD} into: \n{fname}")
    
    
    os.makedirs( Bdiro , exist_ok=True )
    
    lat1,lon1 = GrU.latlon( grid='fv1x1' )
    lev=A.X.lev.values

    dims   = ["lon","lat","time","lev","ilev"]
    coords = dict( 
        lon  = ( ["lon"],lon1 ),
        lat  = ( ["lat"],lat1 ),
        lev  = ( ["lev"],lev),
        time = ( ["time"],  A.X.time.values ) )
                
    Xout = xr.Dataset( coords=coords  )
    
    Dar = xr.DataArray( data=zeta_x1 , dims=('time','lev','lat','lon',),
                                attrs=dict( description='relative_vorticity',units='s-1',) ,) 
    Xout['zeta'] = Dar

    Dar = xr.DataArray( data=zeta2_x1 , dims=('time','lev','lat','lon',),
                                attrs=dict( description='enstrophy',units='s-2',) ,) 
    Xout['zeta2'] = Dar

    Xout.to_netcdf( fname )
        

#def driver():
#    global A, RgOb_ne240_x_ne16, RgOb_ne16_x_ne240  # Declare as global to make accessible in procField()
def vorticity(write_to_file=True):

    with open('configure_gwana_genl.yaml', 'r') as file:
        cfg = yaml.safe_load(file)
    
    print( f"In regrid_gwana_wrt", flush=True )
    ######################################################
    
    x='A_xy_vort'
    exp, subd, Src, Hkey, Dst, useri = cfg[x]['name'] , cfg[x]['subdir'] , cfg[x]['SrcGrid'] , cfg[x]['Hkey'] , cfg[x]['DstGrid'] , cfg[x]['user'] 
    ymdPatLs = cfg[x]['ymdPat']
    hsPat = cfg[x]['hsPat']
    print( exp, subd, Src, Hkey, Dst, useri , flush=True )
    print( ymdPatLs , flush=True )

    #created_RegridObjs = False
    #####################################
    # Initialize regrid-object library
    RgObLib={}
    RgOb_ne240_x_llOxO = GrU.regrid_object_lib(RgOb=RgObLib, src='ne240pg3', dst='latlonOxO', RegridMethod='BILINEAR' ) 
    RgOb_llOxO_x_fv1x1 = GrU.regrid_object_lib(RgOb=RgObLib, src='latlonOxO', dst='fv1x1', RegridMethod='BILINEAR' ) 
    
    for ymdPat in ymdPatLs:
        print( f"\n ################## \n Gonna do - {ymdPat}" , flush=True )
    
        A = uti.MakeDict4Exp( exp=exp , user=useri, subd=subd , 
                             hsPat=hsPat , ymdPat=ymdPat,verbose=True, open_dataset=True )

        target_plevs = [ 932.0, 856., 499., 227.0,  58., 2.8]
        plevs = A.X.lev.values
        lev_sel = uti.find_nearest_plev_indices(plev=plevs , target_levels=target_plevs )
        A.X = A.X.isel( lev=lev_sel )
        print( f'Selected out levels {A.X.lev.values}' , flush=True )

        uu=A.X.U.values
        print( f'Extracted U' , flush=True )
        vv=A.X.V.values
        print( f'Extracted V' , flush=True )

        #uu_xO=RgF.Horz(xfld_Src=uu , Src='ne240pg3', Dst='latlonOxO' , RegridObj_In=  RgOb_ne240_x_llOxO  ) 
        #vv_xO=RgF.Horz(xfld_Src=vv , Src='ne240pg3', Dst='latlonOxO' , RegridObj_In=  RgOb_ne240_x_llOxO ) 
        uu_xO=RgF.Horz(xfld_Src=uu , RegridObj_In=  RgOb_ne240_x_llOxO  ) 
        vv_xO=RgF.Horz(xfld_Src=vv , RegridObj_In=  RgOb_ne240_x_llOxO ) 
        lat_xO,lon_xO =GrU.latlon(grid='latlonOxO' )
        
        nt,nz,ny,nx=np.shape( uu_xO )
        zeta_xO=np.zeros((nt,nz,ny,nx) )
        for t in np.arange( nt ):
            for z in np.arange(nz):
                if ( (t%10)==0):
                    print( f'time {t}' ) 
                zeta_xO [t,z,:,:] = nuti.Sphere_Curl2( uu_xO[t,z,:,:] , vv_xO[t,z,:,:] , lat_xO, lon_xO, wrap=True )

        zeta2_xO = zeta_xO**2
        print( f'Squared zeta' , flush=True )

        #zeta_xO_x1  = RgF.Horz(xfld_Src=zeta_xO , Src='latlonOxO' , Dst='fv1x1', RegridObj_In=  RgOb_llOxO_x_fv1x1  ) 
        #zeta2_xO_x1 = RgF.Horz(xfld_Src=zeta2_xO , Src='latlonOxO' , Dst='fv1x1', RegridObj_In=  RgOb_llOxO_x_fv1x1  ) 
        zeta_xO_x1  = RgF.Horz(xfld_Src=zeta_xO , RegridObj_In=  RgOb_llOxO_x_fv1x1  ) 
        zeta2_xO_x1 = RgF.Horz(xfld_Src=zeta2_xO, RegridObj_In=  RgOb_llOxO_x_fv1x1  ) 

        if (write_to_file==False):
            return zeta_xO,zeta_xO_x1,zeta2_xO_x1
        else:            
            writeFile(A, zeta_xO_x1, zeta2_xO_x1)
        

if __name__ == "__main__":
    ################################################

    vorticity()
    

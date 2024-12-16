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
#def driver():
#    global A, RgOb_ne240_x_ne16, RgOb_ne16_x_ne240  # Declare as global to make accessible in procField()
def PBLq(write_to_file=True):

    with open('configure_gwana_genl.yaml', 'r') as file:
        cfg = yaml.safe_load(file)
    
    print( f"In regrid_gwana_wrt", flush=True )
    ######################################################
    
    x='A_xy_PBLq'
    exp, subd, Src, Hkey, Dst, useri = cfg[x]['name'] , cfg[x]['subdir'] , cfg[x]['SrcGrid'] , cfg[x]['Hkey'] , cfg[x]['DstGrid'] , cfg[x]['user'] 
    ymdPatLs = cfg[x]['ymdPat']
    hsPat = cfg[x]['hsPat']
    print( exp, subd, Src, Hkey, Dst, useri , flush=True )
    print( ymdPatLs , flush=True )

    # process an estimate of TS variance ... ot NOT
    do_var_ts=True

    #created_RegridObjs = False
    #####################################
    # Initialize regrid-object library
    RgObLib={}
    RgOb_ne240_x_fv1x1 = GrU.regrid_object_lib(RgOb=RgObLib, src='ne240pg3', dst='fv1x1') 
    RgOb_fv1x1_x_ne240 = GrU.regrid_object_lib(RgOb=RgObLib, src='fv1x1',dst='ne240pg3',RegridMethod='BILINEAR') 




    for ymdPat in ymdPatLs:
        print( f"\n ################## \n Gonna do - {ymdPat}" , flush=True )
        A = uti.MakeDict4Exp( exp=exp , user=useri, subd=subd , 
                             hsPat='cam.h3i' , ymdPat=ymdPat,verbose=True, open_dataset=True )
        B = uti.MakeDict4Exp( exp=exp , user=useri, subd=subd , 
                             hsPat='cam.h2i' , ymdPat=ymdPat,verbose=True, open_dataset=True )
    
        target_plevs = [ 990., 950., 900., 850.,800.,700.]  
        plevs = A.X.ilev.values
        ilev_sel = uti.find_nearest_plev_indices(plev=plevs , target_levels=target_plevs )
        A.X = A.X.isel( ilev=ilev_sel )
        plevs = A.X.lev.values
        lev_sel = uti.find_nearest_plev_indices(plev=plevs , target_levels=target_plevs )
        A.X = A.X.isel( lev=lev_sel )
        
        print( f'Selected out levels {A.X.lev.values} and ilevels {A.X.ilev.values}' , flush=True )
    
        B.X = B.X.sel( time= A.X.time )

        # 'PS', 'PTTEND', 'THLP2_CLUBB', 'UPWP_CLUBB', 'VPWP_CLUBB', 'WP2_CLUBB', 
        T = xr.open_dataset( A.X.topography_file )
        
        thlp2=A.X.THLP2_CLUBB.values
        print( f'Extracted THLP2' , flush=True )
        wp2=A.X.WP2_CLUBB.values
        print( f'Extracted WP2' , flush=True )
        upwp=A.X.UPWP_CLUBB.values
        print( f'Extracted UPWP' , flush=True )
        vpwp=A.X.VPWP_CLUBB.values
        print( f'Extracted VPWP' , flush=True )
        pttend=A.X.PTTEND.values
        print( f'Extracted PTTEND' , flush=True )

        shflx=B.X.SHFLX.values
        print( f'Extracted SHFLX' , flush=True )
        lhflx=B.X.LHFLX.values
        print( f'Extracted LHFLX' , flush=True )
        precc=B.X.PRECC.values
        print( f'Extracted PRECC' , flush=True )
        precl=B.X.PRECL.values
        print( f'Extracted PRECL' , flush=True )

        ts=B.X.TS.values
        print( f'Extracted TS' , flush=True )

        phis=T.PHIS.values
        sgh=T.SGH.values
        sgh30=T.SGH30.values
        print( f'Extracted Topo fields', flush=True )
        

        wp2_x1=RgF.Horz(xfld_Src=wp2 , RegridObj_In=  RgOb_ne240_x_fv1x1  ) 
        thlp2_x1=RgF.Horz(xfld_Src=thlp2 , RegridObj_In=  RgOb_ne240_x_fv1x1 ) 
        upwp_x1=RgF.Horz(xfld_Src=upwp , RegridObj_In=  RgOb_ne240_x_fv1x1 ) 
        vpwp_x1=RgF.Horz(xfld_Src=upwp , RegridObj_In=  RgOb_ne240_x_fv1x1 ) 
        pttend_x1=RgF.Horz(xfld_Src=pttend , RegridObj_In=  RgOb_ne240_x_fv1x1 ) 

        shflx_x1=RgF.Horz(xfld_Src=shflx , RegridObj_In=  RgOb_ne240_x_fv1x1 ) 
        lhflx_x1=RgF.Horz(xfld_Src=lhflx , RegridObj_In=  RgOb_ne240_x_fv1x1 ) 
        precc_x1=RgF.Horz(xfld_Src=precc , RegridObj_In=  RgOb_ne240_x_fv1x1 ) 
        precl_x1=RgF.Horz(xfld_Src=precl , RegridObj_In=  RgOb_ne240_x_fv1x1 ) 
        ts_x1=RgF.Horz(xfld_Src=ts , RegridObj_In=  RgOb_ne240_x_fv1x1 ) 

        phis_x1=RgF.Horz(xfld_Src=phis , RegridObj_In=  RgOb_ne240_x_fv1x1 ) 
        sgh_x1=RgF.Horz(xfld_Src=sgh , RegridObj_In=  RgOb_ne240_x_fv1x1 ) 
        sgh30_x1=RgF.Horz(xfld_Src=sgh30 , RegridObj_In=  RgOb_ne240_x_fv1x1 ) 



        # Var TS calcukatuin
        if (do_var_ts==True ):
            ts_x1_cO = RgF.Horz(xfld_Src=ts_x1 , RegridObj_In=  RgOb_fv1x1_x_ne240) 
            dts = ts - ts_x1_cO
            dts2 = dts**2
            dts2_x1 = RgF.Horz(xfld_Src=dts2 , RegridObj_In=  RgOb_ne240_x_fv1x1 )
        else:
            dts2_x1=None


        

        if (write_to_file==False):
            return wp2_x1
        else:            
            writeFile(A, thlp2_x1, wp2_x1, upwp_x1, vpwp_x1, pttend_x1,
                     shflx_x1, lhflx_x1, precc_x1, precl_x1, ts_x1, dts2_x1,
                     phis_x1, sgh_x1, sgh30_x1 )


def writeFile(A, thlp2, wp2, upwp, vpwp, pttend, 
              shflx, lhflx, precc, precl, ts , dts2,
              phis, sgh, sgh30 ):

    FLD='PBLq'
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
    ilev=A.X.ilev.values

    dims   = ["lon","lat","time","lev","ilev"]
    coords = dict( 
        lon  = ( ["lon"],lon1 ),
        lat  = ( ["lat"],lat1 ),
        lev  = ( ["lev"],lev),
        ilev  = ( ["ilev"],ilev),
        time = ( ["time"],  A.X.time.values ) )
                
    Xout = xr.Dataset( coords=coords  )
    
    Dar = xr.DataArray( data=thlp2 , dims=('time','ilev','lat','lon',),
                                attrs=dict( description='thlp_var',units='K+2',) ,) 
    Xout['thlp2'] = Dar

    Dar = xr.DataArray( data=wp2 , dims=('time','ilev','lat','lon',),
                                attrs=dict( description='wp_var',units='m+2 s-2',) ,) 
    Xout['wp2'] = Dar

    Dar = xr.DataArray( data=upwp , dims=('time','ilev','lat','lon',),
                                attrs=dict( description='momflux_x',units='m+2 s-2',) ,) 
    Xout['upwp'] = Dar

    Dar = xr.DataArray( data=vpwp , dims=('time','ilev','lat','lon',),
                                attrs=dict( description='momflux_y',units='m+2 s-2',) ,) 
    Xout['vpwp'] = Dar

    Dar = xr.DataArray( data=pttend , dims=('time','lev','lat','lon',),
                                attrs=dict( description='total_Phys_tend',units='K s-1',) ,) 
    Xout['pttend'] = Dar

    Dar = xr.DataArray( data=shflx , dims=('time','lat','lon',),
                                attrs=dict( description='sens_heat_flux',units='W m-2',) ,) 
    Xout['shflx'] = Dar

    Dar = xr.DataArray( data=lhflx , dims=('time','lat','lon',),
                                attrs=dict( description='latent_heat_flux',units='W m-2',) ,) 
    Xout['lhflx'] = Dar

    Dar = xr.DataArray( data=precc , dims=('time','lat','lon',),
                                attrs=dict( description='conv_prec',units='m s-1',) ,) 
    Xout['precc'] = Dar

    Dar = xr.DataArray( data=precl , dims=('time','lat','lon',),
                                attrs=dict( description='large_scale_prec',units='m s-1',) ,) 
    Xout['precl'] = Dar

    Dar = xr.DataArray( data=ts , dims=('time','lat','lon',),
                                attrs=dict( description='surface_temp',units='K',) ,) 
    Xout['ts'] = Dar

    Dar = xr.DataArray( data=phis , dims=('lat','lon',),
                                attrs=dict( description='surface_geopot_ht',units='m+2',) ,) 
    Xout['phis'] = Dar

    Dar = xr.DataArray( data=sgh , dims=('lat','lon',),
                                attrs=dict( description='topo_variance',units='m+2',) ,) 
    Xout['sgh'] = Dar

    if (dts2 is not None):
        Dar = xr.DataArray( data=dts2 , dims=('time','lat','lon',),
                                    attrs=dict( description='surface_temp_var',units='K+2',) ,) 
        Xout['var_ts'] = Dar


    Xout.to_netcdf( fname )
        

        

if __name__ == "__main__":
    ################################################

    PBLq()
    

#workdir_ = '/glade/work/juliob/'
import sys
#######################################
# Leave this for now. But it should change to better
# method as here:
import os
This_module_path = os.getcwd()  #os.path.dirname(os.path.abspath(__file__))
workdir_ = os.path.join(This_module_path, '../../' )

# The usual
from datetime import date
import numpy as np
import xarray as xr

# Cartopy for pretty maps
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Some other useful packages 
import importlib
import copy
import time
import cftime

########################################
# Now import your own stuff
########################################

#From in here
import gw_common as GWc
import gw_movmtn as GWmm

sys.path.append(workdir_ + 'myPythonTools/Utils/')
import utils as uti
import numerical_utils as nuti
import AveragingUtils as Av
import shr_const
from shr_const import ShrConst as Cs

sys.path.append(workdir_ + 'PyRegridding/Regridder/')
import var_A_x_B as vAB

import yaml



def tau_prof( date=None,year=None,month=None,day=None,hour=None, SourceMethod="vort500"):

    if (date != None):
        year,month,day,hour = date
    
    seconds = hour * 3600
    timetag = f'{year:04d}-{month:02d}-{day:02d}-{seconds:05d}'
    
    #----------------------------------------------------
    #  Get CAM experiment
    #-----------------------------------------------------
    with open('gw_driver.yaml', 'r') as file:
        cfg = yaml.safe_load(file)

    expA = cfg['A1']['casename'] 
    subd = cfg['A1']['grid'] 
    user = cfg['A1']['user'] 
    A1 = uti.MakeDict4Exp( exp=expA ,
                         user=user ,
                         subd=subd , 
                         hsPat=cfg['A1']['hsPat'],
                         ymdPat= timetag ,
                         verbose=True, open_dataset=True )
    A3 = uti.MakeDict4Exp( exp=expA ,
                         user=user ,
                         subd=subd ,  
                         hsPat=cfg['A3']['hsPat'],
                         ymdPat= timetag ,
                         verbose=True, open_dataset=True )


    #####################################################################
    # Before subselecting to SH high lats:
    # Need to do some calcualtions on full global grid for new sources
    #####################################################################
    if (SourceMethod=="vort500"):
        k5000=68 # ~500hPa surface in L93
        u_se = A1.X.UEGW[0,k5000,:].values
        v_se = A1.X.VEGW[0,k5000,:].values        
        u_fv,lat_fv,lon_fv = vAB.Hregrid( avar=u_se, agrid='ne30pg3', bgrid='fv1x1' )
        v_fv,lat_fv,lon_fv = vAB.Hregrid( avar=v_se, agrid='ne30pg3', bgrid='fv1x1' )
        zeta_fv = nuti.Sphere_Curl2( u_fv , v_fv , lat_fv, lon_fv, wrap=True )
        zeta_se,lat_se,lon_se = vAB.Hregrid( avar=zeta_fv , agrid='fv1x1' ,bgrid='ne30pg3')
        zeta_se = np.expand_dims( zeta_se, axis=0 )
        print("Shape of zeta_se:",np.shape(zeta_se))
        Dar = xr.DataArray( data=zeta_se, 
                        dims=('time','ncol',),
                        attrs=dict( description='Vorticity near 500hPa',units='Pa',) ,) 
        A1.X['ZETA500'] = Dar
 
    #####################################################################
    # Subselecting SH high-lats
    #####################################################################
    latitude_range = (-90, -45)  # Example latitude range
    
    latitude_mask = np.where(
        (A1.X['lat'].values >= latitude_range[0]) & (A1.X['lat'].values <= latitude_range[1]) )
    
    
    print( len(latitude_mask[0]))
    # Select data using these indices
    #selected_data = ds.isel(ncol=latitude_indices)
    
    latitude_indices=latitude_mask[0]
    
    X1 = A1.X.isel(ncol=latitude_indices )
    X3 = A3.X.isel(ncol=latitude_indices )


    print( np.shape(A1.X.UEGW))
    print( np.shape(X1.UEGW))

    


    #######################
    # get constants to use here
    ######################
    cpair=Cs.CPDAIR

    pifc = X1.PINT[0,:,:].values.T
    te = X1.TEGW[0,:,:].values.T
    uu = X1.UEGW[0,:,:].values.T
    vv = X1.VEGW[0,:,:].values.T
    zm = X3.Z3[0,:,:].values.T

    if (SourceMethod=="SavedCLUBBmomflux"):
        xpwp = X3.XPWP_SRC_MOVMTN[0,:].values.T
        taudesc="TAU_xpwp_SV_b"
        
    elif (SourceMethod=="vort500"):
        taudesc="TAU_vort500"
        xpwp = 100.*np.abs(X1.ZETA500[0,:].values.T)
        
    elif (SourceMethod=="uniform"):
        taudesc="TAU_uniform"
        xpwp = 0.0*X3.XPWP_SRC_MOVMTN[0,:].values.T + 0.005
        
    else:
        xpwp = X3.XPWP_SRC_MOVMTN[0,:].values.T
        taudesc="TAU_xpwp_SV"
    

    #######################
    # Make derived types
    ######################
    band_movmtn = GWc.BandType( ngwv=0, dc=5., kwv=2*np.pi/(100_000.) )
    PP=GWc.PType(ifc=pifc)
   
    # get dimensions
    ncol,pver =np.shape( te )
    ngwv = band_movmtn.ngwv
    pcnst=1

    rhoi, nm, ni =GWc.gw_prof(ncol, PP, cpair, te)


    ####
    netdt=np.zeros( (ncol,pver) )
    netdt_shcu=np.zeros( (ncol,pver) )
    src_level = np.zeros(ncol, dtype=int)
    tend_level = np.zeros(ncol, dtype=int)
    ubi=np.zeros( (ncol,pver+1) )
    ubm=np.zeros( (ncol,pver) )
    tau=np.zeros( (ncol, 2*ngwv+1, pver+1) )
    c = np.zeros((ncol, 2*ngwv+1) )
    xv= np.zeros( ncol )
    yv= np.zeros( ncol )
    hdepth= np.zeros( ncol )
    
    aack= GWmm.gw_movmtn_src(ncol , band_movmtn , uu, vv, netdt, netdt_shcu, xpwp , zm, src_level, tend_level, tau, ubm, ubi, xv, yv, c, hdepth)

    ###########################################################
    # Complete init of arguments to gw_drag_prof 
    # '#'-commented indicates this was init somewhere above
    ###########################################################
        
    # Initialize the real arrays (using float64 for real(r8))
    dt = 1.0
    #te = np.zeros((ncol, pver))
    piln = np.log( PP.ifc )  # Log of interface pressures.
    effgw = np.zeros(ncol) + 1.0
    
    kvtt = np.zeros((ncol, pver+1))
    
    # Assuming q has some third dimension, let's set it to 5 for now
    q = np.zeros((ncol, pver, pcnst))
    dse = np.zeros((ncol, pver))
    vramp = np.zeros(1)  # Assuming vramp is a pointer to a single-element array
    
    #tau = np.zeros((ncol, 2*ngwv+1, pver+1))
    utgw = np.zeros((ncol, pver))
    vtgw = np.zeros((ncol, pver))
    ttgw = np.zeros((ncol, pver))
    
    # Assuming qtgw has the same third dimension as q
    qtgw = np.zeros((ncol, pver, pcnst))
    egwdffi = np.zeros((ncol, pver+1))
    gwut = np.zeros((ncol, pver, 2*ngwv+1))
    
    dttdf = np.zeros((ncol, pver))
    dttke = np.zeros((ncol, pver))
    
    # Optional arguments
    ro_adjust = None # np.zeros((ncol, 2*ngwv+1, pver+1))
    kwvrdg = None #np.zeros(ncol)
    satfac_in = 2.0  # Assuming 2.0 for backward compatibility
    lapply_effgw_in = False
    lapply_vdiff = False
    tau_diag = np.zeros((ncol, pver+1))
    
    # Example function call
    gw_calc = GWc.gw_drag_prof(ncol, band_movmtn, PP , src_level, tend_level, dt, 
                    te, vramp,   
                    piln, rhoi,    nm,   ni,  ubm,  ubi,  xv,    yv,   
                    effgw,      c, kvtt, q,   dse,  tau,  utgw,  vtgw, 
                    ttgw, qtgw, egwdffi,   gwut, dttdf, dttke, ro_adjust, 
                    kwvrdg, satfac_in, lapply_effgw_in, lapply_vdiff, tau_diag ,
                    perform_second_half=False, tau_0_ubc=False, do_vertical_diffusion=False )


    nt,nz,ncol = np.shape( A1.X.TEGW )
    newtauX = np.zeros( (nt,nz+1,ncol) )
    newtau = np.zeros( (nz+1,ncol) )

    # 3. Reshape B_array to (94, 7020)
    tau_reshaped = np.transpose(np.squeeze(tau, axis=1), (1, 0))

    print('\n')
    print( 'rehsaped tau ', np.shape( tau_reshaped) )
    print( 'newtau', np.shape( newtau) )
    
    #newtau[0,:,latitude_indices] = tau_reshaped
    newtau[:,latitude_indices] = tau_reshaped

    #newtau = np.expand_dims(newtau, axis=0)
    newtauX[0,:,:] = newtau
    print( 'EXPANDED newtau', np.shape( newtauX) )

    write_file(A1=A1, tau=newtauX ,timetag=timetag,taudesc=taudesc)

    return 0

def write_file(A1,tau,timetag,taudesc):
    
    nt,nz,ncol = np.shape( A1.X.TEGW )
    dims   = ["ncol","time","lev","ilev","nbnd"]
    
    coords = dict( 
        lon  = ( ["ncol"], A1.X.lon.values  ),
        lat  = ( ["ncol"], A1.X.lat.values  ),
        lev  = ( ["lev"],  A1.X.lev.values  ),
        ilev = ( ["ilev"], A1.X.ilev.values ),
        time = ( ["time"], A1.X.time.values ), 
        nbnd = ( ["nbnd"], [0,1] ), 
        )
    
    print( '\n' )
    print( nt,nz,ncol )
    
    Wds = xr.Dataset( coords=coords  )


    Wds['date']    = A1.X.date
    Wds['datesec'] = A1.X.datesec

    if ( 'time_bnds' in A1.X ):
        Wds['time_bnds'] = A1.X['time_bnds']
    if ( 'time_bounds' in A1.X ):
        Wds['time_bounds'] = A1.X['time_bounds']

    Dar = xr.DataArray( data=A1.X.hyai.values, dims=('ilev',),
                        attrs=dict( description='interface hybrid eta coordinate A-coeff ',units='1',) ,) 
    Wds['hyai'] = Dar

    Dar = xr.DataArray( data=A1.X.hybi.values, dims=('ilev',),
                        attrs=dict( description='interface hybrid eta coordinate B-coeff ',units='1',) ,) 
    Wds['hybi'] = Dar

    Dar = xr.DataArray( data=A1.X.hyam.values, dims=('lev',),
                        attrs=dict( description='mid-layer hybrid eta coordinate A-coeff ',units='1',) ,) 
    Wds['hyam'] = Dar

    Dar = xr.DataArray( data=A1.X.hybm.values, dims=('lev',),
                        attrs=dict( description='mid-layer hybrid eta coordinate B-coeff ',units='1',) ,) 
    Wds['hybm'] = Dar

    Dar = xr.DataArray( data=A1.X.area.values, dims=('ncol',),
                        attrs=dict( description='Cell area',units='Steradians',) ,) 
    Wds['area'] = Dar
    
    Dar = xr.DataArray( data=tau, 
                        dims=('time','ilev','ncol',),
                        attrs=dict( description='GW momentum flux',units='Pa',) ,) 
    Wds['TAU_redo'] = Dar




    C = uti.MakeDict4Exp( exp=A1.exp ,
                     user=A1.user,
                     subd=A1.subd, 
                     hsPat=f'cam.h1i.{taudesc}',
                     ymdPat= timetag ,
                     verbose=True, open_dataset=False )

    Wds.to_netcdf( C.path )






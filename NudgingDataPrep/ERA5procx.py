# Import packages 
import sys

import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.tri as tri
import matplotlib.colors as colors

from scipy.io import FortranFile
from scipy import interpolate as intr

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import ESMF as E

import importlib
import glob
import copy
import time

import scripGen as SG
import esmfRegrid as erg
import VertRegrid as vrg

# import modules in other directories
sys.path.append('../Utils/')
import MakePressures as MkP
import humiditycalcs as hum

# Reload local packages that are under
# development
importlib.reload( erg )
importlib.reload( vrg )
importlib.reload( SG )
importlib.reload( MkP )
importlib.reload( hum )

# Physical Constants
Rgas = 287.0 # J K-1 kg-1


#-------------------------------------------------------------
#  Naming conventions
#-------------------------------------------------------------
# aaa_{CAM,ERA} 
# Indicates the immediate provenance of a variable, e.g.,
#      phis_CAM ==> phis from CAM on CAM's grid
#      phis_ERA ==> phis from ERA on ERA's grid
# lower case 'phis' indicates this is ndarray-like 
#
# aaa_{CAM,ERA}_x{ERA,CAM}
# Indicates variable has been remapped horizontall. So, e.g.
#      phis_ERA_xCAM ==> ERA phis remapped horizontally to the CAM grid 
#
# aaa_{CAM,ERA}_xz{ERA,CAM}
# Indicates variable has been remapped horizontally AND vertically. So, e.g.
#      te_ERA_xzCAM ==> ERA temperature remapped horizontally to the CAM horizontal grid 
#                       and then also vertically interpoated to the CAM vertical grid
#-------------------------------------------------------------

def prep():
    #---------------------------------------------
    # This function sets-up variables and objects 
    # that are need for horizontal and vertical 
    # regridding of ERA reanalyses.
    #---------------------------------------------

    global regrd,srcf,dstf
    global phis_CAM, phis_ERA
    global amid_CAM,bmid_CAM,aint_CAM,bint_CAM
    global lon_CAM,lat_CAM,area_CAM

    # Grid keys for remapping
    global srcHkey,dstHkey,srcTHkey,dstTHkey,srcZHkey,dstZHkey,srcTZHkey,dstTZHkey


    #------- 
    # Begin
    #-------

    tic = time.perf_counter()

    # Set grid keys for Src ERA5 renalysis
    srcHkey = 'yx'
    srcTHkey  = 't'  + srcHkey
    srcZHkey  = 'z'  + srcHkey
    srcTZHkey = 'tz' + srcHkey

    # Set grid keys for Dst CAM-SE
    dstHkey   = 'c'
    dstTHkey  = 't'  + dstHkey
    dstZHkey  = 'z'  + dstHkey
    dstTZHkey = 'tz' + dstHkey


 
    # Get all the met and topo data we will use
    bnd_topo = \
               '/glade/p/cgd/amp/juliob/bndtopo/latest/ne30pg3_gmted2010_modis_bedmachine_nc3000_Laplace0100_20230105.nc'
    dsTopo_CAM=xr.open_dataset( bnd_topo )
    phis_CAM = dsTopo_CAM['PHIS'].values
    lon_CAM  = dsTopo_CAM['lon'].values
    lat_CAM  = dsTopo_CAM['lat'].values
    area_CAM = dsTopo_CAM['area'].values

    print(  dsTopo_CAM['PHIS'].dims )

    # Read in ERA5 topography
    eratopof='/glade/work/juliob/ERA5-proc/ERA5interp/phis/ERA5_phis.nc'
    dsTopo_ERA=xr.open_dataset( eratopof )
    phis_ERA=dsTopo_ERA['Z_GDS4_SFC'].values


    # ----------------------------------------------
    #  Setp regridding machinery
    # ----------------------------------------------

    # Scrip file for ERA5 created by ERA5scrip.ipynb
    ERA5scrip = '/glade/work/juliob/ERA5-proc/ERA5interp/grids/ERA5_640x1280_scrip.nc'
    
    src_scrip=ERA5scrip
    src_type='grid'

    # ------- CAM SE ne30pg3 Scrip file
    scripdir='/glade/p/cesmdata/cseg/inputdata/share/scripgrids/'
    ne30scrip  = scripdir +  "ne30pg3_scrip_170611.nc"
    
    dst_scrip=ne30scrip
    dst_type='mesh'
    
    griddir = "/glade/work/juliob/ERA5-proc/ERA5interp/grids/"
    wgts_file_Con = griddir + "ERA5_ne30pg3_Conserv_wgts.nc"
    
    # Make object for Conservative regridding from ERA5
    # grid to CAM target. Scrip files need to be provided even 
    # when a weight file is used
    regrd, srcf, dstf = erg.Regrid( srcScrip = src_scrip , 
                                    srcType  = src_type  ,
                                    dstScrip = dst_scrip ,
                                    dstType  = dst_type  ,
                                    write_weights = False ,
                                    read_weights = True ,
                                    weights_file = wgts_file_Con )
    


    # Read in CAM L58 vertical grid
    vCAMfile = '/glade/work/juliob/ERA5-proc/CAM-grids/Vertical/'+\
               'GRID_48_taperstart10km_lowtop_BL10_v3p1_beta1p75.nc'
    vCAM=xr.open_dataset( vCAMfile )
    amid_CAM = vCAM['hyam'].values
    bmid_CAM = vCAM['hybm'].values
    aint_CAM = vCAM['hyai'].values
    bint_CAM = vCAM['hybi'].values

    toc = time.perf_counter()
    pTime = f"Prepping for ERA5 proc took  {toc - tic:0.4f} seconds"
    print(pTime)
 
    code = 1
    return code

def proc_1_ERA5_time( year=2022, month=11, day=1, hour0=0):

    global pdTime_ERA
    global ps_CAM, phis_ERA_xCAM, phis_CAM
    global te_ERA_xzCAM
    global q_ERA_xzCAM
    global u_ERA_xzCAM
    global v_ERA_xzCAM
    global w_ERA_xzCAM

    hour1=hour0+5

    monStr=str( year ).zfill(4)+str(month).zfill(2)

    ymdh0=str( year ).zfill(4)+str(month).zfill(2)+str(day).zfill(2)+str(hour0).zfill(2)
    ymdh1=str( year ).zfill(4)+str(month).zfill(2)+str(day).zfill(2)+str(hour1).zfill(2)
    ymdh=ymdh0+'_'+ymdh1

    print( "Time tags for ERA5 files ...")
    print(monStr)
    print(ymdh) 

    print(dstTZHkey)

    era5dir = "/glade/collections/rda/data/ds633.6/e5.oper.an.ml/"
    wrkdir=era5dir+monStr+"/"

    tic = time.perf_counter()

    spfile= wrkdir + 'e5.oper.an.ml.128_134_sp.regn320sc.'+ymdh+'.nc'
    dsPS_ERA  = xr.open_dataset( spfile )
    ps_ERA = dsPS_ERA['SP'].values

    tfile = wrkdir + 'e5.oper.an.ml.0_5_0_0_0_t.regn320sc.'+ymdh+'.nc'
    dsT_ERA   = xr.open_dataset( tfile )
    te_ERA = dsT_ERA['T'].values
    pdTime_ERA = pd.to_datetime( dsT_ERA['time'].values )

    qfile = wrkdir + 'e5.oper.an.ml.0_5_0_1_0_q.regn320sc.'+ymdh+'.nc'
    dsQ_ERA   = xr.open_dataset( qfile )
    q_ERA = dsQ_ERA['Q'].values

    ufile = wrkdir + 'e5.oper.an.ml.0_5_0_2_2_u.regn320uv.'+ymdh+'.nc'
    dsU_ERA   = xr.open_dataset( ufile )
    u_ERA = dsU_ERA['U'].values

    vfile = wrkdir + 'e5.oper.an.ml.0_5_0_2_3_v.regn320uv.'+ymdh+'.nc'
    dsV_ERA   = xr.open_dataset( vfile )
    v_ERA = dsV_ERA['V'].values

    wfile = wrkdir + 'e5.oper.an.ml.0_5_0_2_8_w.regn320sc.'+ymdh+'.nc'
    dsW_ERA   = xr.open_dataset( wfile )
    w_ERA = dsW_ERA['W'].values

    toc = time.perf_counter()
    pTime = f"Reading one set of ERA5 vars took  {toc - tic:0.4f} seconds"
    print(pTime)

    tic = time.perf_counter()
    # Horz remap 
    #########################
    phis_ERA_xCAM = erg.HorzRG( aSrc = phis_ERA , 
                                regrd = regrd , 
                                srcField=srcf , 
                                dstField=dstf , 
                                srcGridkey='yx',  # srcHkey
                                dstGridkey='c' )  # dstHkey
    
    
    #Calculate difference between phis's
    Dphis = phis_ERA_xCAM - phis_CAM
    
    
    
    #Horz remap. 
    #########################
    ps_ERA_xCAM    = erg.HorzRG( aSrc = ps_ERA , 
                                 regrd = regrd , 
                                 srcField=srcf , 
                                 dstField=dstf , 
                                 srcGridkey='tyx', # srcTHkey
                                 dstGridkey='c' )  # dstHkey
    
    
    
    #try function
    #########################
    te_ERA_xCAM    = erg.HorzRG( aSrc = te_ERA , 
                                 regrd = regrd , 
                                 srcField=srcf , 
                                 dstField=dstf , 
                                 srcGridkey='tzyx', # srcTZHkey
                                 dstGridkey='c' )   # dstHkey
    
    
    #-------------------------------------------------------------------
    #                    "CAM surface pressure"
    #-------------------------------------------------------------------
    # We don't actually have ps from CAM, so we make a guess based on
    # ps_ERA_xCAM and te_bot asdescribed above.  In a sense this is a 
    # vertical remapping, so we could call this variable ps_ERA_xzCAM, 
    # but we'll just call it ps_CAM ...
    #-------------------------------------------------------------------
    nt,nz,ncol = np.shape(te_ERA_xCAM)
    L_bot = nz-1
    
    te_bot = te_ERA_xCAM[:,L_bot,:]
    ps_CAM = np.zeros( (nt , ncol) )
    for i in np.arange( nt ):
        ps_CAM[i,:] = ps_ERA_xCAM[i,:] * np.exp( Dphis / (Rgas*te_bot[i,:]) )
        
    #-----------------------------------------------------------------------------------------------
    # Now we creat full 4(3)D pressure fields on the ERA vertical grid. These are "CAM pressures" 
    # as they are based on the "CAM surface pressure" derived above.  These are used for 
    # vertical interpolation below.  I'm not sure about this.  The WO2015 document seems to be
    # instructing this ... but you could also use ERA surface press on xCAM.
    #-----------------------------------------------------------------------------------------------
    amid=dsT_ERA['a_model'].values 
    bmid=dsT_ERA['b_model'].values
    aint=dsT_ERA['a_half'].values 
    bint=dsT_ERA['b_half'].values
    
    p_00 = 1.0 # a_model and a_half appear to include * P_00 in their defintion so set to 1.0 here
    pmid_CAM_zERA, pint_CAM_zERA, delp_CAM_zERA \
        = MkP.Pressure (am=amid ,
                        bm=bmid ,
                        ai=aint ,
                        bi=bint ,
                        ps=ps_CAM ,
                        p_00=p_00 )



    p_00=100_000.0
    pmid_CAM,pint_CAM,delp_CAM \
        = MkP.Pressure (am=amid_CAM ,
                        bm=bmid_CAM ,
                        ai=aint_CAM ,
                        bi=bint_CAM ,
                        ps=ps_CAM ,
                        p_00=p_00 )



    #-----------------------------------------------------------
    # Log-pressure is preferred for vertical interpolation
    # per Williamson&Olson
    #-----------------------------------------------------------
    lnpint_CAM = -7_000. * np.log( pint_CAM / p_00 )
    lnpmid_CAM = -7_000. * np.log( pmid_CAM / p_00 )
    lnpmid_CAM_zERA = -7_000. * np.log( pmid_CAM_zERA /p_00 )
    
    
    te_ERA_xzCAM = vrg.VertRG( a_x  = te_ERA_xCAM ,
                               zSrc = lnpmid_CAM_zERA ,
                               zDst = lnpmid_CAM ,
                               Gridkey ='tzc' ,   # dstTZHkey
                               kind = 'quadratic' )

    #---------------------
    #  End of special T,PS,PMID interaction



    q_ERA_xzCAM = fullRegrid ( a_ERA = q_ERA ,
                               zSrc = lnpmid_CAM_zERA ,
                               zDst = lnpmid_CAM )
        
    qx = SaturateQ( q=q_ERA_xzCAM , 
                    te=te_ERA_xzCAM ,
                    p=pmid_CAM)

    q_ERA_xzCAM =  qx

    u_ERA_xzCAM = fullRegrid ( a_ERA = u_ERA ,
                               zSrc = lnpmid_CAM_zERA ,
                               zDst = lnpmid_CAM )

    v_ERA_xzCAM = fullRegrid ( a_ERA = v_ERA ,
                               zSrc = lnpmid_CAM_zERA ,
                               zDst = lnpmid_CAM )

    w_ERA_xzCAM = fullRegrid ( a_ERA = w_ERA ,
                               zSrc = lnpmid_CAM_zERA ,
                               zDst = lnpmid_CAM )

    toc = time.perf_counter()
    pTime = f"Subsequent Horz+Vert regridding of ERA5 vars took  {toc - tic:0.4f} seconds"
    print(pTime)


    


    plt.plot(   te_ERA_xCAM[0,:,20000] , lnpmid_CAM_zERA[0,:,20000] )
    plt.plot(   te_ERA_xzCAM[0,:,20000] , lnpmid_CAM[0,:,20000] )

    plt.show()
        
    return te_ERA_xzCAM

def fullRegrid( a_ERA,  zSrc ,  zDst ):
    
    a_ERA_xCAM    = erg.HorzRG( aSrc = a_ERA , 
                                 regrd = regrd , 
                                 srcField=srcf , 
                                 dstField=dstf , 
                                 srcGridkey='tzyx',  # srcTZHkey
                                 dstGridkey='c' )

    a_ERA_xzCAM = vrg.VertRG( a_x  = a_ERA_xCAM ,
                              zSrc = zSrc ,
                              zDst = zDst ,
                              Gridkey ='tzc' ,  # dstTZHkey
                              kind = 'linear' )

    return a_ERA_xzCAM

def SaturateQ ( q , te , p):

    qsat = hum.qsat( p=p, T=te )
    qx=np.minimum( q , qsat )

    return qx

def write_netcdf():
    ntime = np.shape(pdTime_ERA)[0]
    print(ntime)
    
    ilev = (aint_CAM+bint_CAM)* 100_000.
    lev  = (amid_CAM+bmid_CAM)* 100_000.
    if (dstTZHkey == 'tzc' ):
        for itim in np.arange( ntime ):
            dims   = ["ncol","time","lev","ilev"]
            coords = dict( 
                lon  = ( ["ncol"],lon_CAM ),
                lat  = ( ["ncol"],lat_CAM ),
                lev  = ( ["lev"],lev),
                ilev = ( ["ilev"],ilev),
                #time = ( ["time"],pdTime_ERA[itim]),
            )
        
            Wds = xr.Dataset( coords=coords  )
            Wds["time"] = pd.to_datetime( pdTime_ERA[itim] )
            Wds["P_00"] = 100_000.
        
            Dar = xr.DataArray( data=aint_CAM, dims=('ilev',),
                                attrs=dict( description='interface hybrid eta coordinate A-coeff ',units='1',) ,) 
            Wds['hyai'] = Dar

            Dar = xr.DataArray( data=bint_CAM, dims=('ilev',),
                                attrs=dict( description='interface hybrid eta coordinate B-coeff ',units='1',) ,) 
            Wds['hybi'] = Dar
        
            Dar = xr.DataArray( data=area_CAM, dims=('ncol',),
                                attrs=dict( description='Cell area',units='Steradians',) ,) 
            Wds['area'] = Dar

            Dar = xr.DataArray( data=phis_CAM, dims=('ncol',),
                                attrs=dict( description='Surface Geopotential Height',units='m+2 s-2',) ,) 
            Wds['PHIS'] = Dar

            Dar = xr.DataArray( data=phis_ERA_xCAM, dims=('ncol',),
                                attrs=dict( description='ERA Surface Geopotential Height',units='m+2 s-2',) ,) 
            Wds['PHIS_ERA'] = Dar

            Dar = xr.DataArray( data=ps_CAM[itim,:], dims=('ncol',),
                                attrs=dict( description='Surface Pressure',units='Pa',) ,) 
            Wds['PS'] = Dar
    
            Dar = xr.DataArray( data=te_ERA_xzCAM[itim,:,:], dims=('lev','ncol',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['T'] = Dar

            Dar = xr.DataArray( data=q_ERA_xzCAM[itim,:,:], dims=('lev','ncol',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['Q'] = Dar
        
            Dar = xr.DataArray( data=u_ERA_xzCAM[itim,:,:], dims=('lev','ncol',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['U'] = Dar

            Dar = xr.DataArray( data=v_ERA_xzCAM[itim,:,:], dims=('lev','ncol',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['V'] = Dar

            Dar = xr.DataArray( data=w_ERA_xzCAM[itim,:,:], dims=('lev','ncol',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['W'] = Dar

        
            yymmdd = str(pdTime_ERA[itim])[0:10]
            hr=str(pdTime_ERA[itim])[11:13]
            ss = str(int(hr)*3600).zfill(5)
            timetag =  yymmdd+'-'+ss
            foo="/glade/scratch/juliob/ERA5_x_ne30pg3_L58."+ timetag+ ".nc"
            print( foo )
            Wds.to_netcdf( foo )

    code = 1
    return code

def Driver():

    rcode0 = prep()

    rcode1 = proc_1_ERA5_time()
    rcode2 = write_netcdf()

    code = 1
    return code

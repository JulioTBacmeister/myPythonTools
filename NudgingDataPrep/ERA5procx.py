#!/usr/bin/env python
# Import packages 
import sys
import argparse as arg

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

import dask
import dask.array as da

#import VertRegrid as vrg

# "ChatGPI version" --- 
import VertRegridFlex as vrg
print( "Using Flexible parallel/serial VertRegrid ")

# import modules in other directories
sys.path.append('../Utils/')
import GridUtils as GrU
import MakePressures as MkP
import humiditycalcs as hum
import MyConstants as Con

# Reload local packages that are under
# development
importlib.reload( erg )
importlib.reload( vrg )
importlib.reload( SG )
importlib.reload( MkP )
importlib.reload( hum )
importlib.reload( GrU )
importlib.reload( Con )

# Physical Constants
Rgas = Con.Rdry() # 

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

def prep(Dst = 'ne30pg3', DstVgrid='L58',  Src='ERA5'):
    #---------------------------------------------
    # This function sets-up variables and objects 
    # that are need for horizontal and vertical 
    # regridding of ERA reanalyses.
    #---------------------------------------------

    global MyDst, MyDstVgrid, MySrc
    global regrd,srcf,dstf
    global phis_CAM, phis_ERA
    global amid_CAM,bmid_CAM,aint_CAM,bint_CAM
    global lon_CAM,lat_CAM,area_CAM
    global area_ERA5

    # Grid keys for remapping
    global srcHkey,dstHkey,srcTHkey,dstTHkey,srcZHkey,dstZHkey,srcTZHkey,dstTZHkey


    #------- 
    # Begin
    #-------
    tic = time.perf_counter()
    MyDst,MyDstVgrid,MySrc = Dst,DstVgrid,Src

    print( f"In prep Src= {Src} to Dst={Dst} " )
    if (Dst == 'ne30pg3'):
        dstHkey = 'c'
        dst_type='mesh'
        dst_scrip = '/glade/p/cesmdata/cseg/inputdata/share/scripgrids/ne30pg3_scrip_170611.nc'
        dst_TopoFile = '/glade/p/cgd/amp/juliob/bndtopo/latest/ne30pg3_gmted2010_modis_bedmachine_nc3000_Laplace0100_20230105.nc'

    if ((Dst == 'fv0.9x1.25') or (Dst=='fv1x1')):
        dstHkey = 'yx'
        dst_type='grid'
        dst_scrip = '/glade/p/cesmdata/cseg/inputdata/share/scripgrids/fv0.9x1.25_141008.nc'
        dst_TopoFile='/glade/p/cgd/amp/juliob/bndtopo/latest/fv_0.9x1.25_gmted2010_modis_bedmachine_nc3000_Laplace0100_20220708.nc'

    if (Src == 'ERA5'):
        srcHkey = 'yx'
        src_type='grid'
        src_scrip = '/glade/work/juliob/ERA5-proc/ERA5interp/grids/ERA5_640x1280_scrip.nc'
        src_TopoFile = '/glade/work/juliob/ERA5-proc/ERA5interp/phis/ERA5_phis.nc'


    # Look for pre-computed weights file
    # If none, set params to create weights file
    if ( (Src == 'ERA5') and (Dst == 'ne30pg3') ):
        griddir = "/glade/work/juliob/ERA5-proc/ERA5interp/grids/"
        wgts_file_Con = griddir + "ERA5_ne30pg3_Conserv_wgts.nc"
        write_weights = False 
        read_weights = True 
    else:
        wgts_file_Con = "N/A"
        write_weights = False 
        read_weights = False 

    # Get vertical grid from a file.
    # These should be small files but 
    # they aren't always.
    if (DstVgrid == 'L58' ):
        # Read in CAM L58 vertical grid
        dstVgridFile = '/glade/work/juliob/ERA5-proc/CAM-grids/Vertical/GRID_48_taperstart10km_lowtop_BL10_v3p1_beta1p75.nc'

    if (DstVgrid == 'L32' ):
        dstVgridFile = '/glade/p/cesmdata/cseg/inputdata/atm/cam/inic/se/f.e22.FC2010climo.ne30pg3_ne30pg3_mg17.cam6_2_022.002.cam.i.0020-01-01-00000_c200610.nc'

    # Set grid keys for Src ERA5 reanalysis
    srcTHkey  = 't'  + srcHkey
    srcZHkey  = 'z'  + srcHkey
    srcTZHkey = 'tz' + srcHkey

    # Set grid keys for Dst CAM-SE
    dstTHkey  = 't'  + dstHkey
    dstZHkey  = 'z'  + dstHkey
    dstTZHkey = 'tz' + dstHkey
 
    # Get all topo data we will use
    # Read in CAM topography. Also get
    # lon and lat and area for CAM (Dst)
    # grid.
    dsTopo_CAM=xr.open_dataset( dst_TopoFile )
    varsCAM  = list( dsTopo_CAM.variables )
    phis_CAM = dsTopo_CAM['PHIS'].values
    lon_CAM  = dsTopo_CAM['lon'].values
    lat_CAM  = dsTopo_CAM['lat'].values
    if ('area' in varsCAM):
        area_CAM = dsTopo_CAM['area'].values
    else:
        area_CAM = GrU.area2d( lon=lon_CAM, lat=lat_CAM )

    # Read in ERA5 topography
    dsTopo_ERA=xr.open_dataset( src_TopoFile )
    phis_ERA=dsTopo_ERA['Z_GDS4_SFC'].values


    # ----------------------------------------------
    #  Setp regridding machinery
    # ----------------------------------------------
    # Scrip file for ERA5 created by ERA5scrip.ipynb
    dsERAscrip = xr.open_dataset( src_scrip )
    area_ERA5 = np.reshape( dsERAscrip['grid_area'].values , np.shape( phis_ERA ) )
        
    # Make object for Conservative regridding from ERA5
    # grid to CAM target. Scrip files need to be provided even 
    # when a weight file is used
    regrd, srcf, dstf = erg.Regrid( srcScrip = src_scrip , 
                                    srcType  = src_type  ,
                                    dstScrip = dst_scrip ,
                                    dstType  = dst_type  ,
                                    write_weights = write_weights ,
                                    read_weights = read_weights ,
                                    weights_file = wgts_file_Con )
    


    vCAM=xr.open_dataset( dstVgridFile )
    amid_CAM = vCAM['hyam'].values
    bmid_CAM = vCAM['hybm'].values
    aint_CAM = vCAM['hyai'].values
    bint_CAM = vCAM['hybi'].values

    print( f" Src scripfile {src_scrip} " )
    print( f" Dst scripfile {dst_scrip} " )
    print( f" Src topo file {src_TopoFile} " )
    print( f" Dst topo file {dst_TopoFile} " )
    print( f" {DstVgrid} Dst vertical grid from {dstVgridFile} " )


    toc = time.perf_counter()
    pTime = f"Prepping for {Src} to {Dst} proc took  {toc - tic:0.4f} seconds"
    print(pTime)
 
    code = 1
    return code


# Define a function that loads a NetCDF file and returns an xarray dataset
@dask.delayed
def load_file(path):
    ds = xr.open_mfdataset(path ,  data_vars='different', coords='different' )
    return ds


def proc_1_ERA5_time( year=2022, month=11, day=1, hour0=99):

    global pdTime_ERA
    global ps_CAM, phis_ERA_xCAM, phis_CAM
    global te_ERA_xzCAM
    global q_ERA_xzCAM
    global u_ERA_xzCAM
    global v_ERA_xzCAM
    global w_ERA_xzCAM

    monStr=str( year ).zfill(4)+str(month).zfill(2)

    if ( MySrc == 'ERA5'):
        if ( hour0 != 99 ):
            hour1=hour0+5
            ymdh0=str( year ).zfill(4)+str(month).zfill(2)+str(day).zfill(2)+str(hour0).zfill(2)
            ymdh1=str( year ).zfill(4)+str(month).zfill(2)+str(day).zfill(2)+str(hour1).zfill(2)
            ymdh=ymdh0+'_'+ymdh1
        else:
            ymdh0=str( year ).zfill(4)+str(month).zfill(2)+str(day).zfill(2)+"*"
            ymdh=ymdh0

        print( "Time tags for ERA5 files ...")
        print(monStr)
        print(ymdh) 

        era5dir = "/glade/collections/rda/data/ds633.6/e5.oper.an.ml/"
        wrkdir=era5dir+monStr+"/"



    print(dstTZHkey)


    tic_overall = time.perf_counter()
    tic = time.perf_counter()
    
    #Define all file names for later use in dask function
    #-----------------------------------------------------
    spfile= wrkdir + 'e5.oper.an.ml.128_134_sp.regn320sc.'+ymdh+'.nc'
    tfile = wrkdir + 'e5.oper.an.ml.0_5_0_0_0_t.regn320sc.'+ymdh+'.nc'
    qfile = wrkdir + 'e5.oper.an.ml.0_5_0_1_0_q.regn320sc.'+ymdh+'.nc'
    ufile = wrkdir + 'e5.oper.an.ml.0_5_0_2_2_u.regn320uv.'+ymdh+'.nc'
    vfile = wrkdir + 'e5.oper.an.ml.0_5_0_2_3_v.regn320uv.'+ymdh+'.nc'
    wfile = wrkdir + 'e5.oper.an.ml.0_5_0_2_8_w.regn320sc.'+ymdh+'.nc'
    all_ERA_files = [ spfile , tfile, qfile, ufile, vfile, wfile ]

    print( "Using DASK " )
    # Create a list of delayed objects, one for each file
    delayed_datasets = [load_file(path) for path in  all_ERA_files  ]
    # Use dask.compute to load all of the datasets in parallel
    datasets  = dask.compute(*delayed_datasets)
    dsPS_ERA  = datasets[0] 
    dsT_ERA   = datasets[1] 
    dsQ_ERA   = datasets[2] 
    dsU_ERA   = datasets[3] 
    dsV_ERA   = datasets[4] 
    dsW_ERA   = datasets[5] 

    """
    #Serial read
    #--------------
    dsPS_ERA  = xr.open_mfdataset( spfile, data_vars='different', coords='different' )
    dsT_ERA   = xr.open_mfdataset( tfile , data_vars='different', coords='different')
    dsQ_ERA   = xr.open_mfdataset( qfile,  data_vars='different', coords='different' )
    dsU_ERA   = xr.open_mfdataset( ufile , data_vars='different', coords='different')
    dsV_ERA   = xr.open_mfdataset( vfile , data_vars='different', coords='different')
    dsW_ERA   = xr.open_mfdataset( wfile , data_vars='different', coords='different')
    """

    ps_ERA = dsPS_ERA['SP'].values
    te_ERA = dsT_ERA['T'].values
    q_ERA  = dsQ_ERA['Q'].values
    u_ERA  = dsU_ERA['U'].values
    v_ERA  = dsV_ERA['V'].values
    w_ERA  = dsW_ERA['W'].values
    
    #---------------------------
    # Get shape of ERA data
    # nt,nL,nx,ny are set here 
    # and used in rest of code
    #-----------------------------
    nt,nx,ny = np.shape( ps_ERA )
    print( "Ps shape ",nt,nx,ny )
    nt,nL,nx,ny = np.shape( te_ERA )
    print( "Te shape ",nt,nL,nx,ny )

    #    Have a look at ERA global mean surface pressures
    for n in np.arange(nt):
        globPS = np.sum( area_ERA5*ps_ERA[n,:,:] ) / np.sum( area_ERA5 )
        print( "ERA5 Global mean surface pressure=",globPS )
        
    #----------------------------------
    # Create a time array from on of
    # the ERA datasets
    #----------------------------------
    pdTime_ERA = pd.to_datetime( dsT_ERA['time'].values )

    #-----------------------------------------------
    # Get hybrid eta-coordinate coefficients for ERA5
    #-----------------------------------------------
    amid=dsT_ERA['a_model'].values 
    bmid=dsT_ERA['b_model'].values
    aint=dsT_ERA['a_half'].values 
    bint=dsT_ERA['b_half'].values
    print( "shape of a_model ", np.shape( amid ) ) 
    

    toc = time.perf_counter()
    pTime = f"Reading one set of ERA5 vars took  {toc - tic:0.4f} seconds"
    print(pTime)

    # To make this code general there should be a split here,
    # i.e., after reading in all the data to process

    tic = time.perf_counter()
    #-----------------------------------------
    # Horz remap of PHIS - No time dimension
    #-----------------------------------------
    phis_ERA_xCAM = erg.HorzRG( aSrc = phis_ERA , 
                                regrd = regrd , 
                                srcField=srcf , 
                                dstField=dstf , 
                                srcGridkey=srcHkey ,
                                dstGridkey=dstHkey )
    
    
    #-----------------------------------------
    # Calculate difference between phis's
    #-----------------------------------------
    Dphis = phis_ERA_xCAM - phis_CAM
    
        
    #-----------------------------------------
    # Horz remap of PS
    #-----------------------------------------
    ps_ERA_xCAM    = erg.HorzRG( aSrc = ps_ERA , 
                                 regrd = regrd , 
                                 srcField=srcf , 
                                 dstField=dstf , 
                                 srcGridkey=srcTHkey ,
                                 dstGridkey=dstHkey )
    
    #    Have a look at ERA global mean surface pressures
    #    after remapping
    for n in np.arange(nt):
        globPS = np.sum( area_CAM*ps_ERA_xCAM[n,:] ) / np.sum( area_CAM )
        print( "Remapped ERA5 Global mean surface pressure=",globPS )

    
    #-----------------------------------------
    # Horz remap of temperature
    #-----------------------------------------
    te_ERA_xCAM    = erg.HorzRG( aSrc = te_ERA , 
                                 regrd = regrd , 
                                 srcField=srcf , 
                                 dstField=dstf , 
                                 srcGridkey=srcTZHkey ,
                                 dstGridkey= dstHkey )
    

    #-----------------------------------------------------------
    # Find "T_bot" and "P_bot" defined in Williamson & Olson 
    # as occurring at "the first level above 150m ... " 
    # This uses ERA temperature and surface pressure remapped  
    # to CAM horz grid (adapted from Isla's function)
    #-----------------------------------------------------------
    
    tic_150 = time.perf_counter()
    p_00 = 1.0 # a_model and a_half appear to include * P_00 in their defintion so set to 1.0 here
    te_150, pmid_150 = MkP.Pressure_TandP150 ( am=amid ,
                                               bm=bmid ,
                                               ai=aint ,
                                               bi=bint ,
                                               ps=ps_ERA_xCAM ,
                                               te=te_ERA_xCAM , 
                                               p_00=p_00 , 
                                               Gridkey = dstTZHkey )


    toc_150 = time.perf_counter()
    pTime = f"Finding Te150 and P150 took  {toc_150 - tic_150:0.4f} seconds"
    print(pTime)

    
    #-------------------------------------------------------------------
    #                    "CAM surface pressure"
    #-------------------------------------------------------------------
    # We don't actually have ps from CAM, so we make a guess based on
    # ps_ERA_xCAM and te_bot asdescribed above.  In a sense this is a 
    # vertical remapping, so we could call this variable ps_ERA_xzCAM, 
    # but we'll just call it ps_CAM ...
    #-------------------------------------------------------------------
    if (dstTZHkey == 'tzc'):
        nt,nz,ncol = np.shape(te_ERA_xCAM)
        L_bot = nz-1
        te_bot = te_150 # te_ERA_xCAM[:,L_bot,:]
        ps_CAM = np.zeros( (nt , ncol) )
        for i in np.arange( nt ):
            ps_CAM[i,:] = ps_ERA_xCAM[i,:] * np.exp( Dphis / (Rgas*te_bot[i,:]) )
        for n in np.arange(nt):
            globPS = np.sum( area_CAM*ps_CAM[n,:] ) / np.sum( area_CAM )
            print( "Remapped+corrected ERA5 Global mean surface pressure=",globPS )

    if (dstTZHkey == 'tzyx'):
        nt,nz,ny,nx = np.shape(te_ERA_xCAM)
        L_bot = nz-1
        te_bot = te_150 # te_ERA_xCAM[:,L_bot,:]
        ps_CAM = np.zeros( (nt , ny, nx) )
        for i in np.arange( nt ):
            ps_CAM[i,:] = ps_ERA_xCAM[i,:,:] * np.exp( Dphis / (Rgas*te_bot[i,:,:]) )
        for n in np.arange(nt):
            globPS = np.sum( area_CAM*ps_CAM[n,:,:] ) / np.sum( area_CAM )
            print( "Remapped+corrected ERA5 Global mean surface pressure=",globPS )

        
    #-----------------------------------------------------------------------------------------------
    # Now we creat full 4(3)D pressure fields on the ERA vertical grid. These are "CAM pressures" 
    # as they are based on the "CAM surface pressure" derived above.  These are used for 
    # vertical interpolation below.  I'm not sure about this.  The WO2015 document seems to be
    # instructing this ... but you could also use ERA surface press on xCAM.
    #-----------------------------------------------------------------------------------------------
    
    p_00 = 1.0 # a_model and a_half appear to include * P_00 in their defintion so set to 1.0 here
    pmid_CAM_zERA, pint_CAM_zERA, delp_CAM_zERA \
        = MkP.Pressure (am=amid ,
                        bm=bmid ,
                        ai=aint ,
                        bi=bint ,
                        ps=ps_CAM ,
                        p_00=p_00 ,
                        Gridkey = dstTZHkey )



    p_00=100_000.0
    pmid_CAM,pint_CAM,delp_CAM \
        = MkP.Pressure (am=amid_CAM ,
                        bm=bmid_CAM ,
                        ai=aint_CAM ,
                        bi=bint_CAM ,
                        ps=ps_CAM ,
                        p_00=p_00 , 
                        Gridkey = dstTZHkey )



    #-----------------------------------------------------------
    # Log-pressure is preferred for vertical interpolation
    # per Williamson&Olson
    #-----------------------------------------------------------
    lnpint_CAM = -7_000. * np.log( pint_CAM / p_00 )
    lnpmid_CAM = -7_000. * np.log( pmid_CAM / p_00 )
    lnpmid_CAM_zERA = -7_000. * np.log( pmid_CAM_zERA /p_00 )

    toc2 = time.perf_counter()
    pTime = f"Creating 3D P-fields etc., took   {toc2 - tic:0.4f} seconds"
    print(pTime)
    
    print(" going into vertical regrid of T " )
    te_ERA_xzCAM = vrg.VertRG( a_x  = te_ERA_xCAM ,
                               zSrc = lnpmid_CAM_zERA ,
                               zDst = lnpmid_CAM ,
                               Gridkey =dstTZHkey ,
                               kind = 'quadratic' )

    #---------------------
    #  End of special T,PS,PMID interaction



    print(" going into horz+vertical regrid of Q " )
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


    pTime = f"Overall time in this function  {toc - tic_overall:0.4f} seconds"
    print(pTime)
        
    rcode =1 
    return rcode

def fullRegrid( a_ERA,  zSrc ,  zDst , kind='linear' ):
    
    print("Horz RG in fullRegrid " )
    a_ERA_xCAM    = erg.HorzRG( aSrc = a_ERA , 
                                 regrd = regrd , 
                                 srcField=srcf , 
                                 dstField=dstf , 
                                 srcGridkey=srcTZHkey,
                                 dstGridkey=dstHkey )

    print("Vert RG in fullRegrid " )
    a_ERA_xzCAM = vrg.VertRG( a_x  = a_ERA_xCAM ,
                              zSrc = zSrc ,
                              zDst = zDst ,
                              Gridkey=dstTZHkey,
                              kind = kind )

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

            Dar = xr.DataArray( data=amid_CAM, dims=('lev',),
                                attrs=dict( description='mid-level hybrid eta coordinate A-coeff ',units='1',) ,) 
            Wds['hyam'] = Dar

            Dar = xr.DataArray( data=bmid_CAM, dims=('lev',),
                                attrs=dict( description='mid-level hybrid eta coordinate B-coeff ',units='1',) ,) 
            Wds['hybm'] = Dar
        
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
            #foo="/glade/scratch/juliob/ERA5_x_ne30pg3_L58_v6."+ timetag+ ".nc"
            foo="/glade/scratch/juliob/ERA5_x_"+ MyDst + "_"+ MyDstVgrid +"_v6."+ timetag+ ".nc"
            print( foo )
            Wds.to_netcdf( foo )

    if (dstTZHkey == 'tzyx' ):
        for itim in np.arange( ntime ):
            dims   = ["lon","lat","time","lev","ilev"]
            coords = dict( 
                lon  = ( ["lon"],lon_CAM ),
                lat  = ( ["lat"],lat_CAM ),
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

            Dar = xr.DataArray( data=amid_CAM, dims=('lev',),
                                attrs=dict( description='mid-level hybrid eta coordinate A-coeff ',units='1',) ,) 
            Wds['hyam'] = Dar

            Dar = xr.DataArray( data=bmid_CAM, dims=('lev',),
                                attrs=dict( description='mid-level hybrid eta coordinate B-coeff ',units='1',) ,) 
            Wds['hybm'] = Dar
        
            Dar = xr.DataArray( data=area_CAM, dims=('lat','lon',),
                                attrs=dict( description='Cell area',units='Steradians',) ,) 
            Wds['area'] = Dar

            Dar = xr.DataArray( data=phis_CAM, dims=('lat','lon',),
                                attrs=dict( description='Surface Geopotential Height',units='m+2 s-2',) ,) 
            Wds['PHIS'] = Dar

            Dar = xr.DataArray( data=phis_ERA_xCAM, dims=('lat','lon',),
                                attrs=dict( description='ERA Surface Geopotential Height',units='m+2 s-2',) ,) 
            Wds['PHIS_ERA'] = Dar

            Dar = xr.DataArray( data=ps_CAM[itim,:,:], dims=('lat','lon',),
                                attrs=dict( description='Surface Pressure',units='Pa',) ,) 
            Wds['PS'] = Dar
    
            Dar = xr.DataArray( data=te_ERA_xzCAM[itim,:,:,:], dims=('lev','lat','lon',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['T'] = Dar

            Dar = xr.DataArray( data=q_ERA_xzCAM[itim,:,:,:], dims=('lev','lat','lon',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['Q'] = Dar
        
            Dar = xr.DataArray( data=u_ERA_xzCAM[itim,:,:,:], dims=('lev','lat','lon',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['U'] = Dar

            Dar = xr.DataArray( data=v_ERA_xzCAM[itim,:,:,:], dims=('lev','lat','lon',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['V'] = Dar

            Dar = xr.DataArray( data=w_ERA_xzCAM[itim,:,:,:], dims=('lev','lat','lon',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['W'] = Dar

        
            yymmdd = str(pdTime_ERA[itim])[0:10]
            hr=str(pdTime_ERA[itim])[11:13]
            ss = str(int(hr)*3600).zfill(5)
            timetag =  yymmdd+'-'+ss
            #foo="/glade/scratch/juliob/ERA5_x_ne30pg3_L58_v6."+ timetag+ ".nc"
            foo="/glade/scratch/juliob/ERA5_x_"+ MyDst + "_"+ MyDstVgrid +"_v6."+ timetag+ ".nc"
            print( foo )
            Wds.to_netcdf( foo )

    code = 1
    return code

#def Driver():
def main(year,month,day):
    import calendar
    
    tic_total = time.perf_counter()
    days_in_month = calendar.monthrange(year,month)[1]

    print( f"About to process {year:n}-{month:n}-{day:n}")
    rcode0 = prep()
    rcode1 = proc_1_ERA5_time( year=year, month=month, day=day)
    rcode2 = write_netcdf()

    code = 1
    toc_total = time.perf_counter()

    pTime = f"Total processing time was  {toc_total - tic_total:0.4f} seconds"
    print(pTime)

if __name__ == "__main__":
    # argument: indir -> get all nc files in this directory
    # argument: map -> the offlinemap file already prepared
    # argument: outdir -> directory where remapped files should go
    # my_parser = arg.ArgumentParser()
    # my_parser.add_argument("--month", type=int)
    # my_parser.add_argument("--year", type=int)
    # args = my_parser.parse_args()
    
    my_parser = arg.ArgumentParser()
    my_parser.add_argument("--month", type=int, default=1)
    my_parser.add_argument("--year", type=int, default=2010)
    my_parser.add_argument("--day", type=int, default=1)
    args = my_parser.parse_args()
    main(args.year, args.month, args.day )

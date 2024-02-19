##########
# get validation data
##########
workdir_ = '/glade/work/juliob/'
import sys
import os
sys.path.append(workdir_ + 'myPythonTools/GenlTools/')
sys.path.append(workdir_ + 'myPythonTools/Utils/')
sys.path.append(workdir_ + 'PyRegridding/Regridder/')
sys.path.append(workdir_ + 'PyRegridding/Utils/')

# Own local packages
import AveragingUtils as Av
import MakePressures as MkP
import var_A_x_B as vAB
import GridUtils as GrU


import numpy as np
import xarray as xr

import importlib

importlib.reload(vAB)

def data(fld,season,months=-999,**kwargs):

    ##################################
    # Parse out logic etc of kwargs
    ##################################
    if 'mgrid' in kwargs: 
        if kwargs['mgrid'] in [False,None]:
            regrid_x_mgrid = False
        else:
            regrid_x_mgrid = True
    else:
        regrid_x_mgrid = False

    if 'Years' in kwargs:
        yearsA = kwargs['Years']
        ERA5dir = '/glade/campaign/cgd/amp/juliob/ERA5/native_monthly/'
    else:
        yearsA = '1979_2022'
        ERA5dir = '/glade/campaign/cgd/amp/juliob/ERA5/native_monthly_climo/'
    
    if (fld in ('U','V','T','Q','PS','OMEGA') ):
        # get ERA5 pl 
        path_C = ERA5dir + 'ERA5.native.time.' + yearsA +'-*.nc'
        Dc = xr.open_mfdataset( path_C ,data_vars='different', coords='different' )
        print( f" Validation data {fld} has dims {Dc[fld].dims}")
        
        if (Dc.hybm.dims[0] == 'time'):
            hybm = Dc.hybm.values[0,:]
        else:
            hybm = Dc.hybm.values
        if (Dc.hyam.dims[0] == 'time'):
            hyam = Dc.hyam.values[0,:]
        else:
            hyam = Dc.hyam.values
        if (Dc.hybi.dims[0] == 'time'):
            hybi = Dc.hybi.values[0,:]
        else:
            hybi = Dc.hybi.values
        if (Dc.hyai.dims[0] == 'time'):
            hyai = Dc.hyai.values[0,:]
        else:
            hyai = Dc.hyai.values

        lev = 1_000.*hyam[:] +1_000.*hybm[:] 
        lon =Dc.longitude.values
        lat =Dc.latitude.values
        aa = Av.Seasonal( ds=Dc, season=season , fld=fld)

        
        if (regrid_x_mgrid ==True ):
            # Put on model grid for comparisons
            #-----------------------------------
            # mgrid needs the following components
            #   mgrid={ps:ps_m, hyam:hyam_m, hybm:hybm_m, hyai:hyai_m, hybi:hybi_m, hgrid:hgrid_m }
            ps = Av.Seasonal( ds=Dc, season=season , fld='PS')
            mgrid  = kwargs['mgrid']
            obsgrid={'ps':ps, 'hyam':hyam, 'hybm':hybm, 'hyai':hyai, 'hybi':hybi, 'hgrid':'ERA5' }
            dic = to_mgrid( aa, mgrid=mgrid, obsgrid=obsgrid ) 
            if 'zlev' in kwargs:
                dic['lev'] = -7. * np.log( dic['lev'] /1_000. )
        else:
            if 'zlev' in kwargs:
                lev = -7. * np.log( lev /1_000. )
            #----- Pack output into a 'dict'
            dic={'aa':aa,'lev':lev,'lat':lat,'lon':lon,
                 'hyai':hyai,'hyam':hyam,'hybi':hybi,'hybm':hybm,
                 'years':yearsA, 'data_path':path_C,'rcode':0}
        
    else:
        dic = {'aa':-999e10,'rcode':-99}
        
            
    return dic  # aa,lev,lat,lon

def to_mgrid(aa=None,mgrid=None,obsgrid=None):

    print(f" In to_mgrid {mgrid.keys()}")  # Output will be dict_keys(['a', 'b', 'c'])
    
    #-----------------------------------------
    # Make 3D prssure field on Obs grid
    #-----------------------------------------
    pmid_o, pint_o, delp_o \
        = MkP.Pressure (am=obsgrid['hyam'],
                        bm=obsgrid['hybm'] ,
                        ai=obsgrid['hyai'] ,
                        bi=obsgrid['hybi'] ,
                        ps=obsgrid['ps'] ,
                        p_00=100_000. ,
                        Gridkey = 'zyx' )
    print(f" In to_mgrid shap pmid_o {np.shape(pmid_o)}")

    """
    if ( Dst == 'ne30pg3'):
    griddir = "/glade/work/juliob/ERA5-proc/ERA5interp/grids/"
    wgts_file_Con = griddir + "ERA5_ne30pg3_Conserv_wgts.nc"
    write_weights = False 
    read_weights = True 
    """

    pmid_o_Hm , lat_o_Hm , lon_o_Hm \
        = vAB.Hregrid(avar=pmid_o,
                        agrid=obsgrid['hgrid'],
                        akey='zyx',
                        bgrid=mgrid['hgrid'] )

    print(f" In to_mgrid shap pmid_o_Hm {np.shape(pmid_o_Hm)}")
 
    aa_Hm , lat_o_Hm , lon_o_Hm \
        = vAB.Hregrid(avar=aa,
                        agrid=obsgrid['hgrid'],
                        akey='zyx',
                        bgrid=mgrid['hgrid'] )

    print(f" In to_mgrid shap aa_Hm {np.shape(aa_Hm)}")
    print(f" In to_mgrid shap mgrid ps {np.shape(mgrid['ps'])}")

    
    pmid_m, pint_m, delp_m \
        = MkP.Pressure (am=mgrid['hyam'],
                        bm=mgrid['hybm'] ,
                        ai=mgrid['hyai'] ,
                        bi=mgrid['hybi'] ,
                        ps=mgrid['ps'] ,
                        p_00=100_000. ,
                        Gridkey = 'zyx' )

    p_00 = 100_000. # Here we just use the sensible value of p_00
    lnpmid_o_Hm   = -7_000. * np.log( pmid_o_Hm / p_00 )
    lnpmid_m = -7_000. * np.log( pmid_m / p_00 )
    
    aa_m = vAB.VertRG( a_x = aa_Hm ,
                zSrc = lnpmid_o_Hm ,
                zDst = lnpmid_m ,
                Gridkey ='zyx' ,
                kind = 'quadratic' ) #linea


    GridInfo = GrU.gridInfo( grid=mgrid['hgrid'] )
    lat,lon = GrU.latlon( scrip=GridInfo['scrip'],Hkey=GridInfo['Hkey'] )

    lev = 1_000.*mgrid['hyam'][:] +1_000.*mgrid['hybm'][:] 

    oDict = {'aa':aa_m, 'lev':lev, 'lat':lat,'lon':lon }
    
    return oDict





   
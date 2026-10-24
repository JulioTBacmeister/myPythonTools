################################################
# New style 
################################################
import sys
import os
import time

workdir_ = '/glade/work/juliob'
if ( workdir_ not in sys.path ):
    sys.path.append(workdir_)
    print( f" a path to {workdir_} added in {__name__} ")

from PyRegridding.Utils import GridUtils as GrU
from PyRegridding.Drivers import RegridField as RgF

import xarray as xr
import numpy as np

def regrid_X( X, Src='ne30', Dst='fv1x1' ):
    ############################
    # Initialize dict for regridding
    ############################
    RgObs={}

    
    DstInfo = GrU.gridInfo(Dst) #,Vgrid=DstVgrid)
    lat_a,lon_a = GrU.latlon( scrip= DstInfo['scrip'], Hkey=DstInfo['Hkey'] )
    RegridObj = GrU.regrid_object_lib(RgOb=RgObs, src=Src, dst=Dst)
    x_dims = X.dims

    coords = dict( 
        lon  = ( ["lon"], lon_a ),
        lat  = ( ["lat"], lat_a ),
    )
    
    Xo = xr.Dataset( coords=coords  )

    dimas=['time','lev','ilev']
    for dima in dimas:
        if dima in x_dims:
            Xo.coords[dima]=X.coords[dima]

    excl_var = ['lat','lon',]
    for var in X.variables:
        var_dims = X[var].dims
        if ( ( 'ncol' in var_dims) and ( var not in excl_var) ):
            # replace 'ncol' in var_dims with 'lat','lon'
            new_dims = list(var_dims)
            i = new_dims.index('ncol')
            new_dims[i:i+1] = ['lat','lon']
            new_dims = tuple(new_dims)
            print( f"extracting {var} " )
            tic_total_0 = time.perf_counter()
            A_c = X[var].values
            tic_total_1 = time.perf_counter()
            print( f"extracted {var} shape {np.shape(A_c)} took {tic_total_1 - tic_total_0:0.4f} seconds " )
            A_yx = RgF.Horz(xfld_Src=A_c, Src=Src, Dst=Dst , RegridObj_In= RegridObj )
            tic_total_2 = time.perf_counter()
            print( f"regridded {var} shape {np.shape(A_c)} took {tic_total_2 - tic_total_1:0.4f} seconds " )
            Dar = xr.DataArray( data=A_yx, 
                                dims=new_dims, )
            Xo[var]=Dar

    return Xo





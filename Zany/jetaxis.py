#########################
# Imports 
#########################
from PyRegridding.Utils import GridUtils as GrU
from PyRegridding.Drivers import RegridField as RgF


# The usual
from datetime import date
import numpy as np
import xarray as xr

def lat_z(u=None,lat=None,zlev=None,LatLims=[-90,0]  ):
    # expects slices of zonal mean winds nt,nz,ny or nz,ny

    ndims = len( u.shape )
    if ndims == 2:
        nz,ny=np.shape(u)
    elif ndims==3:
        nt,nz,ny = np.shape(u)
        
    lat0,lat1=LatLims[0],LatLims[1]
    L0=np.argmin( abs(lat - lat0 ) )    

    L1=np.argmin( abs(lat - lat1 ) )

    
    if ndims==2:
        jetmax=np.zeros( nz )
        for z in np.arange( nz ):
            y_u_max = np.argmax( u[ z, L0:L1 ] )
            jetmax[ z ] = lat[ L0+ y_u_max ] 
    if ndims==3:
        jetmax=np.zeros( (nt,nz) )
        for t in np.arange( nt ):
            for z in np.arange( nz ):
                y_u_max = np.argmax( u[ t, z, L0:L1 ] )
                jetmax[ t, z ] = lat[ L0+ y_u_max ] 
                            
    
    return jetmax


def make_axes_mean_range( X=None , Src=None, Dst='fv1x1' , season=[6,7,8] ):

    RgObs={}

    AAa = X.U.values

    nt = np.shape( AAa )[0]
    nyears = nt // 12
    
    DstInfo = GrU.gridInfo(Dst) #,Vgrid=DstVgrid)
    lat_a,lon_a = GrU.latlon( scrip= DstInfo['scrip'], Hkey=DstInfo['Hkey'] )

    if Src not in ('fv1x1','fv0.9x1.25'):
        #####################################
        # For unstructured or higher res data
        #####################################
        RegridObj = GrU.regrid_object_lib(RgOb=RgObs, src=Src, dst=Dst)
        AAa_yx=RgF.Horz(xfld_Src=AAa , Src=Src, Dst=Dst , RegridObj_In= RegridObj )
    else:
        AAa_yx = AAa
    
    AAa_zav = np.average( AAa_yx , axis=3)
    iseason = np.array( season ) - 1
    AAa_zav_s = np.stack([
        np.average(AAa_zav[y*12 + iseason , :, :], axis=0)
        for y in np.arange(nyears)
    ])
    
    jaxes=lat_z( u=AAa_zav_s, lat=lat_a )

    mean = np.mean(jaxes, axis=0)
    min_ = np.min(jaxes, axis=0)
    max_ = np.max(jaxes, axis=0)
    std  = np.std(jaxes, axis=0)

    return jaxes,mean,min_,max_,std,nyears





import sys
sys.path.append('../Plotting/')
""" Now you can imprt modules in ../Plotting"""
import xyp_plot as xyp

import numpy as np
import xarray as xr
from scipy.interpolate import LinearNDInterpolator as Li
from scipy.interpolate import NearestNDInterpolator as Ni
from scipy.spatial import Delaunay as Dl


def files():

    tag='ndg04'
    fil1='/glade/scratch/juliob/SAMwrf_'+tag+'_ML_super_v2.nc'
    tag='ndg05'
    fil2='/glade/scratch/juliob/SAMwrf_'+tag+'_ML_super_v2.nc'

    return fil1, fil2

def press(**kwargs):

    if 'file' in kwargs:
        file=kwargs['file']
        d=xr.open_dataset( file )
        hybm=d['hybm']
        hyam=d['hyam']
        ps=d['PS']
    else:
        hybm=kwargs['hybm']
        hyam=kwargs['hyam']
        ps=kwargs['PS']

    dimp=np.shape( ps )
    dimh=np.shape( hyam )
    ndimh=len(dimh)
    """
    If len(dimh) is 2 then dimh[0] is probalbly time
    """
    
    if (ndimh==1):
        nlev=dimh[0]
        ncol=dimp[0]
        p3=np.zeros( (nlev, ncol) )
        for L in np.arange( nlev ):
            p3[L,:] = hyam[L]*100000. + hybm[L]*ps[:]
    elif (ndimh==2):
        ntim=dimh[0]
        nlev=dimh[1]
        ncol=dimp[1]
        p3=np.zeros( (ntim, nlev, ncol) )
        for n in np.arange( ntim ):
            for L in np.arange( nlev ):
                p3[n,L,:] = hyam[n,L]*100000. + hybm[n,L]*ps[n,:]

    return p3

def corr_utn_utgw( fil1, fil2 ):

    d1=xr.open_dataset( fil1 )
    d2=xr.open_dataset( fil2 )

    utn1=d1['UTEND_NDG'].values
    utgw2=d2['UTEND_GWDTOT'].values
    lon=d1['lon'].values
    lat=d1['lat'].values
    lev=d1['lev'].values

    sh=np.asarray(np.shape(utn1))
    ntime=sh[0]
    nlev=sh[1]
    ncol=sh[2]
    corro=np.zeros( (nlev,ncol) )

    for L in np.arange( sh[1]-1,0,-1 ):
        #for L in np.arange( 31,30,-1 ):
        print("Level= ",L)
        for i in np.arange( sh[2] ):
            poo = np.corrcoef( x=utn1[:,L,i], y=utgw2[:,L,i] )
            corro[L,i]=poo[0,1]
            
    corro=np.nan_to_num(corro,nan=0.0)



    return corro   #,lev,lat,lon



def c_o_xy(idata,lon,lat,dx=1.,dy=1.,lonr=[270.,340.],latr=[-60.,20.],verbose=False ):


    # Create grid values first.
    nlon=int( (lonr[1]-lonr[0])/dx )
    nlat=int( (latr[1]-latr[0])/dy )
    xi = np.linspace(270. , 340. , nlon )
    yi = np.linspace(-60., 20., nlat )
    Xi, Yi = np.meshgrid(xi, yi)

    # Calculate Delaunay traingulation
    # Note, need to make a weird array input
    # from lons and lats
    llz=np.c_[lon,lat]
    triang = Dl( llz )

    
    #Determine shape of idata
    dims = np.shape( idata )
    ndims = len(dims)

    if (ndims==2):
        nlev=dims[0]
        odata = np.zeros( [nlev, nlat , nlon ] )
        for L in np.arange(nlev):
            if(verbose==True):
                print(" interpolating Level =",L)
            odata[L,:,:] =LiNi(triang, idata[L,:] ,Xi,Yi)
  
    else:
        odata = np.zeros( [nlat, nlon ] )
        odata[:,:] = LiNi(triang, idata[:],Xi,Yi  )
 

    return odata,Xi,Yi

def LiNi(triang,idatas,Xi,Yi):
    """ 
    This code takes care of problem with linear interpolation 
    using triangulation. When 'hull' extends beyond valid points NaN or fill_value 
    must be used. Nearest neighbor doesn't have this issue. So, we do both. 
    Where linear fails we use the nearest neighbor value.
    """
    
    odataXY_L = Li(triang, idatas ,fill_value=-999999.)
    odataXY_N = Ni(triang, idatas )
    odata_L =odataXY_L( Xi, Yi )
    odata_N =odataXY_N( Xi, Yi )
    odata_O =  np.where( odata_L==-999999. , odata_N   ,   odata_L      )

    return odata_O

def wrt_to_nc(ofile_root,fld,data,**kwargs):

    if 'levlatlon' in kwargs:
        lon=kwargs['lon']
        lat=kwargs['lat']
        lev=kwargs['lev']
        d = { 
            'lon':{'dims':('lon'), 'data':lon },
            'lat':{'dims':('lat'), 'data':lat },
        }

        dou = xr.Dataset.from_dict(d)
        Dar = xr.DataArray( data=data , dims=['lev','lat','lon'] , coords=(lev,lat,lon) , attrs=dict( description=fld,units='N/A',) ,) 
        dou[fld] = Dar

    ofile=ofile_root+'_'+fld+'.nc'
    dou.to_netcdf(ofile)

def irun():

    f1,f2=files()
    corro,lev,lat,lon = corr_utn_utgw(f1,f2)

    d1=xr.open_dataset( f1 )
    ps=d1['PS']
    hyam=d1['hyam']
    hybm=d1['hybm']
    plev=d1['lev']

    print("SHAPES: " )
    print("   PS   -",ps.shape)
    print("   plev -",plev.shape)

    gps=np.average( ps, axis=0 )
    ghya=np.average( hyam , axis=0 )
    ghyb=np.average( hybm , axis=0 )
    gplv=np.average( plev , axis=0 )

    p3=press(PS=gps,hybm=ghyb,hyam=ghya )

    corrx,xlon,xlat= c_o_xy(idata=corro,lon=lon,lat=lat,dx=0.25,dy=0.25)
    p3x,xlon,xlat  = c_o_xy(idata=p3 ,lon=lon,lat=lat ,dx=0.25,dy=0.25)
    print("   xlat -",xlat.shape)

    lon=xlon[0,:]
    lat=xlat[:,0]

    #xyp.pltxp( a=corrx, p=p3x, x=xlon, plev=gplv, j0=30 )

    ofile_root = 'SAMwrf_latlon_timeavg'
    fld='CORR_UTN_UTGW'
    wrt_to_nc(ofile_root=ofile_root ,fld=fld ,data=corrx, levlatlon=True, lon=lon,lat=lat,lev=plev)

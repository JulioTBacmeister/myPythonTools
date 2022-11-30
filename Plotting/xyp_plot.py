import matplotlib.pyplot as plt
import numpy as np


def pltxp( a ,p, x, **kwargs):
    """
    Input: 
    a - 2D (lev,{lat,lon}), 3D (lev,lat,lon) or 4D (time,lev,lat,lon) array to plotted
    p - 2D (lev,{lat,lon}), 3D (lev,lat,lon) or 4D (time,lev,lat,lon) array of pressures
    x is a 1D (lon) array
    """
    if 'lats' in kwargs and 'lat0' in kwargs:
        lats=kwargs['lats']
        lat0=kwargs['lat0']
        j0=np.argmin( abs(lat-lat0) )
    
    if 'j0' in kwargs:
        j0=kwargs['j0']
    
    if 'plev' in kwargs:
        plev=kwargs['plev']

    if 'clevs' in kwargs:
        clevs=kwargs['clevs']
    else:
        clevs=30
   
    if 'ylim' in kwargs:
        ylim=kwargs['ylim']
        plt.ylim( ylim[0],ylim[1] )

    if 'xlim' in kwargs:
        xlim=kwargs['xlim']
        plt.xlim( xlim[0],xlim[1] )

    nlon=x.size
    nlev=plev.size


    #plt.contourf(xx,pp2[0,:,100,:], uu2[0,:,100,:]  ,levels=ulv2*5  ,cmap='RdBu_r')

    if p.ndim==4:
        pp= p[0,:,j0,:]
    if p.ndim==3:
        pp= p[:,j0,:]
    if p.ndim==2:
        pp= p
    if a.ndim==4:
        aa= a[0,:,j0,:]
    if p.ndim==3:
        aa= a[:,j0,:]
    if p.ndim==2:
        aa= a

    zees=np.ones(nlev)
    xx,zz = np.meshgrid(x,zees)


    plt.contourf(xx,pp,aa  ,levels=clevs  ,cmap='RdBu_r')

import matplotlib.pyplot as plt
import numpy as np


def pltxp( a ,p, x, **kwargs):
    """
    Input: 
    a - 2D (lev,{lat,lon}), 3D (lev,lat,lon) or 4D (time,lev,lat,lon) array to plotted
    p - 2D (lev,{lat,lon}), 3D (lev,lat,lon) or 4D (time,lev,lat,lon) array of pressures
    x - 2D (lev,{lat,lon}), 3D (lev,lat,lon) or 4D (time,lev,lat,lon) array of "X" distances

    a and p need to conform
    """
    assert(np.shape(a)==np.shape(p))

    if 'lats' in kwargs and 'lat0' in kwargs:
        lats=kwargs['lats']
        lat0=kwargs['lat0']
        j0=np.argmin( abs(lats-lat0) )
    
    if 'j0' in kwargs:
        j0=kwargs['j0']
    print("j0=",j0)

    if 't0' in kwargs:
        t0=kwargs['t0']
    else:
        t0=0

    if 'plev' in kwargs:
        plev=kwargs['plev']

    if 'clevs' in kwargs:
        clevs=kwargs['clevs']
    else:
        clevs=30
   
    if 'zlim' in kwargs:
        zlim=kwargs['zlim']
        plt.ylim( zlim[0],zlim[1] )

    if 'xlim' in kwargs:
        xlim=kwargs['xlim']
        plt.xlim( xlim[0],xlim[1] )

    #plt.contourf(xx,pp2[0,:,100,:], uu2[0,:,100,:]  ,levels=ulv2*5  ,cmap='RdBu_r')

    if (p.ndim==4):
        pp = p[t0,:,j0,:]
        aa = a[t0,:,j0,:]
        nlev = p.shape[1]

    if (p.ndim==3):
        pp = p[:,j0,:]
        aa = a[:,j0,:]
        nlev = p.shape[0]

    if (p.ndim==2):
        pp = p
        aa = a
        nlev = p.shape[0]

    if (x.ndim==2):
        xv = x[j0,:]

    if (x.ndim==1):
        xv = x


    zees=np.ones(nlev)
    xx,zz = np.meshgrid(xv,zees)

    print('why dont I plot')
    print("xv",xv.shape," nlev",nlev)
    print("pp",pp.shape)
    print("aa",aa.shape)

    plt.contourf(xx,pp,aa  ,levels=clevs  ,cmap='bwr')
    plt.colorbar()
    #plt.show()


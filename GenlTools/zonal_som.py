import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator as Li
from scipy.interpolate import NearestNDInterpolator as Ni
from scipy.spatial import Delaunay as Dl

def gzonal( xc, yc, mask, fc ):

    xc=xc.flatten()
    yc=yc.flatten()
    #ooo=np.transpose( np.asarray([xc,yc]) )  #np.c_[xc,yc]
    ooo=np.c_[xc,yc]

    dlo=Dl(ooo)
    nlat=361
    nlon=720

    X = np.linspace(0., 360. , num=nlon ) # linspace defaults to 50 samples
    Y = np.linspace(-90., 90. ,num=nlat )

    XX,YY=np.meshgrid(X,Y)

    # Interpolote ocean mask to lat lon
    inx = Ni( dlo, mask.flatten() )
    xmask = inx( XX, YY )

    # This makes an array that is 1 over ocean and 0 over land
    # Note, actually only mask=0 is land. mask<0 are inland seas and
    # apparently bhave values in ocean model
    good=np.where( xmask>0, xmask*0.+1., xmask*0.)
    print("lat lon shape", np.shape(good))

    fc = np.nan_to_num(fc , nan=0.)
    inx = Ni( dlo, fc.flatten() )
    fx = inx( XX, YY)
    """ fx is the input array x lat-lon """

    gfx = np.zeros( nlat )

    for j in np.arange( nlat ):
        ggood  = np.sum( good[j,:] )
        if (ggood>0.):
            gfx[j] = np.sum( fx[j,:]*good[j,:] )/ggood
        else:
            gfx[j]=0.

    zafx=fx*0.
    for j in np.arange( nlat ):
        zafx[j,:]= gfx[j]
    zafx=np.where( good==0. , zafx*0.-999999., zafx )


    """ Reverse trasnofrm """
    ooo=np.c_[ XX.flatten(),YY.flatten() ]
    dlbk=Dl(ooo)

    bkinx = Li( dlbk , zafx.flatten() )
    fxc = bkinx( xc, yc )
    
    fxc=np.reshape( fxc , [384,320] )

    """
    
    fig,axs=plt.subplots(2,2)

    ax=axs[0,0]
    im=ax.set_title('Input field on input grid')
    im=ax.tricontourf( xc , yc, fc.flatten() ,   cmap="RdBu_r" )
    plt.colorbar(im,ax=ax)

    ax=axs[1,0]
    im=ax.set_title('Field on lat lon')
    im=ax.contourf( XX , YY, zafx ,   cmap="RdBu_r" , levels=np.linspace(0,100,num=21 ) )
    plt.colorbar(im,ax=ax)
    #fig.colorbar()

    ax=axs[0,1]
    im=ax.set_title('Field on lat lon back to input grid')
    im=ax.tricontourf( xc , yc, fxc.flatten() ,   cmap="RdBu_r" , levels=np.linspace(0,100,num=21 ) )
    plt.colorbar(im,ax=ax)
    #fig.colorbar()

    plt.show()
    """

    return fxc


foo="/glade/p/cesmdata/cseg/inputdata/ocn/docn7/SOM/pop_frc.b.e11.B1850LENS.f09_g16.pi_control.002.20190923.nc"
sm=xr.open_dataset(foo)
zsm=sm

xc=sm['xc'].values
yc=sm['yc'].values

print(xc.shape)

fcq=sm['hblt']
zfcq=fcq

for imo in np.arange(12):
    zfcq[ imo ,:,:] = gzonal( xc=xc, yc=yc, mask=sm['mask'].values , fc=fcq[ imo ,:,:]  )

print(" Out of interpolation" )
print(zfcq.shape)


zsm['hblt']=zfcq
zsm.to_netcdf("goo.nc")

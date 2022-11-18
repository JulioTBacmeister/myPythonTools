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

    X = np.linspace(-1., 361. , num=724 ) # linspace defaults to 50 samples
    Y = np.linspace(-90., 90. ,num=361 )

    XX,YY=np.meshgrid(X,Y)

    inx = Ni( dlo, mask.flatten() )
    xmask = inx( XX, YY )
    good=np.where( xmask>0, xmask*0.+1., xmask*0.)

    fc = np.nan_to_num(fc , nan=0.)
    inx = Ni( dlo, fc.flatten() )

    fx = inx( XX, YY)

    """ Reverse trasnofrm """
    ooo=np.c_[ XX.flatten(),YY.flatten() ]
    dlbk=Dl(ooo)

    bkinx = Li( dlbk , fx.flatten() )
    fxc = bkinx( xc, yc )
    
    fxc=np.reshape( fxc , [384,320] )

    fig,axs=plt.subplots(2,2)

    ax=axs[0,0]
    ax.set_title('Input field on input grid')
    ax.tricontourf( xc , yc, fc.flatten() ,   cmap="RdBu_r" )

    ax=axs[1,0]
    ax.set_title('Field on lat lon')
    ax.contourf( XX , YY, fx ,   cmap="RdBu_r" )

    ax=axs[0,1]
    ax.set_title('Field on lat lon back to input grid')
    ax.tricontourf( xc , yc, fxc.flatten() ,   cmap="RdBu_r" )

    plt.show()


foo="/glade/p/cesmdata/cseg/inputdata/ocn/docn7/SOM/pop_frc.b.e11.B1850LENS.f09_g16.pi_control.002.20190923.nc"
sm=xr.open_dataset(foo)

xc=sm['xc'].values
yc=sm['yc'].values

print(xc.shape)

gzonal( xc=xc, yc=yc, mask=sm['mask'].values , fc=sm['hblt'][0,:,:]  )


print(xc.shape)

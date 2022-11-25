from datetime import date
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator as Li
from scipy.interpolate import NearestNDInterpolator as Ni
from scipy.spatial import Delaunay as Dl

def gzonal( xc, yc, mask, fc, maskRepVal=0. ):

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
    # apparently have values in ocean model
    good=np.where( np.abs(xmask)>0, xmask*0.+1., xmask*0.)
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

    zafx=np.where( good==0. , zafx*0.-9999999., zafx )


    """ Reverse trasnofrm """
    ooo=np.c_[ XX.flatten(),YY.flatten() ]
    dlbk=Dl(ooo)

    bkinx = Ni( dlbk , zafx.flatten() )
    fxc = bkinx( xc, yc )
    
    fxc=np.reshape( fxc , [384,320] )

    # Repair fxc so that no bad values remain over ocean
    # np.where sets questionable points to ZERO. OK for hblt and qdp,
    # but not a good general approach
    #                                         
    #                       ocean True     badval True       over bad    over NOT bad
    fxc_good = np.where( ( abs(mask)>0.)&(abs(fxc)>1000.)    , fxc*0.+maskRepVal   ,   fxc      )
    fxc = fxc_good
    

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

today  = date.today()
yymmdd = today.strftime("%Y%m%d")

#/glade/p/cesmdata/cseg/inputdata/ocn/docn7/SOM/pop_frc.b.e20.B1850.f09_g17.pi_control.all.297.20180523.nc
#/glade/p/cesmdata/cseg/inputdata/ocn/docn7/SOM/mom_frc_b.cesm3_cam058_mom_e.B1850WscMOM.ne30_L58_t061.camdev_cice5.026c_50-99_c20221111.nc

case="cesm2"

if case=="cesm1":
    ifi="/glade/p/cesmdata/cseg/inputdata/ocn/docn7/SOM/pop_frc.b.e11.B1850LENS.f09_g16.pi_control.002.20190923.nc"
    ofi="pop_frc.b.e11.B1850LENS.f09_g16.pi_control.002.ZONAV2."+yymmdd+".nc"
if case=="cesm2":
    ifi="/glade/p/cesmdata/cseg/inputdata/ocn/docn7/SOM/pop_frc.b.e20.B1850.f09_g17.pi_control.all.297.20180523.nc"
    ofi="pop_frc.b.e20.B1850.f09_g17.pi_control.all.297.ZONAV2."+yymmdd+".nc"
if case=="cesm3":
    ifi="/glade/p/cesmdata/cseg/inputdata/ocn/docn7/SOM/mom_frc_b.cesm3_cam058_mom_e.B1850WscMOM.ne30_L58_t061.camdev_cice5.026c_50-99_c20221111.nc"
    ofi="/mom_frc_b.cesm3_cam058_mom_e.B1850WscMOM.ne30_L58_t061.camdev_cice5.026c_50-99.ZONAV2."+yymmdd+".nc"

sm=xr.open_dataset(ifi)
zsm=sm

xc=sm['xc'].values
yc=sm['yc'].values

flds0=list( sm.variables )
flds=[]

for ifld in flds0:
    if sm[ifld].ndim==3:
        flds.append(ifld)

""" 
Unclear whether we should zonally
average all fields. Maybe just
hblt and qdp
"""
flds.remove('dhdx')
flds.remove('dhdy')
flds.remove('S')
flds.remove('T')
flds.remove('U')
flds.remove('V')


print("Fields in input SOM ",flds0)
print("will only zonavg ",flds)
print("will output to ",ofi)

for ifld in flds:
    print("  .., doing ",ifld)
    if (ifld=='hblt'):
        RepVal=20.
    else:
        RepVal=0.

    fcq=sm[ ifld ]
    zfcq=fcq

    for imo in np.arange(12):
        zfcq[ imo ,:,:] = gzonal( xc=xc, yc=yc, mask=sm['mask'].values , fc=fcq[ imo ,:,:] , maskRepVal=RepVal )

    print(" Out of interpolation for ",ifld )
    print(zfcq.shape)

    zsm[ ifld ]=zfcq



zsm.to_netcdf( ofi  )

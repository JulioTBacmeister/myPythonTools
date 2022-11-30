import numpy as np
import xarray as xr

def files():

    tag='ndg04'
    fil1='/glade/scratch/juliob/SAMwrf_'+tag+'_ML_super_v2.nc'
    tag='ndg05'
    fil2='/glade/scratch/juliob/SAMwrf_'+tag+'_ML_super_v2.nc'

    return fil1, fil2

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



    return corro,lev,lat,lon

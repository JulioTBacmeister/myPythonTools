# Load external packages
import numpy as np
import xarray as xr

def cube( xp , dir , year, fld ):

    time=np.arange(12)
    
    for m in range(12):
        mm = str(m+1).zfill(2)
        yy = str(year).zfill(4)
        dtag = yy+'-'+mm
        filn = xp+'.cam.h0.'+dtag+'.nc'
        filn = dir+'/'+xp+'/atm/hist/'+filn
        print(filn)
        a=xr.open_dataset( filn )
        aa=a[fld] #.isel(time=0)
        nx=a.lon.size
        ny=a.lat.size

        if m==0:
            lon=a.lon
            lat=a.lat
            dummy=np.zeros( nx*ny*12)
            dummy=dummy.reshape( nx, ny, 12 )
            cu=xr.DataArray( dummy , coords=[lon,lat,time], dims=['lon','lat','time'] )

        print(aa.size,nx,ny,nx*ny)
        print(aa.dims)
        aa=aa.transpose('lon','lat','time')
        print(aa.dims)
        cu.values[:,:,m] = aa[:,:,0]

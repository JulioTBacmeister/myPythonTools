# Load external packages
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs

# should ahve run these
#   module load conda/latest
#   conda activate npl

poo="/glade/scratch/hannay/archive/b.cesm3_cam058_mom_c.B1850WscMOM.ne30_L58_t061.009/atm/hist/b.cesm3_cam058_mom_c.B1850WscMOM.ne30_L58_t061.009.cam.h0.0008-11.nc"


a=xr.open_dataset( poo )

ts=a['TS'].isel(time=0)

nx=a.lon.size
ny=a.lat.size
ntime=100
time=np.zeros( ntime )

lon=a.lon
lat=a.lat

dummy=np.zeros( nx * ny * ntime )
dummy=dummy.reshape( ny, nx ,100)

#dts = xr.DataArray( dummy , coords=[lon,lat,time] , dims=["lon","lat","time"] )
dts = xr.DataArray( dummy , coords=[lat,lon,time] , dims=["lat","lon","time"] )


lons, lats = np.meshgrid(lon, lat)

dts.values[:,:,0] = ts   #.transpose("lon","lat")

wks, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree()})

#CS = ax.contourf(lons, lats, ts, cmap='nipy_spectral', transform=ccrs.PlateCarree ())
CS = ax.contourf(lons, lats, dts.values[:,:,1]  , cmap='nipy_spectral', transform=ccrs.PlateCarree ())



ax.coastlines()



cb=wks.colorbar(CS , ax=ax, orientation='horizontal',fraction=0.1)

plt.show()
#wks.savefig( "poopoopoo.png" )










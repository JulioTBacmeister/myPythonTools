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


lon=a['lon']
lat=a['lat']

lons, lats = np.meshgrid(lon, lat)


wks, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree()})

CS = ax.contourf(lons, lats, ts, cmap='nipy_spectral', transform=ccrs.PlateCarree ())



ax.coastlines()



cb=wks.colorbar(CS , ax=ax, orientation='horizontal',fraction=0.1)

plt.show()
#wks.savefig( "poopoopoo.png" )










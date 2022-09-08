import one_year as oy
import xarray as xr

year=30
xp='b.cesm3_cam058_mom_c.B1850WscMOM.ne30_L58_t061.009'
fld='TS'
dir='/glade/scratch/hannay/archive/'

oy.cube(fld=fld,dir=dir,year=year,xp=xp)

print("Out in test: TS dims ")
#print(ts.dims)

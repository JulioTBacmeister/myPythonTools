import xarray as xr
import numpy as np

def get_cam_time(tim):
    ymd = get_ymd(tim)
    sec = get_seconds_str(tim)
    return ["-".join([ymd[i].item(),s.item()]) for i,s in enumerate(sec)]

def get_ymd(tim):
    return tim.dt.strftime("%Y-%m-%d")

def get_seconds_str(tim):
    seconds = tim.dt.hour*3600 + tim.dt.minute*60 + tim.dt.second
    seconds_str = [f'{s.item():05d}' for s in seconds]
    return xr.DataArray(seconds_str, dims='time', coords={'time':tim})


f='f.e22r.SAMwrf01.f09.L70.NODEEP_2016_01.cam.h1.2015-01-29-03600.nc'

a=xr.open_dataset(f)

tees = get_cam_time( a.time )
time = a.time

for i in np.arange(24):
    g=a.sel( time= time[i] )
    foo='goo_'+tees[i]+'.nc'
    print(foo)
    g.to_netcdf( foo)




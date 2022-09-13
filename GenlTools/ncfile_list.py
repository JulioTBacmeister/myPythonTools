import numpy as np
import xarray as xr


s="/project/amp/juliob/ERAI/f09_omega/L58/2010/ERAI_fv09_L58.cam2.i.2010-12-31-00000.nc"
ds=xr.open_dataset( s, decode_times=False)

lvar=list( ds.variables)

print("\n")

print("This will give a list of variables in: \n "+s)

print("\n List of variables=  ")
print(lvar)

print("\n Number of variables=  " )
print(len(lvar))




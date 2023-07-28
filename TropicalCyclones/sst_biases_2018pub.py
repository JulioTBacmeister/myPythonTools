import xarray as xr
import numpy as np
import pandas as pd

import importlib

def original():
    
    lnd_file = '/glade/p/cesmdata/cseg/inputdata/atm/cam/topo/fv0.9x1.25-gmted2010_modis-smooth_cam.nc'
    sst_file = '/glade/p/cesmdata/cseg/inputdata/atm/cam/sst/sst_HadOIBl_bc_0.9x1.25_1850_2017_c180507.nc'
    dLnd = xr.open_dataset(lnd_file) 
    lfrac=dLnd.LANDFRAC.values
    
    dS_HadSST = xr.open_dataset(sst_file ) 
    hadsst=dS_HadSST.SST_cpl.values
    dates=dS_HadSST.date
    
    
    b40sst='/glade/collections/cdg/data/cmip5/output1/NSF-DOE-NCAR/CESM1-CAM5/historical/mon/atmos/Amon/r3i1p1/files/ts_20130302/'
    b40sst=b40sst+'ts_Amon_CESM1-CAM5_historical_r3i1p1_185001-200512.nc'
    
    dS_b40 = xr.open_dataset(b40sst)
    
    
    iiii=1584-12
    Ny=20
    print(dS_b40.time_bnds[iiii,:].values[0])
    print(dS_b40.time_bnds[iiii+Ny*12,:].values[0])

    print(dates[iiii].values)

    subsst_Had=hadsst[iiii:iiii+Ny*12,:,:].reshape(Ny,12,192,288) + 273.15
    subsst_b40=dS_b40.ts.values[iiii:iiii+Ny*12,:,:].reshape(Ny,12,192,288)
    print(np.shape(subsst_b40))

    lat=dS_b40.lat.values
    lon=dS_b40.lon.values

    monthly_sst_bias=np.zeros( (12,192,288) )
    for im in np.arange(12):
        monthly_sst_bias[im,:,:] =np.average( subsst_b40[:,im,:,:]-subsst_Had[:,im,:,:], axis=0 )


    return monthly_sst_bias,lat,lon
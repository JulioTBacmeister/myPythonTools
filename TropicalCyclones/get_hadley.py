import sys
# import modules in other directories
sys.path.append('../Utils/')

import numpy as np
import xarray as xr

import pandas as pd

import copy
import importlib
import get_lens1_rcp85 as lens1
import sst_biases_2018pub as sstbias

import trax_util as trx
import ibtracs as IBT


importlib.reload( sstbias )
importlib.reload( lens1 )
importlib.reload( trx )
importlib.reload( IBT )




#===========================================================
# Class to allow things to be accessed via dict.thing syntax
# as well as dict['thing'] syntax. They are equivalent.
#===========================================================
class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


def hadley1x1(limit_dates_for_pub=False):
    
    lnd_file = '/glade/p/cesmdata/cseg/inputdata/atm/cam/topo/fv0.9x1.25-gmted2010_modis-smooth_cam.nc'
    sst_file = '/glade/p/cesmdata/cseg/inputdata/atm/cam/sst/sst_HadOIBl_bc_0.9x1.25_1850_2017_c180507.nc'
    dLnd = xr.open_dataset(lnd_file) 
    lfrac=dLnd.LANDFRAC.values

    dS_HadSST = xr.open_dataset(sst_file ) 
    hadsst=dS_HadSST.SST_cpl.values  + 273.15
    dates=dS_HadSST.date

    lat=dS_HadSST.lat.values
    lon=dS_HadSST.lon.values
    
    hadsst=hadsst.reshape( 2016//12, 12, 192,288)
    dates=dates.values.reshape(2016//12, 12)

    
    if (limit_dates_for_pub == True):
        hadsst = hadsst[-38:-5,:,:,:]
        dates  = dates[-38:-5,:]

    
    return hadsst,lfrac,dates,lat,lon
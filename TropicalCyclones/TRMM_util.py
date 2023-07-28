import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import trax_util as trx
import ibtracs as IBT

import importlib

import time as TimeUtils

import sys
# import modules in other directories
sys.path.append('/glade/work/juliob/PyRegridding/Regridder/')

import esmfRegrid as erg

def shiftlon( lon, precip2d ):
    lon0=np.roll(lon,shift=720,axis=0)
    precip0=np.roll(precip2d ,shift=720, axis=1)
    oo=np.where(lon0<0)
    lon0[oo[0]]=lon0[oo[0]]+360.
    return lon0,precip0 

def timemean( precip, time=-999 ):
    
    IsDataArray = isinstance( precip , xr.DataArray )
    nt,ny,nx = np.shape( precip )
    
    precip_av = np.zeros((ny,nx))
    for t in np.arange( nt ):
        precip_av = precip_av + precip[t,:,:].values

    precip_av = precip_av/nt
    
    return precip_av


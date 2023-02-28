# Import packages 
import sys
sys.path.append('../Plotting/')
""" Now you can imprt modules in ../Plotting"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate as intr


import ESMF as E

import importlib
import glob
import copy
import time

import scripGen as SG
import esmfRegrid as erg


importlib.reload( erg )
importlib.reload(SG)

def VertRG( a_x , zSrc, zDst, Gridkey , fill_value='extrapolate', kind='linear' ):

    # Assumes shapes of a_x, zSrc are the same and conformable with the shape of zDst'
    
    if (Gridkey == 'zc'):
        nzS,ncol = np.shape( zSrc )
        nzD,ncol = np.shape( zDst )
        a_xz = np.zeros( (nzD,ncol) )
        
        tic = time.perf_counter()
        for i in np.arange(ncol):
            fint=intr.interp1d( x = zSrc[:,i], y=a_x[:,i] , 
                                fill_value=fill_values, kind=kind  )
            a_xz[ :, i] = fint(   zDst[:, i ] )
                                   
        toc = time.perf_counter()
        IntrTime = f"Vertical int  {toc - tic:0.4f} seconds"
        print(IntrTime)

    if (Gridkey == 'tzc'):
        nt,nzS,ncol = np.shape( zSrc )
        nt,nzD,ncol = np.shape( zDst )
        a_xz = np.zeros( (nt,nzD,ncol) )
        
        tic = time.perf_counter()
        for n in np.arange(nt):
            print(n,end=',')
            for i in np.arange(ncol):
                fint=intr.interp1d( x = zSrc[n,:,i], y=a_x[n,:,i] , 
                                   fill_value=fill_value, kind=kind  )
                a_xz[n, :, i] = fint(   zDst[n, :, i ] )
                                   
        toc = time.perf_counter()
        IntrTime = f"Vertical int {toc - tic:0.4f} seconds"
        print(IntrTime)
        
        
    return a_xz

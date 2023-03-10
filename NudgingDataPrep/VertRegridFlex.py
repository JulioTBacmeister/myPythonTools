# Import packages 
import sys
sys.path.append('../Plotting/')
""" Now you can imprt modules in ../Plotting"""
import os

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

import multiprocessing 
#from multiprocessing import Pool, cpu_count  #()


def interpolate_column(zSrcT_col, a_xT_col, zDstT_col, fill_value, kind):
    """Interpolate a single column of data."""
    fint = intr.interp1d(x=zSrcT_col, y=a_xT_col, fill_value=fill_value, kind=kind)
    return fint(zDstT_col)


def VertRG( a_x , zSrc, zDst, Gridkey , fill_value='extrapolate', kind='linear' ,npools=1):

    # set number of pools to create
    #nworkers = cpu_count()
    nworkers = len(os.sched_getaffinity(0))
    print('Flexible VertRegrid using reshaped ARRAYS ')
    print('nworkers available ',nworkers)
    if (nworkers>1):
        print('full loop pll-ized according to ChatGPT')
    else:
        print('Using serial code which is faster for 1 worker')
    # setting npools to nworkers might not be considerate of others...
    npools = nworkers

    # Assumes shapes of a_x, zSrc are the same and conformable with the shape of zDst'
    
    if (Gridkey == 'zc'):
        tic = time.perf_counter()
        nzS,ncol = np.shape( zSrc )
        nzD,ncol = np.shape( zDst )
        a_xz = np.zeros( (nzD,ncol) )
        if (nworkers > 1):
            # parallelized loop
            """Interpolate all columns of data in parallel."""
            with multiprocessing.Pool(nworkers) as pool1:
                results = pool1.starmap(
                    interpolate_column,
                    [(zSrc[:, i], a_x[:, i], zDst[:, i], fill_value, kind)
                     for i in range(zSrc.shape[1])] )
                a_xz = np.column_stack(results) 
        else:
            # Serial loop (Faster if nworkers = 1)
            for i in np.arange(ncol):
                fint=intr.interp1d( x = zSrc[:,i], y=a_x[:,i] , 
                                    fill_value=fill_value, kind=kind  )
                a_xz[:, i] = fint(   zDst[:, i ] )
                                  
        toc = time.perf_counter()
        IntrTime = f"Pll'zed Vertical int  {toc - tic:0.4f} seconds"
        print(IntrTime)

    if (Gridkey == 'tzc'):
        # Reshape arrays
        tic = time.perf_counter()
        nt,nzS,ncol = np.shape( zSrc )
        nt,nzD,ncol = np.shape( zDst )
        a_xT        = np.reshape( np.transpose( a_x , (1,0,2) ) , (nzS,nt*ncol) )
        zSrcT       = np.reshape( np.transpose( zSrc , (1,0,2) ) , (nzS,nt*ncol) )
        zDstT       = np.reshape( np.transpose( zDst , (1,0,2) ) , (nzD,nt*ncol) )
        nzS,ntcol = np.shape( zSrcT )
        nzD,ntcol = np.shape( zDstT )
        a_xzT = np.zeros( (nzD,ntcol) )
        
        if (nworkers > 1):
            # parallelized loop
            # Interpolate all columns of data in parallel
            with multiprocessing.Pool(nworkers) as pool1:
                results = pool1.starmap(
                    interpolate_column,
                    [(zSrcT[:, i], a_xT[:, i], zDstT[:, i], fill_value, kind)
                     for i in range(zSrcT.shape[1])] )
                a_xzT = np.column_stack(results)
        else:
            # Serial loop (Faster if nworkers = 1)
            for i in np.arange(ntcol):
                fint=intr.interp1d( x = zSrcT[:,i], y=a_xT[:,i] , 
                                    fill_value=fill_value, kind=kind  )
                a_xzT[:, i] = fint(   zDstT[:, i ] )


        a_xz = np.transpose( np.reshape( a_xzT, (nzD,nt,ncol) ), (1,0,2) )
        toc = time.perf_counter()
        IntrTime = f"Vertical int {toc - tic:0.4f} seconds with {nworkers:n} nworkers"
        print(IntrTime)

    if (Gridkey == 'tzyx'):
        # Reshape arrays
        tic = time.perf_counter()
        nt,nzS,ny,nx = np.shape( zSrc )
        nt,nzD,ny,nx = np.shape( zDst )
        a_xT        = np.reshape( np.transpose( a_x , (1,0,2,3) ) , (nzS,nt*nx*ny) )
        zSrcT       = np.reshape( np.transpose( zSrc , (1,0,2,3) ) , (nzS,nt*nx*ny) )
        zDstT       = np.reshape( np.transpose( zDst , (1,0,2,3) ) , (nzD,nt*nx*ny) )
        nzS,ntcol = np.shape( zSrcT )
        nzD,ntcol = np.shape( zDstT )
        a_xzT = np.zeros( (nzD,ntcol) )
        
        if (nworkers > 1):
            # parallelized loop
            # Interpolate all columns of data in parallel
            with multiprocessing.Pool(nworkers) as pool1:
                results = pool1.starmap(
                    interpolate_column,
                    [(zSrcT[:, i], a_xT[:, i], zDstT[:, i], fill_value, kind)
                     for i in range(zSrcT.shape[1])]  )
                a_xzT = np.column_stack(results) 
        else:
            # Serial loop (Faster if nworkers = 1)
            for i in np.arange(ntcol):
                fint=intr.interp1d( x = zSrcT[:,i], y=a_xT[:,i] , 
                                    fill_value=fill_value, kind=kind  )
                a_xzT[:, i] = fint(   zDstT[:, i ] )

        a_xz = np.transpose( np.reshape( a_xzT, (nzD,nt,ny,nx) ), (1,0,2,3) )
        toc = time.perf_counter()
        IntrTime = f"Pll'zd Vertical int {toc - tic:0.4f} seconds"
        print(IntrTime)
        
        
    return a_xz

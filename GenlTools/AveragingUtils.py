import sys
sys.path.append('../Plotting/')
""" Now you can imprt modules in ../Plotting"""
sys.path.append('../SAMwrf/')
""" Now you can imprt modules in ../Plotting"""


import xyp_plot as xyp
import ana as a

from datetime import date
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# Some useful packages 
import importlib
import copy
import time as ttime
import cftime

def MonthsSeason( season ):
    if (season.lower()=='djf'):
        monthsx=[12,1,2]
    if (season.lower()=='jja'):
        monthsx=[6,7,8]
    if (season.lower()=='mam'):
        monthsx=[3,4,5]
    if (season.lower()=='son'):
        monthsx=[9,10,11]
    return monthsx

def SeasonalZonal( ds, season, **kwargs ):
    
    if 'fld' in kwargs:
        fld = kwargs['fld']
        print( 'Field name ', fld, ' was supplied as argument')
        A = ds[fld]
        #print( ' values method complete')
    elif 'data' in kwargs:
        A = kwargs['data']
    else:
        assert A is not None, "A needs to gien "
    
    if 'dims' in kwargs:
        VarDims = kwargs['dims']
    else:
        VarDims = 'tzyx'
    
    
    
    #A=ds[fld].values
    time_bnds = ds['time_bnds']
    #date=ds['date']
    time=ds['time']
    lats=ds['lat']

    print(np.shape(time))
    print(time_bnds[0].values[0])
    time0=[]
    months=[]
    years=[]
    
    for ixtime in time_bnds:
        time0.append( ixtime.values[0] )
        years.append( ixtime.values[0].year )
        months.append( ixtime.values[0].month )
    
    imos = MonthsSeason( season=season )
    imonths=np.asarray(months)
    season = np.where( ( imonths== imos[0] ) | ( imonths==imos[1] ) | ( imonths==imos[2] ) )
    print( "Indices of months", season[0] )
    nmos=len(season[0]) 
    
    
    tic_mmm_zon = ttime.perf_counter()
    """
    if (VarDims=='levlatlon'):
        A_mmm = np.average( A[ season[0] ,:,:,:].values, axis=0 )
        A_mmm_zon = np.average( A_mmm , axis=2 )
        print(" 3D lev-lat-lon variable " )
    """
    print("New algo")
    if (VarDims=='tzyx'):
        nt,nz,ny,nx=np.shape( A )
        A_zon = np.zeros( (nmos,nz,ny) )
        for n in np.arange( nmos ):
            t = season[0][n]
            print( t ,end=',')
            A_zon[n,:,:] = np.average( A[ t,:,:,:].values, axis=2 )
        A_mmm_zon = np.average( A_zon , axis=0 )
        print(" 3D lev-lat-lon variable " )
    toc_mmm_zon = ttime.perf_counter()
    pTime = f"aVERAGING took  {toc_mmm_zon - tic_mmm_zon:0.4f} seconds"
    print(pTime)
    
    return A_mmm_zon

def Seasonal( ds, season, **kwargs ):
    
    if 'fld' in kwargs:
        fld = kwargs['fld']
        A = ds[fld]
    elif 'data' in kwargs:
        A = kwargs['data']
    else:
        assert A is not None, "A needs to gien "
    
    if 'dims' in kwargs:
        VarDims = kwargs['dims']
    else:
        VarDims = 'tzyx'
    
    
    
    time_bnds = ds['time_bnds']
    #date=ds['date']
    time=ds['time']
    lats=ds['lat']

    print(np.shape(time))
    print(time_bnds[0].values[0])
    time0=[]
    months=[]
    years=[]
    
    for ixtime in time_bnds:
        time0.append( ixtime.values[0] )
        years.append( ixtime.values[0].year )
        months.append( ixtime.values[0].month )
    
    imos = MonthsSeason( season=season )
    imonths=np.asarray(months)
    season = np.where( ( imonths== imos[0] ) | ( imonths==imos[1] ) | ( imonths==imos[2] ) )
    print( "Indices of months", season[0] )
    nmos=len(season[0]) 
    
    if (VarDims=='tzyx'):
        nt,nz,ny,nx = np.shape( A )
        A_mmm = np.zeros( (nz,ny,nx) )
        for n in np.arange( nmos ):
            t = season[0][n]
            print( t ,end=',')
            A_mmm = A_mmm + A[ t,:,:,:].values/nmos
            
        #A_mmm = np.average( A[ season[0] ,:,:,:], axis=0 )
        print(" 3D lev-lat-lon variable " )
    if (VarDims=='tyx'):
        A_mmm = np.average( A[ season[0] ,:,:], axis=0 )
        print(" 2D lev-lat-lon variable " )

        
    return A_mmm

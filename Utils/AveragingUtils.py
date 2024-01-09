workdir_ = '/glade/work/juliob/'
import sys
sys.path.append(workdir_ + 'myPythonTools/GenlTools/')
sys.path.append(workdir_ + 'myPythonTools/Plotting/')
sys.path.append(workdir_ + 'PyRegridding/Regridder/')


#import xyp_plot as xyp
#import ana as a

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
    if (season.lower()=='jan'):
        monthsx=[1]
    if (season.lower()=='feb'):
        monthsx=[2]
    if (season.lower()=='mar'):
        monthsx=[3]
    if (season.lower()=='apr'):
        monthsx=[4]
    if (season.lower()=='may'):
        monthsx=[5]
    if (season.lower()=='jun'):
        monthsx=[6]
    if (season.lower()=='jul'):
        monthsx=[7]
    if (season.lower()=='aug'):
        monthsx=[8]
    if (season.lower()=='sep'):
        monthsx=[9]
    if (season.lower()=='oct'):
        monthsx=[10]
    if (season.lower()=='nov'):
        monthsx=[11]
    if (season.lower()=='dec'):
        monthsx=[12]
    return monthsx

def ListMatch( list1 , list2 ):
    
    theMonth=[]  # This can be any 'element' 
    IndexIn1=[]
    IndexIn2=[]
    # Check each element in list1
    icommon=0
    for idx1, elem1 in enumerate(list1):
        for idx2, elem2 in enumerate(list2):
            if elem1 == elem2:
                theMonth.append(elem1)
                IndexIn1.append(idx1)
                IndexIn2.append(idx2)
                icommon += 1
    
    return IndexIn1, IndexIn2, theMonth 

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
    #time_bnds = ds['time_bnds']
    #date=ds['date']
    time=ds['time']
    lats=ds['lat']

    #print(np.shape(time))
    #print(time_bnds[0].values[0])
    
    if ('time_bnds' in ds):
        print(f" time_bnds present in DataSet " )
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
            
    elif ('time' in ds) :
        if (ds['time'].dtype == 'datetime64[ns]' ):
            print(f" time of type datetime64[ns] present in DataSet " )
            years = ds['time'].values.astype('datetime64[Y]' ).astype(int)
            months = ds['time'].values.astype('datetime64[M]' ).astype(int) %12 +1

        
    imos = MonthsSeason( season=season )
    imonths=np.asarray(months)
    #Iseason = np.where( ( imonths== imos[0] ) | ( imonths==imos[1] ) | ( imonths==imos[2] ) )
    Iseason = ListMatch( list1=imonths, list2=imos )
    print( "Indices of months", Iseason[0] )
    nmos=len(Iseason[0]) 
    
    
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
            t = Iseason[0][n]
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
    
    
    if ('time_bnds' in ds):
        print(f" time_bnds present in DataSet " )
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
            
    elif ('time' in ds) :
        if (ds['time'].dtype == 'datetime64[ns]' ):
            print(f" time of type datetime64[ns] present in DataSet " )
            years = ds['time'].values.astype('datetime64[Y]' ).astype(int)
            months = ds['time'].values.astype('datetime64[M]' ).astype(int) %12 +1

    
    imos = MonthsSeason( season=season )
    imonths=np.asarray(months)
    #Iseason = np.where( ( imonths== imos[0] ) | ( imonths==imos[1] ) | ( imonths==imos[2] ) )
    Iseason = ListMatch( list1=imonths, list2=imos )
    print( "Indices of months", Iseason[0] )
    nmos=len(Iseason[0]) 
    print( f" length of Iseason {nmos} ")
    
    if (VarDims=='tzyx'):
        nt,nz,ny,nx = np.shape( A )
        A_mmm = np.zeros( (nz,ny,nx) )
        for n in np.arange( nmos ):
            t = Iseason[0][n]
            print( t ,end=',')
            A_mmm = A_mmm + A[ t,:,:,:].values/nmos
            
        #A_mmm = np.average( A[ season[0] ,:,:,:], axis=0 )
        print(" 3D lev-lat-lon variable " )
    if (VarDims=='tyx'):
        A_mmm = np.average( A[ Iseason[0] ,:,:], axis=0 )
        print(" 2D lev-lat-lon variable " )

        
    return A_mmm

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

def timemean( precip, time=-999 , Gridkey='tyx', verby=False ):
    
    IsDataArray = isinstance( precip, xr.DataArray )
    if (Gridkey=='tyx'):
        nt,ny,nx = np.shape( precip )

        precip_av = np.zeros((ny,nx))
        for t in np.arange( nt ):
            precip_av = precip_av + precip[t,:,:].values
            if((t%100==0)and(verby==True)):
                print(f" {t} of {nt} "  )

    if (Gridkey=='tc'):
        nt,nc = np.shape( precip )

        precip_av = np.zeros((nc))
        for t in np.arange( nt ):
            precip_av = precip_av + precip[t,:].values
            if((t%100==0)and(verby==True)):
                print(f" {t} of {nt} "  )

    precip_av = precip_av/nt
    
    return precip_av

def monthlymean( precip, time , month, Gridkey='tyx', verby=False ):
    
    IsDataArray = isinstance( precip, xr.DataArray )
    if (Gridkey=='tyx'):
        nt,ny,nx = np.shape( precip )
        ntimes=0
        precip_av = np.zeros((ny,nx))
        for t in np.arange( nt ):
            OurTime=time[t].values.item()
            if(month == OurTime.month): 
                precip_av = precip_av + precip[t,:,:].values
                ntimes=ntimes+1
                if((t%100==0)and(verby==True)):
                    print(f" {t} of {nt} "  )

    if (Gridkey=='tc'):
        nt,nc = np.shape( precip )
        ntimes=0
        precip_av = np.zeros((nc))
        for t in np.arange( nt ):
            OurTime=time[t].values.item()
            if(month == OurTime.month): 
                precip_av = precip_av + precip[t,:].values
                ntimes=ntimes+1
                if((t%100==0)and(verby==True)):
                    print(f" {t} of {nt}: Month={OurTime.month}"  )

    precip_av = precip_av/ntimes
    
    return precip_av

def exceedTime( precip , scale=1.0 ,gridkey='tyx'):
    
    # precip*scale must be in mm d-1
    
    tholds = [50.,100.,300.,500.,1000.]
    nthold = len( tholds )
    
    IsDataArray = isinstance( precip , xr.DataArray )
    if (gridkey=='tyx'):
        nt,ny,nx = np.shape( precip )
        #precipR = np.reshape( precip, (nt,ny*nx) )
        print(f"About to calculate exceed time")
        ExceedTime = np.zeros((nthold, ny*nx))
        for t in np.arange( nt ):
            precval = precip[t,:,:].values.reshape( (ny*nx) ) * scale
            for x in np.arange( nthold ):
                oox = np.where( (precval >= tholds[x]) , 1., 0. )
                ExceedTime[x,:] = ExceedTime[x,:]+oox
            if (t % 200 ==0 ):
                print( t,nt,' in exceedTime ') 

        ExceedTime = np.reshape(ExceedTime,(nthold, ny,nx))
    if (gridkey=='tc'):
        nt,nc = np.shape( precip )
        #precipR = np.reshape( precip, (nt,ny*nx) )
        print(f"About to calculate exceed time")
        ExceedTime = np.zeros((nthold, nc))
        for t in np.arange( nt ):
            precval = precip[t,:].values * scale
            for x in np.arange( nthold ):
                oox = np.where( (precval >= tholds[x]) , 1., 0. )
                ExceedTime[x,:] = ExceedTime[x,:]+oox
            if (t % 200 ==0 ):
                print( t,nt,' in exceedTime ') 

    return ExceedTime


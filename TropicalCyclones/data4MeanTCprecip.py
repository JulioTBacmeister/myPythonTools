import sys
# import modules in other directories
sys.path.append('../Utils')

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import trax_util as trx
import ibtracs as IBT
import TRMM_util as trmm
import Precip_util as prc

import importlib

importlib.reload ( trmm )
importlib.reload ( prc )


##
year0,year1 = 1998,2012
nyrs = 0
prec_grand_av = 0.
for Yr in np.arange( start=year0 , stop=year1+1 ):
    yA = str( Yr ).zfill(4)
    ftrmm = "/glade/p/cgd/amp/juliob/TRMM/3B42_TCs/timeseries/3B42_3hrly_TCPRECT_xTS." + yA + ".nc"
    ds=xr.open_dataset( ftrmm )
    prec_av = trmm.timemean( ds.TCPRECT )
    prec_grand_av = prec_grand_av + prec_av
    print(f"File = {ftrmm} ")
    nyrs=nyrs+1

TRMM_tcprec_grand_av = 24.*prec_grand_av / nyrs
lono=ds.lon.values
lato=ds.lat.values

##
# CAM part 
# input and output Directories
drc_o='/glade/campaign/cgd/amp/juliob/TC-cesm1/precip/'

nyrs = 0
prec_grand_av = 0.
for e in np.arange( start=1,stop=4 ):
    BaseName = trx.pdfname(ens=e,justBaseName=True) 
    print( f"BaseName is  {BaseName} " )
    basename_o = BaseName + 'cam.h4.TCPRECT_xTS.'
    for Yr in np.arange( start=year0 , stop=year1+1 ):
        yA = str( Yr ).zfill(4)
        fname =drc_o+basename_o+yA+'010100Z-'+yA+'123121Z.nc'
        ds=xr.open_dataset( fname )
        prec_av = trmm.timemean( ds.TCPRECT )
        prec_grand_av = prec_grand_av + prec_av
        print(f"File = {fname} ")
        nyrs=nyrs+1

PD_tcprec_grand_av = 1000. * 86400. * prec_grand_av / nyrs
lonc=ds.lon.values
latc=ds.lat.values

##
# CAM part 
BaseName = trx.rcp85fname(sst='sst7',justBaseName=True) 
print( f"BaseName is  {BaseName} " )
basename_o = BaseName + 'cam.h4.TCPRECT_xTS.'
# input and output Directories
drc_o='/glade/campaign/cgd/amp/juliob/TC-cesm1/precip/'


year0,year1 = 2070,2099 #2012
nyrs = 0
prec_grand_av = 0.
for Yr in np.arange( start=year0 , stop=year1+1 ):
    yA = str( Yr ).zfill(4)
    fname =drc_o+basename_o+yA+'010100Z-'+yA+'123121Z.nc'
    ds=xr.open_dataset( fname )
    prec_av = trmm.timemean( ds.TCPRECT )
    prec_grand_av = prec_grand_av + prec_av
    print(f"File = {fname} ")
    nyrs=nyrs+1

SST7_tcprec_grand_av = 1000. * 86400. * prec_grand_av / nyrs


npSaveFile = '/glade/work/juliob/NumPySaves/MeanTC_xTSprecip.npz'
print(f"Writing = {npSaveFile} ")
 
np.savez( npSaveFile , 
         TRMM_tcprec_grand_av=TRMM_tcprec_grand_av, 
         PD_tcprec_grand_av=PD_tcprec_grand_av, 
         SST7_tcprec_grand_av=SST7_tcprec_grand_av, 
         lato=lato, lono=lono, latc=latc, lonc=lonc )

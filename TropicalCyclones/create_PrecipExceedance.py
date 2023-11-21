#!/usr/bin/env python
import argparse as arg


import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import trax_util as trx
import ibtracs as IBT
import TRMM_util as trmm
import Precip_util as prc

import importlib

import dask

import time as TimeUtils

importlib.reload ( prc )


def main(CAM=0,TRMM=0,ens=1,sstA='sst1',year0=-9999,year1=-9999):

    tic_overall = TimeUtils.perf_counter()
    ##
    etag = 'exceedance_xTS'
   
    ###############
    # TRMM part
    ###############
    if (TRMM==1):
        nyrs = 0
        prec_grand_X = 0.

        for Yr in np.arange( start=year0 , stop=year1+1 ):
            yA = str( Yr ).zfill(4)
            if (etag=='exceedance_xTS'):
                ftrmm = "/glade/p/cgd/amp/juliob/TRMM/3B42_TCs/timeseries/3B42_3hrly_TCPRECT_xTS." + yA + ".nc"
            if (etag=='exceedance_o'):
                ftrmm = "/glade/p/cgd/amp/juliob/TRMM/3B42/timeseries/3B42_3hrly." + yA + ".nc"
            ds=xr.open_dataset( ftrmm )
            lon=ds.lon.values
            lat=ds.lat.values
            prec_X = prc.exceedTime( precip=ds.precip , scale=24.)
            prec_grand_X = prec_grand_X + prec_X
            print(f"File = {ftrmm} ")
            npSaveFile = '/glade/work/juliob/NumPySaves/'+etag+'/'+'TRMM-'+etag+'-'+yA+'.npz'
            np.savez( npSaveFile , prec_X=prec_X, lat=lat, lon=lon )
            nyrs=nyrs+1

        toc_overall = TimeUtils.perf_counter()
        pTime = f"Overall time creating prec exceedance   {toc_overall - tic_overall:0.4f} seconds"
        print(pTime)
    


    ###############
    # CAM part
    ##############
    if (CAM==1):
        if (year0>2023):
            BaseName = trx.rcp85fname(sst=sstA,justBaseName=True) 
            print( f"BaseName is  {BaseName} " )
            # Directories and basenames for source data
            if (etag=='exceedance_xTS'):
                drc ='/glade/campaign/cgd/amp/juliob/TC-cesm1/precip/'
                basename = BaseName + 'cam.h4.TCPRECT_xTS.'
            if (etag=='exceedance_o'):
                drc ='/glade/campaign/cgd/ccr/nanr/HPSS/tokeep/TIMESLICE/atm/proc/tseries/hourly3/PRECT/'
                basename = BaseName + 'cam.h4.PRECT.'
        else:
            BaseName = trx.pdfname(ens=ens,justBaseName=True) 
            print( f"BaseName is  {BaseName} " )
            # Directories and basenames for source data
            if (etag=='exceedance_xTS'):
                drc ='/glade/campaign/cgd/amp/juliob/TC-cesm1/precip/'
                basename = BaseName + 'cam.h4.TCPRECT_xTS.'
            if (etag=='exceedance_o'):
                drc ='/glade/campaign/cgd/ccr/nanr/HPSS/tokeep/TIMESLICE/atm/proc/tseries/hourly3/PRECT/'
                basename = BaseName + 'cam.h4.PRECT.'

        if (etag=='exceedance_o'):
            PREC_A = 'PRECT'
            gridkey = 'tc'
        if (etag=='exceedance_xTS'):
            PREC_A = 'TCPRECT'
            gridkey = 'tyx'
           
        nyrs = 0
        prec_grand_X = 0.
        for Yr in np.arange( start=year0 , stop=year1+1 ):
            tic_1yr = TimeUtils.perf_counter()
            yA=str(Yr).zfill(4)
            fname =drc + basename + yA+'010100Z-'+yA+'123121Z.nc'
            ds=xr.open_dataset( fname )
            print(f"Opened file = {fname} ")
            lon=ds.lon.values
            lat=ds.lat.values
            prec_X = prc.exceedTime( precip=ds[PREC_A] , scale=1000.*86400., gridkey=gridkey )
            prec_grand_X = prec_grand_X + prec_X
            print(f"Finished with file = {fname} ")
            npSaveFile = '/glade/work/juliob/NumPySaves/'+etag+'/'+BaseName+'-' +etag+ '-'+yA+'.npz'
            np.savez( npSaveFile , prec_X=prec_X, lat=lat, lon=lon )
            print(f"Wrote file = {npSaveFile} ")
            nyrs=nyrs+1
            toc_1yr = TimeUtils.perf_counter()
            pTime = f"Creating prec exceedance for 1 year took  {toc_1yr - tic_1yr:0.4f} seconds"
            print(pTime)

        toc_overall = TimeUtils.perf_counter()
        pTime = f"Overall time creating prec exceedance   {toc_overall - tic_overall:0.4f} seconds"
        print(pTime)

    
if __name__ == "__main__":

    my_parser = arg.ArgumentParser()
    my_parser.add_argument("--year0",     type=int, default=1979)
    my_parser.add_argument("--year1",     type=int, default=1979)
    my_parser.add_argument("--ens",       type=int, default=1)
    my_parser.add_argument("--sstA",      type=str, default='sst1' )
    my_parser.add_argument("--CAM",       type=int, default=0 )
    my_parser.add_argument("--TRMM",      type=int, default=0 )
    args = my_parser.parse_args()
    main( CAM=args.CAM, TRMM=args.TRMM, sstA=args.sstA, ens=args.ens, year0=args.year0, year1=args.year1 )

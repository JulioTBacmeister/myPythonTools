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


#===========================================================
#               functions
#===========================================================
# ---------------------------
# --- Basic construction of
# --- coincident track and 
# --- TS data
# ---------------------------
def trx_w_ts(sub_basins=False,xNEPacific=False):

    
    if (xNEPacific==True ):
        basins = IBT.basin_outlines(Mediterranean=False,xNEPacific=True)
        XY_xNEPac = np.array(basins[8].xy)
        print(f"You'll be working with xNEPacific sub-basin")
    
    monthly_mean_bias,lat_sst,lon_sst = sstbias.original()
    print( f"Read in monthly mena biases ")
    
    
    
    ts,landf,lat,lon = lens1.FullTS()
    tsC=copy.deepcopy( ts )
    print( f"Read in LENS1 fututre TS ")

    
    tsPD,landfPD,latPD,lonPD = lens1.PresentDayTS()
    #ts=ts[0:30,:,:,:]
    tsPDC=copy.deepcopy( tsPD )
    print( f"Read in LENS1 present day TS ")

    pub_ens_x=[0,8,11,21,24,29,-1]
    SSTlab = ['SST1' , 'SST2' , 'SST3' , 'SST4' , 'SST5' , 'SST6' , 'SST7' ]

    nens,nt,ny,nx = np.shape( ts )
    nYr =  nt//12 
    print(nens,nYr,ny,nx)

    tsC = tsC.reshape( nens ,nYr, 12, ny, nx )
    tsR = np.zeros( (nens ,nYr, 12, ny, nx ) )

    for nen in np.arange( nens ):
        for iy in np.arange( nYr ):
            for im in np.arange( 12 ):
                tsR[nen,iy,im,:,:] = tsC[nen,iy,im,:,:] - monthly_mean_bias[im,:,:]

    nensPD,ntPD,ny,nx = np.shape( tsPD )
    nYrPD =  ntPD//12 
    print(f"Shape of tsPD {nensPD,nYrPD,ny,nx}")

    tsPDC = tsPDC.reshape( nensPD ,nYrPD, 12, ny, nx )
    tsPDR = np.zeros( (nensPD ,nYrPD, 12, ny, nx ) )

    for nen in np.arange( nensPD ):
        for iy in np.arange( nYrPD ):
            for im in np.arange( 12 ):
                tsPDR[nen,iy,im,:,:] = tsPDC[nen,iy,im,:,:] - monthly_mean_bias[im,:,:]

    # This grabs the ensemble members associated with
    # SST1-7 in the paper.
    #tsCpub=tsC[pub_ens_x,:,:,:]
    ts_stack = tsR[pub_ens_x, :-1 ,:,:,:]

    print(f"Shape of TS_stack {np.shape(ts_stack)}")
    nens,nYr,nmo,ny,nx = np.shape(ts_stack)


    alpha_pw = 0.11 # power wind law exponent
    zbot = 60.
    z10m = 10.
    power_wind=(z10m /zbot )**(alpha_pw)
    
    rc1=trx.readtrx( trx.rcp85fname(sst='sst1') )
    rc2=trx.readtrx( trx.rcp85fname(sst='sst2') )
    rc3=trx.readtrx( trx.rcp85fname(sst='sst3') )
    rc4=trx.readtrx( trx.rcp85fname(sst='sst4') )
    rc5=trx.readtrx( trx.rcp85fname(sst='sst5') )
    rc6=trx.readtrx( trx.rcp85fname(sst='sst6') )
    rc7=trx.readtrx( trx.rcp85fname(sst='sst7') )

    BACE_sst1,Yr1,Mon1 =trx.basinACEmonth(wind=rc1.wind,basin=rc1.basin,year=rc1.year,month=rc1.month,power_wind=power_wind)
    BACE_sst2,Yr1,Mon1 =trx.basinACEmonth(wind=rc2.wind,basin=rc2.basin,year=rc2.year,month=rc2.month,power_wind=power_wind)
    BACE_sst3,Yr1,Mon1 =trx.basinACEmonth(wind=rc3.wind,basin=rc3.basin,year=rc3.year,month=rc3.month,power_wind=power_wind)
    BACE_sst4,Yr1,Mon1 =trx.basinACEmonth(wind=rc4.wind,basin=rc4.basin,year=rc4.year,month=rc4.month,power_wind=power_wind)
    BACE_sst5,Yr1,Mon1 =trx.basinACEmonth(wind=rc5.wind,basin=rc5.basin,year=rc5.year,month=rc5.month,power_wind=power_wind)
    BACE_sst6,Yr1,Mon1 =trx.basinACEmonth(wind=rc6.wind,basin=rc6.basin,year=rc6.year,month=rc6.month,power_wind=power_wind)
    BACE_sst7,Yr1,Mon1 =trx.basinACEmonth(wind=rc7.wind,basin=rc7.basin,year=rc7.year,month=rc7.month,power_wind=power_wind)

    
    ###
    BACE_stack = np.stack( (BACE_sst1[:,0:30,:] , BACE_sst2[:,0:30,:] , BACE_sst3[:,0:30,:] , BACE_sst4[:,0:30,:] , 
                          BACE_sst5[:,0:30,:] , BACE_sst6[:,0:30,:] , BACE_sst7[:,0:30,:] ) )
    
    
    
    if (xNEPacific==True ):
        print(f" Setting up xNEPacific sub-basin")
        rc1 = trx.add_sub_basin( rc1, XY_xNEPac[:,0], XY_xNEPac[:,1], sub_basin_number=1 )
        rc2 = trx.add_sub_basin( rc2, XY_xNEPac[:,0], XY_xNEPac[:,1], sub_basin_number=1 )
        rc3 = trx.add_sub_basin( rc3, XY_xNEPac[:,0], XY_xNEPac[:,1], sub_basin_number=1 )
        rc4 = trx.add_sub_basin( rc4, XY_xNEPac[:,0], XY_xNEPac[:,1], sub_basin_number=1 )
        rc5 = trx.add_sub_basin( rc5, XY_xNEPac[:,0], XY_xNEPac[:,1], sub_basin_number=1 )
        rc6 = trx.add_sub_basin( rc6, XY_xNEPac[:,0], XY_xNEPac[:,1], sub_basin_number=1 )
        rc7 = trx.add_sub_basin( rc7, XY_xNEPac[:,0], XY_xNEPac[:,1], sub_basin_number=1 )
        
        sBACE_sst1,Yr1,Mon1 =trx.basinACEmonth(wind=rc1.wind,basin=rc1.sub_basin,year=rc1.year,month=rc1.month,power_wind=power_wind)
        sBACE_sst2,Yr1,Mon1 =trx.basinACEmonth(wind=rc2.wind,basin=rc2.sub_basin,year=rc2.year,month=rc2.month,power_wind=power_wind)
        sBACE_sst3,Yr1,Mon1 =trx.basinACEmonth(wind=rc3.wind,basin=rc3.sub_basin,year=rc3.year,month=rc3.month,power_wind=power_wind)
        sBACE_sst4,Yr1,Mon1 =trx.basinACEmonth(wind=rc4.wind,basin=rc4.sub_basin,year=rc4.year,month=rc4.month,power_wind=power_wind)
        sBACE_sst5,Yr1,Mon1 =trx.basinACEmonth(wind=rc5.wind,basin=rc5.sub_basin,year=rc5.year,month=rc5.month,power_wind=power_wind)
        sBACE_sst6,Yr1,Mon1 =trx.basinACEmonth(wind=rc6.wind,basin=rc6.sub_basin,year=rc6.year,month=rc6.month,power_wind=power_wind)
        sBACE_sst7,Yr1,Mon1 =trx.basinACEmonth(wind=rc7.wind,basin=rc7.sub_basin,year=rc7.year,month=rc7.month,power_wind=power_wind)
        
        sBACE_stack = np.stack( (sBACE_sst1[:,0:30,:] , sBACE_sst2[:,0:30,:] , sBACE_sst3[:,0:30,:] , sBACE_sst4[:,0:30,:] , 
                          sBACE_sst5[:,0:30,:] , sBACE_sst6[:,0:30,:] , sBACE_sst7[:,0:30,:] ) )
       
        print(f" finsihed up sBACE for xNEPacific sub-basin")
    
    
    lnd_file = '/glade/p/cesmdata/cseg/inputdata/atm/cam/topo/fv0.9x1.25-gmted2010_modis-smooth_cam.nc'
    sst_file = '/glade/p/cesmdata/cseg/inputdata/atm/cam/sst/sst_HadOIBl_bc_0.9x1.25_1850_2017_c180507.nc'
    dLnd = xr.open_dataset(lnd_file) 
    lfrac=dLnd.LANDFRAC.values

    dS_HadSST = xr.open_dataset(sst_file ) 
    hadsst=dS_HadSST.SST_cpl.values  + 273.15
    dates=dS_HadSST.date



    pd1=trx.readtrx( trx.pdfname(ens=1) )
    pd2=trx.readtrx( trx.pdfname(ens=2) )
    pd3=trx.readtrx( trx.pdfname(ens=3) )

    BACE_pd1,Yr1,Mon1 =trx.basinACEmonth(wind=pd1.wind,basin=pd1.basin,year=pd1.year,month=pd1.month,power_wind=power_wind)
    BACE_pd2,Yr2,Mon2 =trx.basinACEmonth(wind=pd2.wind,basin=pd2.basin,year=pd2.year,month=pd2.month,power_wind=power_wind)
    BACE_pd3,Yr2,Mon2 =trx.basinACEmonth(wind=pd3.wind,basin=pd3.basin,year=pd3.year,month=pd3.month,power_wind=power_wind)
    
    print(2016./12.)
    hadsstR=hadsst.reshape( 2016//12, 12, 192,288)
    datesR=dates.values.reshape(2016//12, 12)
    print(datesR[-38,:])
    print(datesR[-6,:])

    had_stack = np.stack ( ( hadsstR[-38:-5,:,:,:], hadsstR[-38:-5,:,:,:], hadsstR[-38:-5,:,:,:] ) )
    BACE_pd_stack = np.stack( (BACE_pd1[:,1:,:] , BACE_pd2[:,1:,:] , BACE_pd3[:,1:,:] ))
    print(np.shape(BACE_pd_stack))
    print(np.shape(had_stack))
    
    
    if (xNEPacific != True ):
        return BACE_stack,ts_stack,BACE_pd_stack,had_stack,landf,lat,lon
    else:
        return BACE_stack,ts_stack,BACE_pd_stack,had_stack,landf,lat,lon , sBACE_stack
 

    ###################################################################################

def sst_composite_1(basin,BACE_stack,ts_stack,landf,months,dont_scale=False):
    
    b=basin
    month_idxs = np.array( months )-1

    nens,nYr,nMo,ny,nx = np.shape(ts_stack[:,:,month_idxs,:,:])
    ts_stackR = np.reshape( ts_stack[:,:,month_idxs,:,:]   , (nens*nYr*nMo, ny,nx) )
    tmp = BACE_stack[:,:,:,month_idxs]
    BACE_stackR = np.reshape( tmp[:,b,:,:]   , (nens*nYr*nMo ) )
    
    ACEo = np.average(BACE_stackR)
    
    zoo = np.where( BACE_stackR > ACEo)
    ts0 = np.average( ts_stackR , axis=0 )
    ts1 = np.average( ts_stackR[zoo[0],:,:] , axis=0 )
    dts_over = np.where( landf<.001, ts1 - ts0, np.nan ) 

    zoo = np.where( BACE_stackR < ACEo)
    ts0 = np.average( ts_stackR , axis=0 )
    ts1 = np.average( ts_stackR[zoo[0],:,:] , axis=0 )
    dts_under = np.where( landf<.001, ts1 - ts0, np.nan ) 
    
    if (dont_scale==False):
        # Report back 'monthly mean ACE' scaled to 12 months
        ACEo = ACEo * (nMo/12.)
    
    return ACEo, dts_over, dts_under

#!/usr/bin/env python
workdir_ = '/glade/work/juliob/'
import sys
#######################################
# Leave this for now. But it should change to better
# method as here:
# import os
# module_a_dir = os.path.dirname(os.path.abspath(__file__))
# utils_path = os.path.join(module_a_dir, '..', 'Utils')
# sys.path.append(utils_path)
# print( f" a path added in {__name__} {utils_path} ")
########################################
sys.path.append(workdir_ + 'myPythonTools/GenlTools/')
sys.path.append(workdir_ + 'myPythonTools/Utils/')
#sys.path.append(workdir_ + 'PyRegridding/Regridder/')
sys.path.append(workdir_ + 'PyRegridding/Utils/')

# Own local packages
import AveragingUtils as Av
#import VertRegridFlexLL as Vrg  # This is toxic for some reason
import PlotUtil as Pu
import validation_data as Val
import var_A_x_B as vAB
import MakePressures as MkP
import GridUtils as GrU
import utils as uti

# The usual
from datetime import date
import numpy as np
import xarray as xr



if __name__ == "__main__":
    ################################################
    
    exp='c64_gwR2_ne30pg3_FMTHIST_topfix_rdgres_x01'
    
    for y in [1985,1986]:
        for m in np.arange(start=1,stop=13):
            yyyy=str(y)
            mm=str(m).zfill(2)
            ymdPat =f'{yyyy}-{mm}-*'
            print(ymdPat)            
            A = uti.MakeDict4Exp( exp=exp  , user='juliob', subd='hist' , 
                                     hsPat='cam.h2a' , ymdPat=ymdPat ,verbose=True, open_dataset=True )            
            zAX = A.X.mean(dim='lon')
            
            
            f=f'/glade/derecho/scratch/juliob/archive/c64_gwR2_ne30pg3_FMTHIST_topfix_rdgres_x01/atm/hist/c64_gwR2_ne30pg3_FMTHIST_topfix_rdgres_x01.cam.h2aZonal.{yyyy}-{mm}.nc'
            zAX.to_netcdf( f )


    


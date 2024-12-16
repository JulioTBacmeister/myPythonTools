#!/usr/bin/env python
################################################
# New style 
################################################
import sys
import os

workdir_ = '/glade/work/juliob'
if ( workdir_ not in sys.path ):
    sys.path.append(workdir_)
    print( f" a path to {workdir_} added in {__name__} {utils_path} ")

# Own local packages
from myPythonTools.Utils import AveragingUtils as Av
from myPythonTools.Utils import validation_data as Val
from myPythonTools.Utils import PlotUtil as Pu
from myPythonTools.Utils import utils as uti

from PyRegridding.Utils import GridUtils as GrU
from PyRegridding.Utils import MakePressures as MkP



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


    


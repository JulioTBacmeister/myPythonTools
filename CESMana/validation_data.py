##########
# get validation data
##########
workdir_ = '/glade/work/juliob/'
import sys
sys.path.append(workdir_ + 'myPythonTools/GenlTools/')
sys.path.append(workdir_ + 'myPythonTools/Utils/')
sys.path.append(workdir_ + 'PyRegridding/Regridder/')

# Own local packages
import AveragingUtils as Av


import numpy as np
import xarray as xr

def data(fld,season,months=-999,**kwargs):
    
    if 'Years' in kwargs:
        yearsA = kwargs['Years']
    else:
        yearsA = '*'
    
    if (fld=='U'):
        # get ERA5 pl 
        path_C = '/glade/campaign/collections/rda/data/ds633.1/e5.moda.an.pl/' + yearsA + '/e5.moda.an.pl.128_131_u.ll025uv.*.nc'
        Dc = xr.open_mfdataset( path_C ,data_vars='different', coords='different' )
        ###Dc.time.values.astype('datetime64[Y]' ) 
        lev =Dc.level.values
        lon =Dc.longitude.values
        lat =Dc.latitude.values
        
        aa = Av.Seasonal( ds=Dc, season=season , fld='U')
        if 'zlev' in kwargs:
            lev = -7. * np.log( lev /1_000. )
    else:
        aa,lev,lat,lon = -999e10, -999, -999, -999
            
    return aa,lev,lat,lon
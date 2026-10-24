##### 
# netcdf tools

import numpy as np
import xarray as xr
import os

def wrt_objprec(year=None,exp=None,lat=None,lon=None,objName=None,objmsk=None,objprec=None,prec=None,lwrite=True):
    
    archive_base = f'/glade/derecho/scratch/juliob/archive/moaap_ana'
    ##################################################################
    # Set up dataset for regridded data
    ##################################################################
    Bdiro=f"{archive_base}/{exp}"  #{case}.cam.h1i.{date_tag}.nc"
    os.makedirs( Bdiro , exist_ok=True )
    # f"{year:04}-{month:02}-{d:02}-{3600*h:05}" 
    #fout = f"{Bdiro}/{case}.{sTag}.{hsPat}.{year:04}.nc"
    fout = f"{Bdiro}/{exp}-{objName}.{year:04}.nc"
    
    print( f"Will write {fout}", flush=True  )

    nday,ny,nx = np.shape( prec )

    time_values = np.zeros( nday ).astype(int)
    for d in np.arange(nday):
        time_values[d] = year*1000 + (d+1)

    time_values = np.zeros( nday ).astype(int)
    day_of_year = np.zeros( nday ).astype(int)
    
    for d in np.arange(nday):
        time_values[d] = year*1000 + (d+1)
        day_of_year[d] = (d+1)
    
    if nday==365:
        month_len=[31,28,31,30,31,30,31,31,30,31,30,31]
    elif nday==366:
        month_len=[31,29,31,30,31,30,31,31,30,31,30,31]
    
    cumday=np.cumsum( month_len )
    cumdayx=np.concatenate( ([0],cumday) )
    year_val=np.zeros( nday ).astype(int)
    year_val[:] = year
    month=np.zeros( nday ).astype(int)
    for d in np.arange(nday):
        month[d] = int(np.searchsorted(cumday, day_of_year[d]  ) + 1)
    day=np.zeros( nday ).astype(int)
    for d in np.arange(nday):
        day[d] = day_of_year[d] - cumdayx[ month[d]-1 ] 

    dims = ["time","yc","xc"]

    coords = dict( 
        time = ( ["time"],  time_values.astype(int) ), #pd.to_datetime( pdTime_ERA[itim] ) ),
    )



    Xo = xr.Dataset( coords=coords  )
    
    Dar = xr.DataArray( data=year_val, dims=('time',),
                        attrs=dict( description='month',units='1',) ,) 
    Xo['year'] = Dar

    Dar = xr.DataArray( data=month, dims=('time',),
                        attrs=dict( description='year',units='1',) ,) 
    Xo['month'] = Dar

    Dar = xr.DataArray( data=day, dims=('time',),
                        attrs=dict( description='day',units='1',) ,) 
    Xo['day'] = Dar

    Dar = xr.DataArray( data=lat, dims=('yc','xc',),
                        attrs=dict( description='latitude',units='deg',) ,) 
    Xo['lat'] = Dar

    Dar = xr.DataArray( data=lon, dims=('yc','xc',),
                        attrs=dict( description='longitude',units='deg',) ,) 
    Xo['lon'] = Dar

    Dar = xr.DataArray( data=objmsk, dims=('time','yc','xc',),
                        attrs=dict( description='Obj-frequency per day',units='d-1',) ,) 
    Xo[f"{objName}-freq"] = Dar

    Dar = xr.DataArray( data=objprec, dims=('time','yc','xc',),
                        attrs=dict( description='Obj precipitation',units='mm h-1',) ,) 
    Xo[f"{objName}-PR"] = Dar

    Dar = xr.DataArray( data=prec, dims=('time','yc','xc',),
                        attrs=dict( description='Precipitation',units='mm h-1',) ,) 
    Xo[f"PR"] = Dar

    if (lwrite==True):
        Xo.to_netcdf( fout )

    
    return Xo



def rd_objprec(year=None,exp=None,objName=None,months=None,roll=False):
    
    archive_base = f'/glade/derecho/scratch/juliob/archive/moaap_ana'
    ##################################################################
    # Set up dataset for regridded data
    ##################################################################
    Bdiro=f"{archive_base}/{exp}"  #{case}.cam.h1i.{date_tag}.nc"
    fout = f"{Bdiro}/{exp}-{objName}.{year:04}.nc"
    
    print( f"Will read {fout}", flush=True  )

    X=xr.open_dataset( fout )

    if (months is not None):
        #X = X.where(X["month"].isin(months), drop=True)
        mask = X["month"].isin(months)   # dims: (time,)
        X = X.isel(time=mask)

    if (roll==True):
        nroll=X.sizes['xc']//2
        X = X.roll(xc=nroll, roll_coords=True)
        # convert lon to 0..360
        #X = X.assign_coords(xc=(X.xc % 360))
        X.lon.values = (X.lon.values % 360)
        # (optional but usually good) ensure lon is monotonic increasing
        #X = X.sortby("lon")

    return X






    



    
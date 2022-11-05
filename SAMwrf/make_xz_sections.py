#!/usr/bin/env python

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import argparse as arg



def i_c_o_xy(fili,flds,filo='test.nc',dx=1.,dy=1.,lonr=[270.,340.],latr=[-60.,20.] ):

    a=xr.open_mfdataset( fili )
    nlev=a.dims['lev']
    nilev=a.dims['ilev']
    ncol=a.dims['ncol']
    lon=a['lon']
    lat=a['lat']
    vlist=list(a.variables)

    if ('time' in a.dims ):
        ntime=a.dims['time']
        vtime=a['time'].values
        NeedToPadTime=False
    else:
        ntime=1
        vtime=np.arange(ntime)
        NeedToPadTime=True

    # Create grid values first.
    nlon=int( (lonr[1]-lonr[0])/dx )
    nlat=int( (latr[1]-latr[0])/dy )
    xi = np.linspace(270. , 340. , nlon )
    yi = np.linspace(-60., 20., nlat )
    Xi, Yi = np.meshgrid(xi, yi)
    triang = tri.Triangulation( lon, lat)

    """
    Set up 'dict' for new xarray dataset.  This can then
    turned into basis for dataset, which can then be 
    added to.  'dict' sets up dimensions and puts 
    lat and lon vectors
    """
    d = { 
        'lon':{'dims':('nlon'), 'data':xi },
        'lat':{'dims':('nlat'), 'data':yi },
        'time':{'dims':('ntime'), 'data':vtime },
        'lev':{'dims':('lev'), 'data':a['lev'].values },
        'ilev':{'dims':('ilev'), 'data':a['ilev'].values },
        'nbnd':{'dims':('nbnd'), 'data':[0,1] },
        }
    b = xr.Dataset.from_dict(d)
    b['lat']=xr.DataArray( data=yi ,dims=['nlat'],coords=[b.lat] )
    b['lon']=xr.DataArray( data=xi ,dims=['nlon'],coords=[b.lon] )


    if( 'hyai' in vlist):
        b['hyai']= xr.DataArray( data=a['hyai'].values , dims=['ilev'] , coords=[b.ilev] )
        b['hyam']= xr.DataArray( data=a['hyam'].values , dims=['lev'] , coords=[b.lev] )
        b['hybi']= xr.DataArray( data=a['hybi'].values , dims=['ilev'] , coords=[b.ilev] )
        b['hybm']= xr.DataArray( data=a['hybm'].values , dims=['lev'] , coords=[b.lev] )


    for ifld in flds:
        idata = a[ifld].values
        idims = a[ifld].dims
        if(NeedToPadTime==True):
            shap=(1,)+np.shape(idata)
            idata=np.reshape(idata,shap)

        if ('ilev' in idims):
            odata = np.zeros( [ntime, nilev, nlat , nlon ] )
            for L in np.arange(nilev):
                print(" interpolating Level =",L)
                for n in np.arange(ntime):
                    odataXY = tri.LinearTriInterpolator(triang, idata[n,L,:]  )
                    odata[n, L,:,:] =odataXY( Xi, Yi )
                    oDarra = xr.DataArray( data=odata ,dims=['ntime','ilev','nlat','nlon'],coords=[b.time, b.ilev, b.lat, b.lon] )

        elif ('lev' in idims):
            odata = np.zeros( [ntime, nlev, nlat, nlon ] )
            for L in np.arange(nlev):
                print(" interpolating Level =",L)
                for n in np.arange(ntime):
                    odataXY = tri.LinearTriInterpolator(triang, idata[n,L,:]  )
                    odata[n, L,:,:] =odataXY( Xi, Yi )
                    oDarra = xr.DataArray( data=odata ,dims=['ntime','lev','nlat','nlon'],coords=[b.time, b.lev, b.lat, b.lon] )

        else:
            odata = np.zeros( [ntime, nlat, nlon ] )
            for n in np.arange(ntime):
                odataXY = tri.LinearTriInterpolator(triang, idata[n,:]  )
                odata[n,:,:] =odataXY( Xi, Yi )
                oDarra = xr.DataArray( data=odata ,dims=['ntime','nlat','nlon'],coords=[b.time, b.lat, b.lon] )

        print(oDarra.dims)
        print(np.shape(odata) )

        b[ifld]=oDarra

    b.to_netcdf( filo )

def stuff():

    tag='ndg04'
    fil1='/glade/scratch/juliob/SAMwrf_'+tag+'_ML_super.nc'
    tag='ndg05'
    fil2='/glade/scratch/juliob/SAMwrf_'+tag+'_ML_super.nc'
    
    
    d1=xr.open_dataset( fil1 )
    d2=xr.open_dataset( fil2 )
    
    utn1=d1['UTEND_NDG'].values
    utgw2=d2['UTEND_GWDTOT'].values
    lon=d1['lon'].values
    lat=d1['lat'].values
    
    mos= np.array( [  [2010,6,1,0]
                  ] )
    tag='ndg05'
    xp='c6_3_59.ne30pg3_L32_SAMwrf.'+tag
    dir='/glade/scratch/juliob/archive/'+xp+'/atm/hist/'
    i=0
    fili=dir+xp+'.cam.h1.'+str(mos[i,0]).zfill(4) + '-' + str(mos[i,1]).zfill(2) + '-' + str(mos[i,2]).zfill(2) + '-' +str(mos[i,3]).zfill(5) + '.nc'
    topo=xr.open_dataset( a.attrs['topography_file'] )


def main(month=6,year=2010):

    fo='f.e22r.SAMwrf01.ne30_latlon025.L32.NODEEP_2010_01.cam.h1.2010-06-15-00000.nc'
    fi='/glade/p/cesm/amwg_dev/juliob/SAMwrf/ne30x16/o/f.e22r.SAMwrf01.ne30.L32.NODEEP_2010_01.cam.h1.2010-06-15-00000.nc' 
    i_c_o_xy(fili=fi,filo=fo,flds=['PS','U','V','T'],dx=0.25,dy=0.25)

    fi='/glade/scratch/juliob/archive/c6_3_59.ne30pg3_L32_SAMwrf.ndg05/atm/hist/c6_3_59.ne30pg3_L32_SAMwrf.ndg05.cam.h1.2010-06-15-00000.nc'
    fo='c6_3_59.ne30pg3_L32_latlon025_SAMwrf.ndg05.cam.h1.2010-06-15-00000.nc'
    i_c_o_xy(fili=fi,filo=fo,flds=['PS','U','V','T'],dx=0.25,dy=0.25)


if __name__ == "__main__":
    # argument: indir -> get all nc files in this directory
    # argument: map -> the offlinemap file already prepared
    # argument: outdir -> directory where remapped files should go
    my_parser = arg.ArgumentParser()
    my_parser.add_argument("--month", type=int)
    my_parser.add_argument("--year", type=int)
    args = my_parser.parse_args()
    main(args.year, args.month )

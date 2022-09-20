class hfdata:

    def __init__(self, xp, dir, moniker='cam.h0'):
        self.case    = xp
        self.archdir = dir


    def filename(self, year, month, day, second, moniker='cam.h0'):

        xp=self.case
        dir=self.archdir
        self.moniker=moniker
        self.year=year
        self.month=month
        self.day=day
        self.second=second

        #construct file name
        yy=str( year   ).zfill(4)
        mm=str( month  ).zfill(2)
        dd=str( day    ).zfill(2)
        ss=str( second ).zfill(5)


        file = dir + '/' + xp + '.' + moniker + '.' + yy + '-' + mm + '-' + dd + '-' + ss +'.nc'

        monfile = dir + '/' + xp + '.' + moniker + '.' + yy + '-' + mm +'.nc'
        yearfile = dir + '/' + xp + '.' + moniker + '.' + yy +'.nc'

        self.file=file
        self.monthly_file=monfile
        self.yearly_file=yearfile

        return file

    def point_output(self, ifile , ofile, ilon=180., ilat=0. ):
        import numpy as np
        import xarray as xr

        a=xr.open_dataset( ifile )
        ds=xr.Dataset( coords= a.coords )

        lvar=list( a.variables )
        lons = a.lon.values
        lats = a.lat.values

        ii = np.argmin(  abs(lons-ilon) )
        jj = np.argmin(  abs(lats-ilat) )

        print(lvar)


        for var in lvar:
            aa=a[var]
            print(var,aa.values.dtype )
            if( (aa.values.dtype is np.dtype('float32') ) or (aa.values.dtype is np.dtype('float64') )
              or (aa.values.dtype is np.dtype('int') ) ):
                dims=aa.dims
                coords=aa.coords
                print("Trying zonal mean of "+var )
                print("with Dims",dims )
                if ( ('lon' in dims) and ('lat' in dims) and ('lev' in dims or 'ilev' in dims) ):
                    # 3D var 
                    aa2 = aa[:,:,jj,ii]
                    ds[var]=aa2
                elif ( ('lon' in dims) and ('lat' in dims) and ('lev' not in dims) and ('ilev' not in dims) ):
                    # 2D var 
                    aa2 = aa[:,jj,ii]
                    ds[var]=aa2
                else:
                    # something else ??? slon? slat?
                    print( var, "is not pickable")
            else:
                ds[var] = aa
                print( var, "is not pickable")

        ds.to_netcdf( ofile )
        return ds

    def zonal_mean_lonfill(self, ifile , ofile ):
        import numpy as np
        import xarray as xr

        a=xr.open_dataset( ifile )
        ds=xr.Dataset( coords= a.coords )


        #b=xr.open_dataset( ofile )
        za = a.mean(dim='lon')
        dims=a.dims
        if ('slon' in dims):
            za = za.mean( dim='slon')

        lvar=list( a.variables )

        if('lon' in lvar):
            londar=a['lon']
            ds['lon']=londar
            lvar.remove('lon')
        if('slon' in lvar):
            slondar=a['slon']
            ds['slon']=slondar
            lvar.remove('slon')


        print(lvar)

        lvar = ['lat', 'gw', 'lev', 'hyam', 'hybm', 'P0', 
                'ilev', 'hyai', 'hybi', 'time', 'date', 'datesec', 
                'time_bnds', 'date_written', 'time_written', 'ndbase', 
                'nsbase', 'nbdate', 'nbsec', 'mdt', 'ndcur', 'nscur', 
                'co2vmr', 'ch4vmr', 'n2ovmr', 'f11vmr', 'f12vmr', 'sol_tsi', 
                'nsteph', 'OMEGA', 'PRECC', 'PRECL', 'PS', 'Q', 'T', 
                'U', 'UTEND_CORE', 'UTEND_PHYSTOT', 
                'V', 'VTEND_CORE', 'VTEND_PHYSTOT', 
                'Z3']

        for var in lvar:
            aa=a[var]
            print(var,aa.values.dtype )
            if( (aa.values.dtype is np.dtype('float32') ) or (aa.values.dtype is np.dtype('float64') )
              or (aa.values.dtype is np.dtype('int') ) ):
                aaz=za[var]
                aa2=aa  #*0.
                dims=aa.dims
                coords=aa.coords
                print("Trying zonal mean of "+var )
                print("with Dims",dims )
                if ( ('slon' in dims) and ('lev' in dims or 'ilev' in dims) ):
                    n=aa.slon.size
                    for i in np.arange(n):
                        aa2[:,:,:,i]=aaz
                elif ( ('lon' in dims) and ('lev' in dims or 'ilev' in dims) ):
                    n=aa.lon.size
                    for i in np.arange(n):
                        aa2[:,:,:,i]=aaz
                elif ( ('lon' in dims) and ('lev' not in dims) and ('ilev' not in dims) ):
                    n=aa.lon.size
                    for i in np.arange(n):
                        aa2[:,:,i]=aaz
                elif ( ('slon' in dims) and ('lev' not in dims) and ('ilev' not in dims) ):
                    n=aa.slon.size
                    for i in np.arange(n):
                        aa2[:,:,i]=aaz

                ds[var] = aa2
            else:
                ds[var] = aa
                print( var, "is not averageable")


            ds.to_netcdf( ofile )



    def zonal_mean(self, ifile , ofile ):
        import numpy as np
        import xarray as xr

        a=xr.open_dataset( ifile )

        za = a.mean(dim='lon')
        dims=a.dims
        if ('slon' in dims):
            za = za.mean( dim='slon')



        za.to_netcdf( ofile )

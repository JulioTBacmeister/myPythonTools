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

        self.file=file

        return file


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

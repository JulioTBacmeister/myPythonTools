class camdata:

    import sys
    sys.path.append('../Plotting/')
    
    import numpy as np


    def __init__(self):
        
        self.filenames="filenames"

    def make_pmid( self, a):
        import numpy as np
        import xarray as xr

        """ Initialize possible dimension sizes """
        ncol=-999
        nlat=-999
        nlon=-999
        nlev=-999
        ntime=0

        """ Sort out structure of dataset 'a' """
        dims=a.dims
        if ('time' in dims):
            ntime=dims['time']
        if ('ntime' in dims):
            ntime=dims['ntime']
        if ('lev' in dims):
            nlev=dims['lev']
        if ('nlev' in dims):
            nlev=dims['nlev']
        if ('ncol' in dims):
            ncol=dims['ncol']
        if ('ncol' not in dims and 'lat' in dims):
            nlat=dims['lat']
        if ('ncol' not in dims and 'nlat' in dims):
            nlat=dims['nlat']
        if ('ncol' not in dims and 'lon' in dims):
            nlon=dims['lon']
        if ('ncol' not in dims and 'nlon' in dims):
            nlon=dims['nlon']

        PS=a['PS']
        lp3dims=list(PS.dims)
        print(lp3dims)
        if ('time' in lp3dims or 'ntime' in lp3dims):
            if ('lev' in dims):
                lp3dims.insert(1,'lev')
            if ('nlev' in dims):
                lp3dims.insert(1,'nlev')
        print(lp3dims)
        if ('time' not in lp3dims and 'ntime' not in lp3dims):
            if ('lev' in dims):
                lp3dims.insert(0,'lev')
            if ('nlev' in dims):
                lp3dims.insert(0,'nlev')
        
        hyam=a['hyam']
        hybm=a['hybm']

        lon=a['lon']

        print(dims)
        print(lp3dims)

        ndims=len( lp3dims )
        dimsizes=np.zeros( ndims,dtype='int' )
        for i in np.arange( ndims ):
            dimsizes[i]=int( dims[lp3dims[i]] )

        print(ndims , dimsizes )

        """ 
        Acceptable strutcures for PMID variable 
        ndims (2)  = ['{lev','nlev}, 'ncol']
        ndims (3)  = [{'time','ntime'}, '{lev','nlev}, 'ncol']
        ndims (3)  = ['{lev','nlev}, {'lat,'nlat'}', {'lon','nlon'} ]
        ndims (4)  = [{'time','ntime'}, '{lev','nlev}, {'lat,'nlat'}', {'lon','nlon'} ]
        """

        PMIDv = np.zeros( dimsizes , dtype='float' )
        print("PMID", PMIDv.shape )

        assert( lp3dims[0]=='time' or lp3dims[0]=='ntime' or lp3dims[0]=='lev' or lp3dims[0]=='nlev' )

        if ( lp3dims[0]=='lev' or lp3dims[0]=='nlev' ):
            if (ndims==2):
                for L in np.arange(nlev):
                    PMIDv[L,:] =  hyam[L]*100000. + hybm[L]*PS[:]
            if (ndims==3):
                for L in np.arange(nlev):
                    PMIDv[L,:,:] =  hyam[L]*100000. + hybm[L]*PS[:,:]
        if ( lp3dims[0]=='time' or lp3dims[0]=='ntime' ):
            if (ndims==3):
                for n in np.arange(ntime):
                    for L in np.arange(nlev):
                        PMIDv[n,L,:] =  hyam[L]*100000. + hybm[L]*PS[n,:]
            if (ndims==4):
                for n in np.arange(ntime):
                    for L in np.arange(nlev):
                        PMIDv[n,L,:,:] =  hyam[L]*100000. + hybm[L]*PS[n,:,:]

        PMID = xr.DataArray( data=PMIDv, dims=tuple(lp3dims) )
        #a['PMID']=PMID

        return  PMID

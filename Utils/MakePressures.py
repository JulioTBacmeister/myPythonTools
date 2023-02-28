import numpy as np


def Pressure ( am, bm, ai, bi, ps , p_00=100_000., Gridkey='tzc' ):
    # Should ASSERT that am ,bm, ai, and bi are 1D

    if ( Gridkey == 'tzc' ):
        nt, ncol = np.shape( ps )
        nz = np.shape( am )[0]
        delp = np.zeros( (nt,nz,ncol ) )
        pmid = np.zeros( (nt,nz,ncol ) )
        pint = np.zeros( (nt,nz+1,ncol ) )
        for i in np.arange( nt ):
            for L in np.arange( nz ):
                pmid[i,L,:] = am[L]*p_00 + bm[L]*ps[i,:]
            for L in np.arange( nz+1 ):
                pint[i,L,:] = ai[L]*p_00 + bi[L]*ps[i,:]
            for L in np.arange( nz ):
                delp[i,L,:] = pint[i,L+1,:]-pint[i,L,:]
                        
 
    return pmid,pint,delp
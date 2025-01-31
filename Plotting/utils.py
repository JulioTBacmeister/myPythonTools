import numpy as np
def contour_intervals( fld, **kwargs ):

    if (fld=='U'):
        if 'zonal_mean' in kwargs:
            flev=np.linspace( -60,140,num=21) 
            dlev=np.linspace( -20,20,num=21) 
            scale=1.0
    if (fld=='Nudge_U'):
        if 'zonal_mean' in kwargs:
            flev=np.linspace( -20,20,num=21) 
            dlev=np.linspace( -20,20,num=21) 
            scale=86_400.
    if (fld=='UTEND_GWDTOT'):
        if 'zonal_mean' in kwargs:
            flev=np.linspace( -20,20,num=21) 
            dlev=np.linspace( -20,20,num=21) 
            scale=86_400.
    else:
        flev=21
        dlev=21
        scale = 1.0
        
    return flev,dlev,scale

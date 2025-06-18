import numpy as np
def contour_intervals( fld, **kwargs ):

    tendlev = [-50, -30, -10, -5, -3, -1, -.5, -.3, -.1,  -.05, -.03, -.01, -.005, 0,  .005, .01, .03, .05,  .1, .3, .5, 1, 3, 5, 10, 30, 50 ]

    if (fld=='T'):
        if 'zonal_mean' in kwargs:
            flev=np.linspace( 180,300,num=21) 
            dlev=np.linspace( -20,20,num=21) 
            scale=1.0
    elif (fld=='U'):
        if 'zonal_mean' in kwargs:
            flev=np.linspace( -60,140,num=21) 
            dlev=np.linspace( -20,20,num=21) 
            scale=1.0
        if 'yx' in kwargs:
            if 'plev' in kwargs:
                plev=kwargs['plev']  #This should be in hPa
                if ((plev >=50.) and (plev<500.)):
                    flev=np.linspace( -60,140,num=21) 
                    dlev=np.linspace( -20,20,num=21) 
        scale=1.0
    elif (fld=='V'):
        if 'zonal_mean' in kwargs:
            flev=0.5 * np.linspace(-5,5,num=21)
            dlev=np.linspace( -1,1,num=21) 
            scale=1.0
        if 'yx' in kwargs:
            if 'plev' in kwargs:
                plev=kwargs['plev']  #This should be in hPa
                if ((plev >=50.) and (plev<500.)):
                    flev=np.linspace( -60,140,num=21) 
                    dlev=np.linspace( -20,20,num=21) 
        scale=1.0
    elif (fld[0:2]=='U_'):
        if 'zonal_mean' in kwargs:
            flev=np.linspace( -60,140,num=21) 
            dlev=np.linspace( -20,20,num=21) 
            scale=1.0
        if 'yx' in kwargs:
            if 'plev' in kwargs:
                plev=kwargs['plev']  #This should be in hPa
                if ((plev >=50.) and (plev<500.)):
                    flev=np.linspace( -60,140,num=21) 
                    dlev=np.linspace( -20,20,num=21) 
        scale=1.0
    elif (fld=='PRECT'):
        flev= np.linspace(0,20,num=21)
        dlev=np.linspace(-5,5,num=21)
        scale=86_400. * 1_000.
    elif (fld=='TAUX'):
        flev= np.linspace(-0.3,0.3,num=21)
        dlev=np.linspace(-0.06,0.06,num=21)
        scale=1.
    elif (fld=='U10'):
        flev= np.linspace(0,16,num=21)
        dlev=np.linspace(-2,2,num=21)
        scale=1.
    elif (fld[0:2]=='dU'):
        if 'zonal_mean' in kwargs:
            flev= tendlev
            dlev=tendlev
        scale=86_400.
    elif (fld=='Nudge_U'):
        if 'zonal_mean' in kwargs:
            flev= tendlev
            dlev=tendlev
        scale=86_400.
    elif (fld=='UTEND_GWDTOT'):
        if 'zonal_mean' in kwargs:
            flev=tendlev
            dlev=tendlev
        scale=86_400.
    elif (fld=='UTGWORO'):
        if 'zonal_mean' in kwargs:
            flev=tendlev
            dlev=tendlev
        scale=86_400.
    elif (fld=='UTGW_MOVMTN'):
        if 'zonal_mean' in kwargs:
            flev=tendlev
            dlev=tendlev
        scale=86_400.
    elif (fld=='UTGWSPEC'):
        if 'zonal_mean' in kwargs:
            flev=tendlev
            dlev=tendlev
        scale=86_400.
    elif (fld=='UTEND'):
        if 'zonal_mean' in kwargs:
            flev=tendlev
            dlev=tendlev
            scale=86_400.
    elif (fld[0:5]=='UTEND'):
        if 'zonal_mean' in kwargs:
            flev=tendlev
            dlev=tendlev
            scale=86_400.
    elif (fld[0:5]=='VQ'):
        if 'zonal_mean' in kwargs:
            flev=np.linspace( -0.05,0.05,num=41) 
            dlev=flev/5.
            scale=86_400.
    else:
        print( f" No presets found for {fld}. Using 21's ")
        flev=21
        dlev=21
        scale = 1.0
        
    return flev,dlev,scale

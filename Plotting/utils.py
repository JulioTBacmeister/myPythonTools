import numpy as np
def contour_intervals( fld, **kwargs ):

    tendlev = [-50, -30, -10, -5, -3, -1, -.5, -.3, -.1,  -.05, -.03, -.01, -.005, 0,  .005, .01, .03, .05,  .1, .3, .5, 1, 3, 5, 10, 30, 50 ]
    cmap = 'gist_ncar'

    if (fld=='T'):
        if 'zonal_mean' in kwargs:
            flev=np.linspace( 180,300,num=21) 
            dlev=np.linspace( -20,20,num=21) 
            scale=1.0
    elif (fld=='U') or (fld=='V'):
        if 'zonal_mean' in kwargs:
            flev=np.linspace( -60,140,num=21) 
            dlev=np.linspace( -20,20,num=21) 
        elif 'yx' in kwargs:
            if 'plev' in kwargs:
                plev=kwargs['plev']  #This should be in hPa
                if ((plev >=50.) and (plev<500.)):
                    flev=np.linspace( -60,140,num=21) 
                    dlev=np.linspace( -20,20,num=21) 
            if 'zlev' in kwargs:
                zlev=kwargs['zlev']  #This should be in hPa
                if ((zlev >=20.) and (zlev<90.)):
                    flev=np.linspace( -60,140,num=21) 
                    dlev=np.linspace( -20,20,num=21) 
                elif ((zlev >=6.) and (zlev<20.)):
                    flev=np.linspace( -60,140,num=21) 
                    dlev=np.linspace( -20,20,num=21) 
                elif ((zlev >=1.) and (zlev<6.)):
                    flev=0.5*np.linspace( -60,140,num=21) 
                    dlev=0.5*np.linspace( -20,20,num=21) 
                elif ((zlev >=0.5) and (zlev<1.)):
                    flev=0.25*np.linspace( -60,140,num=21) 
                    dlev=0.25*np.linspace( -20,20,num=21) 
                elif ((zlev >=0.0) and (zlev<0.5)):
                    flev=np.linspace( -7.,13.,num=21) 
                    dlev=np.linspace( -5,5,num=21) 
                else: 
                    flev=np.linspace( -60,140,num=21) 
                    dlev=np.linspace( -20,20,num=21) 
            else:
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
    elif (fld=='RESTOM'):
        flev= np.linspace(-120,120,num=25)
        dlev=np.linspace(-50,50,num=21)
        scale=1.0
    elif (fld=='PRECT'):
        flev= np.linspace(0,20,num=21)
        dlev=np.linspace(-5,5,num=21)
        scale=86_400. * 1_000.
        cmap='Blues_r'
    elif (fld=='LHFLX'):
        flev= np.linspace(-20,240,num=27)
        dlev=np.linspace(-30,30,num=31)
        scale=1.
    elif (fld=='TAUX'):
        flev= np.linspace(-0.2,0.2,num=21)
        dlev=np.linspace(-0.06,0.06,num=21)
        scale=1.
    elif (fld=='TAUY'):
        flev= np.linspace(-0.3,0.3,num=21)
        dlev=np.linspace(-0.06,0.06,num=21)
        scale=1.
    elif (fld=='U10'):
        flev= np.linspace(0,16,num=21)
        dlev=np.linspace(-2,2,num=21)
        scale=1.
    elif (fld=='QRS'):
        if 'zonal_mean' in kwargs:
            flev= 21 #0.1*np.linspace(-20,20,num=21) 
            dlev=21 #0.1*np.linspace(-20,20,num=21) 
        scale=86_400.
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
        else:
            flev = np.linspace(-30,30,num=31)
            dlev=  0.2*flev
        scale=86_400.
    elif (fld=='UTGW_MOVMTN'):
        if 'zonal_mean' in kwargs:
            flev=tendlev
            dlev=tendlev
        elif 'yx' in kwargs:
            if 'plev' in kwargs:
                plev=kwargs['plev']  #This should be in hPa
                flev = np.linspace(-30,30,num=31)
                dlev=  0.2*flev
            if 'zlev' in kwargs:
                zlev=kwargs['zlev']  #This should be in km
                mxv =  2.+0.4*np.exp( zlev / 15.0 )
                mxv = 2* mxv//2
                flev = np.linspace(-mxv,mxv,num=31)
                dlev=  0.2*flev
        else:
            flev = np.linspace(-30,30,num=31)
            dlev=  0.2*flev
        scale=86_400.
    elif (fld=='UTGWSPEC'):
        if 'zonal_mean' in kwargs:
            flev=tendlev
            dlev=tendlev
        elif 'yx' in kwargs:
            if 'plev' in kwargs:
                plev=kwargs['plev']  #This should be in hPa
                flev = np.linspace(-30,30,num=31)
                dlev=  0.2*flev
            if 'zlev' in kwargs:
                zlev=kwargs['zlev']  #This should be in km
                mxv =  2.+0.4*np.exp( zlev / 15.0 )
                mxv = 2* mxv//2
                flev = np.linspace(-mxv,mxv,num=31)
                dlev=  0.2*flev
        else:
            flev = np.linspace(-30,30,num=31)
            dlev=  0.2*flev
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
    
    if ('return_cmap' in kwargs):
        return flev,dlev,scale,cmap
    else:
        return flev,dlev,scale


#####################
#
# to ChatGPT
# SO, again if I have defined a dict, how do I enable the dict.key vs dict['key'] functionality? 
#  ...
#  ... code below
# Your dictionary definitions
# my_dict = AttrDict({"name": "John", "age": 30})

"""
The first AttrDict implementation has a clear separation between the underlying dictionary data (_data) and the class's attributes. This separation can be useful in certain situations where you want to avoid potential conflicts between dictionary keys and class attributes. For example, if you have a dictionary key that matches a method or attribute name of the AttrDict class, the first implementation would prevent any conflicts, as the dictionary data is stored separately in the _data attribute.

However, for most practical purposes, the second implementation is more straightforward and convenient, as it allows direct access to the dictionary keys as attributes while still retaining all the functionality of a regular dictionary. It's also more concise and easier to understand for those who are familiar with Python's built-in dict class.

In summary, the first implementation provides a clear separation between class attributes and dictionary keys, which can be useful in certain edge cases, while the second implementation is generally more convenient and straightforward for everyday use.
"""
"""
# This doesn't allow both dict.key and dict['key'] syntax
class AttrDict:
    def __init__(self, data=None):
        if data is None:
            data = {}
        self.__dict__["_data"] = data

    def __getattr__(self, key):
        if key in self._data:
            return self._data[key]
        raise AttributeError(f"'AttrDict' object has no attribute '{key}'")

    def __setattr__(self, key, value):
        self._data[key] = value

    def __delattr__(self, key):
        if key in self._data:
            del self._data[key]
        else:
            raise AttributeError(f"'AttrDict' object has no attribute '{key}'")
"""
# This allow both dict.key and dict['key'] syntax
class AttrDict(dict):
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(f"'AttrDict' object has no attribute '{key}'")

    def __setattr__(self, key, value):
        self[key] = value

    def __delattr__(self, key):
        try:
            del self[key]
        except KeyError:
            raise AttributeError(f"'AttrDict' object has no attribute '{key}'")


def run(code):
    #import time
    print("exec(open(" +"'"+code+"'"+ ").read())")
    #print("will sleep in ",code)
    #time.sleep(10)

def env():
    import matplotlib.pyplot as plt

def HowManyWorkers(): 
    import os
    # set number of pools to create
    #nworkers = cpu_count()
    nworkers = len(os.sched_getaffinity(0))
    print('Flexible VertRegrid using reshaped ARRAYS ')
    print('nworkers available ',nworkers)

    return nworkers


def MakePath( user='juliob' , exp='YaYa', subd='hist', hsPat='cam.h0' , ymdPat='*' ,verbose=False):

    if (user in ['juliob','pel','tilmes'] ):
        path = f'/glade/derecho/scratch/{user}/archive/{exp}/atm/{subd}/{exp}.{hsPat}.{ymdPat}.nc' 
    elif (user == 'juliob_run' ):
        path = f'/glade/derecho/scratch/juliob/{exp}/run/{exp}.cam.h0.*.nc' 
    elif (user == 'juliob_camp' ):
        path = f'/glade/campaign/cgd/amp/juliob/archive/{exp}/atm/{subd}/{exp}.{hsPat}.{ymdPat}.nc'   
    elif (user == 'CMIP6' ):
        dir='/glade/campaign/collections/cmip/CMIP6/timeseries-cmip6/f.e21.FWHISTBgcCrop.f09_f09_mg17.CMIP6-AMIP-WACCM.001/atm/proc/tseries/month_1/'
        path=dir+ 'f.e21.FWHISTBgcCrop.f09_f09_mg17.CMIP6-AMIP-WACCM.001.cam.h0.U.195001-201412.nc'


    if (verbose==True):
        print( path )


    return path

def MakeDict4Exp( user='juliob' , exp='YaYa', subd='hist', hsPat='cam.h0' , ymdPat='*' ,verbose=False, open_dataset=False, help=False ):
    import xarray as xr

    if (help == True):
        print( f" Possible 'users': 'juliob','pel','tilmes''juliob_run' 'juliob_camp' 'amwg_runs'  'omwg_mom6' 'CMIP6' ")
        return
    
    if (user in ['juliob','pel','tilmes'] ):
        path = f'/glade/derecho/scratch/{user}/archive/{exp}/atm/{subd}/{exp}.{hsPat}.{ymdPat}.nc' 
    elif (user == 'juliob_run' ):
        path = f'/glade/derecho/scratch/juliob/{exp}/run/{exp}.{hsPat}.{ymdPat}.nc' 
    elif (user == 'juliob_camp' ):
        path = f'/glade/campaign/cgd/amp/juliob/archive/{exp}/atm/{subd}/{exp}.{hsPat}.{ymdPat}.nc'   
    elif (user == 'amwg_runs' ):
        path = f'/glade/campaign/cgd/amp/amwg/runs/{exp}/atm/{subd}/{exp}.{hsPat}.{ymdPat}.nc'   
    elif (user == 'omwg_mom6' ):
        path = f'/glade/campaign/cesm/development/omwg/projects/MOM6/{exp}/atm/{subd}/{exp}.{hsPat}.{ymdPat}.nc'   
    elif (user == 'CMIP6' ):
        dir='/glade/campaign/collections/cmip/CMIP6/timeseries-cmip6/f.e21.FWHISTBgcCrop.f09_f09_mg17.CMIP6-AMIP-WACCM.001/atm/proc/tseries/month_1/'
        path=dir+ 'f.e21.FWHISTBgcCrop.f09_f09_mg17.CMIP6-AMIP-WACCM.001.cam.h0.U.195001-201412.nc'

    ###/glade/u/home/gmarques/campaign_MOM6/b.cesm3_cam058_mom_e.B1850WscMOM.ne30_L58_t061.camdev_cice5.026g
    ### /glade/campaign/cesm/development/omwg/projects/MOM6/

    if (verbose==True):
        print( path )

    dex = { 'exp':exp ,
           'user':user, 
           'subd':subd,
           'hsPat':hsPat,
           'ymdPat':ymdPat,
           'path':path }
    
    if (open_dataset==True):
        X = xr.open_mfdataset( path ,data_vars='different', coords='different' )
        dex['X']=X

    # Use class defined above to enable 'dict.key' syntax.
    #--------------------------------
    Adex=AttrDict( dex )
    
    return Adex

def extValues( DX , flds ):

    dex = { 'exp':exp ,
           'user':user, 
           'subd':subd,
           'hsPat':hsPat,
           'ymdPat':ymdPat,
           'path':path }

    fldx = ['lat','lon','lev']+flds

    for fld in fldx:
        dex[ fld ] = DX.X[fld].values

    return dex

def trim_to_year( D , nyr_max=1000 ):
    import numpy as np
    
    nyr=np.minimum( nyr_max*12 , 12*(len(D.X.time)//12) )
    print(f"Discarding last {len(D.X.time)-nyr} months of {D.exp}" )
    D.X=D.X.isel(time=np.arange(nyr))




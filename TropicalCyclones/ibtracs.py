import xarray as xr
import numpy as np

import importlib
import copy

import cftime
import time as time
import datetime

#===========================================================
# Class to allow things to be accessed via dict.thing syntax
# as well as dict['thing'] syntax. They are equivalent.
#===========================================================
class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


#===========================================================
#                 Begin functions 
#===========================================================
def readtrx(Fill_NI_before_1990=True,Add_Precip=False):
    
    f="/glade/work/juliob/IBTrACS/IBTrACS.since1980.v04r00.nc"
    dS=xr.open_dataset( f )
        
    cbasin=dS.basin.values
    timevar=dS.time.values

    timevar_A = np.datetime_as_string( timevar )

    
    nstorm,nt=np.shape(timevar_A)
    year,month,day,hour = np.zeros((nstorm,nt) , dtype='int'),\
                        np.zeros((nstorm,nt) , dtype='int'),\
                        np.zeros((nstorm,nt) , dtype='int'),\
                        np.zeros((nstorm,nt) , dtype='int')

    for s in np.arange(nstorm):
        for t in np.arange( nt ):
            if ( timevar_A[s,t] != 'NaT'):
                year[s,t]  = int(timevar_A[s,t][0:4])
                month[s,t] = int(timevar_A[s,t][5:7])
                day[s,t]   = int(timevar_A[s,t][8:10])
                hour[s,t]  = int(timevar_A[s,t][11:13])
            else:
                year[s,t]  = -99999
                
    
    
    wind_kts_wmo = dS.wmo_wind.values # IBTrACS wind in kts
    wind_kts_USA = dS.usa_wind.values # IBTrACS wind in kts
    
    if ( Fill_NI_before_1990 == True):
        usa_lon=dS.usa_lon.values
        usa_lat=dS.usa_lat.values
        NIbasin=np.where( (usa_lat >0.) & (usa_lat<30.) & (usa_lon<105.) & (usa_lon > 40.) , True, False )
        print( "USING USA winds in NI before 1990" )
        wind_kts = np.where( ( year<1990 )& (NIbasin==True)  
                              , wind_kts_USA, wind_kts_wmo )
    else:
        print( "No winds in NI before 1990")
        wind_kts = wind_kts_wmo
    
    wind=wind_kts/1.944  # IBTrACS wind is in kts. This makes it m/s.

    wind=np.nan_to_num( wind , nan=0.) # MAny W Pac 
    
    # [b'' b'EP' b'NA' b'NI' b'SA' b'SI' b'SP' b'WP']
    
    basin = np.zeros((nstorm,nt) , dtype='int')-1
    basin = np.where( cbasin==b'NA' , 0 , basin)
    basin = np.where( cbasin==b'SA' , 1 , basin)
    basin = np.where( cbasin==b'WP' , 2 , basin)
    basin = np.where( cbasin==b'EP' , 3 , basin)
    basin = np.where( cbasin==b'SP' , 4 , basin)
    basin = np.where( cbasin==b'NI' , 5 , basin)
    basin = np.where( cbasin==b'SI' , 6 , basin)
    
                
    struc = {'year':year,'month':month,'day':day,'hour':hour,
            'basin':basin, 'lon':dS.lon.values , 'lat':dS.lat.values ,
            'wind':wind , 'pres':dS.wmo_pres.values }

    if (Add_Precip==True):
        npSaveFile = '/glade/work/juliob/NumPySaves/PrecTrax/TRMM_IBTrACS-prectrax-1998-2019.npz'
        print(f"File = {npSaveFile} ")
        SaveData = np.load( npSaveFile )
        prectrk = SaveData['prectrk']
        struc['prectrk']=prectrk

    Atruc = AttrDict( struc )
    
    return Atruc

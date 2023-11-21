import xarray as xr
import numpy as np

import importlib
import copy

import trax_util as trx

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
    #f="/glade/work/juliob/IBTrACS/IBTrACS.ALL.v04r00.nc"
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
        ####                      , wind_kts_USA, wind_kts_wmo  )
                              , wind_kts_USA, wind_kts_USA  )
    else:
        print( "No winds in NI before 1990")
        wind_kts = wind_kts_wmo
    
    wind=wind_kts/1.944  # IBTrACS wind is in kts. This makes it m/s.

    wind=np.nan_to_num( wind , nan=0.) # MAny W Pac 
    
    speed = trx.translation_speed( dS.lon.values , dS.lat.values ) 
    intfW,intfP = trx.intensification( wind , dS.wmo_pres.values ) 

    
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
            'wind':wind , 'speed':speed, 'pres':dS.wmo_pres.values ,
            'intfW':intfW, 'intfP':intfP }

    if (Add_Precip==True):
        npSaveFile = '/glade/work/juliob/NumPySaves/PrecTrax/TRMM_IBTrACS-prectrax-1998-2019.npz'
        print(f"File = {npSaveFile} ")
        SaveData = np.load( npSaveFile )
        prectrk = SaveData['prectrk']
        struc['prectrk']=prectrk

    Atruc = AttrDict( struc )
    
    return Atruc

########################################################
def basin_outlines(**kwargs):

    basins=[]
    bb={'name':'None','number':-1,'xy':[-9999,-9999]}
    basins.append( AttrDict (bb) )

    if ('Mediterranean' in kwargs):
        addMed=kwargs['Mediterranean']
    else:
        addMed=False
    if ('xNEPacific' in kwargs):
        addxNEPac=kwargs['xNEPacific']
    else:
        addaNEPac=False

    
    shf=360.
    # NAtl
    xybasin=[]
    xybasin.append([-100 + shf, 60])
    xybasin.append([-100 + shf, 18])
    xybasin.append([-90  + shf, 18])
    xybasin.append([-90  + shf, 16])
    xybasin.append([-85  + shf, 16])
    xybasin.append([-85  + shf, 8])
    xybasin.append([-75  + shf, 8])
    xybasin.append([-75  + shf, 5])
    xybasin.append([ 0  + shf, 5])
    ## Notch in Mediterranean
    if ( addMed==True ):
            xybasin.append([ 0  + shf, 34])
            xybasin.append([ -6  + shf, 34])
            xybasin.append([ -6  + shf, 39])
            xybasin.append([ 0  + shf, 39])
    ###
    xybasin.append([ 0  + shf, 60] )
    xybasin.append([-100 + shf, 60])
    bb={'name':'N Atlantic','number':0,'xy':xybasin, 'labloc':[262,52] }
    #basins.append(bb)
    basins.append( AttrDict (bb) )

    # SAtl
    xybasin=[]
    xybasin.append([ 290  , -60] )
    xybasin.append([ 290 , -5])
    xybasin.append([ 360 , -5])
    xybasin.append([ 360 , -60])
    xybasin.append([ 290 , -60])

    bb={'name':'S Atlantic','number':1,'xy':xybasin , 'labloc':[292,-58] }
    basins.append( AttrDict (bb) )

    
    # NWPac
    xybasin=[]
    xybasin.append([ 180  , 60] )
    xybasin.append([ 100 , 60])
    xybasin.append([ 100 , 5])
    xybasin.append([ 180 , 5])
    xybasin.append([ 180 , 60])
    bb={'name':'NW Pacific','number':2,'xy':xybasin, 'labloc':[102,52] }
    #basins.append(bb)
    basins.append( AttrDict (bb) )
    
    # NEPac
    xybasin=[]
    xybasin.append([-100 + shf, 60])
    xybasin.append([-100 + shf, 18])
    xybasin.append([-90  + shf, 18])
    xybasin.append([-90  + shf, 16])
    xybasin.append([-85  + shf, 16])
    xybasin.append([-85  + shf, 8])
    xybasin.append([-75  + shf, 8])
    xybasin.append([-75  + shf, 5])
    xybasin.append([ -180  + shf, 5])
    xybasin.append([-180  + shf, 60] )
    xybasin.append([-100 + shf, 60])
    bb={'name':'NE Pacific','number':3,'xy':xybasin, 'labloc':[182,52] }
    #basins.append(bb)
    basins.append( AttrDict (bb) )
    
    
    # S Pac
    xybasin=[]
    xybasin.append([ 135  , -5] )
    xybasin.append([ 240  , -5])
    xybasin.append([ 240, -60])
    xybasin.append([ 135, -60])
    xybasin.append([ 135, -5])
    bb={'name':'S Pacific','number':4,'xy':xybasin, 'labloc':[137,-58] }
    #basins.append(bb)
    basins.append( AttrDict (bb) )

    
    # NInd
    xybasin=[]
    xybasin.append([ 100 , 60])
    xybasin.append([ 100 , 5])
    xybasin.append([ 30 , 5])
    ## Notch in Mediterranean
    if ( addMed==True ):
            xybasin.append([ 30 , 28])
            xybasin.append([ 38.5 , 28])
            xybasin.append([ 38.5 , 38.5])
            xybasin.append([ 30 , 38.5])
    ###
    xybasin.append([ 30 , 60])
    xybasin.append([ 100 , 60])
    bb={'name':'N Indian','number':5,'xy':xybasin, 'labloc':[32,52] }
    #basins.append(bb)
    basins.append( AttrDict (bb) )
    
    # S Ind
    xybasin=[]
    xybasin.append([ 135  , -5] )
    xybasin.append([ 10  , -5])
    xybasin.append([ 10, -60])
    xybasin.append([ 135, -60])
    xybasin.append([ 135, -5])
    bb={'name':'S Indian','number':6,'xy':xybasin, 'labloc':[12,-58] }
    #basins.append(bb)
    basins.append( AttrDict (bb) )
    bnumber=6
    
    # Med
    if ( addMed==True ):
        bnumber=bnumber+1
        xybasin=[]
        xybasin.append([ 0  , 34])
        xybasin.append([ -6 , 34])
        xybasin.append([ -6 , 39])
        xybasin.append([ 0  , 39])
        xybasin.append([ 0  , 47])
        xybasin.append([ 30  , 47])
        xybasin.append([ 30  , 38.5])
        xybasin.append([ 38.5  , 38.5])
        xybasin.append([ 38.5  , 28])
        xybasin.append([ 0  , 28])
        xybasin.append([ 0  , 34])
        bb={'name':'Med.','number':7,'xy':xybasin, 'labloc':[5,18] }
        #basins.append(bb)
        basins.append( AttrDict (bb) )

    # extreme NE Pacific
    if ( addxNEPac==True ):
        bnumber=bnumber+1
        xybasin=[]
        xybasin.append([-100 + shf, 24])
        xybasin.append([-100 + shf, 18])
        xybasin.append([-90  + shf, 18])
        xybasin.append([-90  + shf, 16])
        xybasin.append([-85  + shf, 16])
        xybasin.append([-85  + shf, 8])
        xybasin.append([-75  + shf, 8])
        xybasin.append([-75  + shf, 5])
        xybasin.append([ -112  + shf, 5])
        xybasin.append([-112 + shf, 24] )
        xybasin.append([-100 + shf, 24])
        bb={'name':'xNE Pacific','number':31,'xy':xybasin, 'labloc':[360-115,32] }
        #basins.append(bb)
        basins.append( AttrDict (bb) )

    
    #basinX = AttrDict( basins ) 
    
    return basins
    
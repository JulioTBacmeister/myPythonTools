import sys
# import modules in other directories
sys.path.append('../Utils/')

import numpy as np
import xarray as xr

#import MyConstants as Con
import constants as Con
import trax_util as trx


#===========================================================
# Class to allow things to be accessed via dict.thing syntax
# as well as dict['thing'] syntax. They are equivalent.
#===========================================================
class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


#===========================================================
#               functions
#===========================================================

# ---------------------------
# --- Basic reader
# ---------------------------
def parse_storm_file(filename, columns=None ):
    storm_list = []
    current_storm = None
    #,lon,lat,slp,wind,wind850,phis"
    if columns is None:
        columns_ = ['cell_id','lon','lat','slp','wind','wind850','phis','year','month','day','hour',]
    else:
        columns_ = columns
    
    with open(filename, 'r') as file:
        for line in file:
            # Remove leading/trailing whitespaces
            line = line.strip()
            
            if line.startswith('start'):
                # If we encounter 'start', we finish the current storm dict and start a new one
                if current_storm is not None:
                    storm_list.append(current_storm)
                
                # Initialize a new storm dictionary
                parts = line.split()
                current_storm = {}
                for f in columns_:
                    current_storm[f]=[]

                continue
            
            # For non-'start' lines, read the data points
            parts = line.split()


            ncol=len(columns_)
            for i in np.arange(ncol):
                if columns_[i] in ( 'year','month','day','hour', ):
                    current_storm[ columns_[i] ].append( int(parts[i])  )
                else:               
                    current_storm[ columns_[i] ].append( float(parts[i])  )
            
        # Append the last storm
        if current_storm is not None:
            storm_list.append(current_storm)
    
    return storm_list
    
# ---------------------------
# --- translate to 'old' CAM5 style
# ---------------------------
def readtrx( fileN ,columns=None, power_wind=1.0):
    
    dS = parse_storm_file( fileN, columns=columns )
    nstorm=len(dS)

    #--------------------------------------------
    # first let's find the longest storm record
    #--------------------------------------------
    duration=np.zeros( nstorm, dtype=int )
    for i in np.arange( nstorm ):
        duration[i] = len(dS[i]['year'])
    nt = np.max( duration )

    empty=np.zeros( (nstorm, nt) )-999.
    
    #wind = dS.wind.values * power_wind

    struc={}
    # Lists of matching fields in Tempest and CAM% track files (modulo small name changes)
    # e.g. pres==>slp
    tmpxt_flds=['year', 'month', 'day', 'hour', 'lon', 'lat', 'wind', 'slp', ]
    match_flds=['year', 'month', 'day', 'hour', 'lon', 'lat', 'wind', 'pres', ]

    for tfld, mfld in zip(tmpxt_flds, match_flds, strict=True):
        data = create_2d_storm_array ( nt=nt ,nstorm= nstorm, dS=dS, fld=tfld ) 
        struc[mfld]=data

    struc['date'] = struc['year']*10000. + struc['month']*100. + struc['day'] + struc['hour']/24.
    struc['basin'] = np.asarray(empty, dtype=float).astype(int) 
    struc['sub_basin'] = np.asarray(empty, dtype=float).astype(int) 
        
    Atruc = AttrDict( struc )
    
    return Atruc

def create_2d_storm_array ( nt=None ,nstorm=None , dS=None, fld=None  ):

    if fld in ( 'year','month','day','hour', ):
        data=np.zeros( (nstorm, nt) , dtype=int )
    else:
        data=np.zeros( (nstorm, nt) , dtype=float )

    for i in np.arange( nstorm ):
        for t in np.arange( len(dS[i][fld]) ):
            data[i,t] = dS[i][fld][t]

    return data
        



def parse_file(filename):
    storm_list = []
    current_storm = None
    
    with open(filename, 'r') as file:
        for line in file:
            # Remove leading/trailing whitespaces
            line = line.strip()
            
            if line.startswith('start'):
                # If we encounter 'start', we finish the current storm dict and start a new one
                if current_storm is not None:
                    storm_list.append(current_storm)
                
                # Initialize a new storm dictionary
                parts = line.split()
                current_storm = {
                    'lat': [],
                    'lon': [],
                    'slp': [],
                    'wind': [],
                    'surface_pressure': [],
                    'timestamp': []
                }
                continue
            
            # For non-'start' lines, read the data points
            parts = line.split()

            # Extract the values for the storm dict
            lon = float(parts[1])
            lat = float(parts[2])
            slp = float(parts[3])
            wind = float(parts[4])
            surface_pressure = float(parts[5])
            year = int(parts[6])
            month = int(parts[7])
            day = int(parts[8])
            hour = int(parts[9])
            
            # Append to the current storm data
            current_storm['lon'].append(lon)
            current_storm['lat'].append(lat)
            current_storm['slp'].append(slp)
            current_storm['wind'].append(wind)
            current_storm['surface_pressure'].append(surface_pressure)
            current_storm['timestamp'].append(f'{year}-{month:02d}-{day:02d} {hour:02d}:00')
        
        # Append the last storm
        if current_storm is not None:
            storm_list.append(current_storm)
    
    return storm_list

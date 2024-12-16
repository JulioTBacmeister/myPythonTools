#!/usr/bin/env python
################################################
# New style 
################################################
import sys
import os

workdir_ = '/glade/work/juliob'
if ( workdir_ not in sys.path ):
    sys.path.append(workdir_)
    print( f" a path to {workdir_} added in {__name__} {utils_path} ")

from myPythonTools.Utils import utils as uti
from myPythonTools.Utils import numerical_utils as nuti
from myPythonTools.Utils import AveragingUtils as Av
from myPythonTools.Utils import validation_data as Val
from myPythonTools.Utils import PlotUtil as Pu
from myPythonTools.Plotting import LatLonMaps as LL 
from myPythonTools.CASutils import filter_utils as fu 

from PyRegridding.Utils import GridUtils as GrU
from PyRegridding.Utils import MakePressures as MkP
from PyRegridding.Drivers import RegridField as RgF



# The usual
from datetime import date
import numpy as np
import xarray as xr

# Some other useful packages 
import copy
import time
import cftime
import yaml
#from box import Box


def hemi(lat):
    if lat==0 :
        Hemi='Eq'
    elif lat<0:
        Hemi=f'{int(-lat)}S'
    elif lat>0:
        Hemi=f'{int(lat)}N'

    return Hemi

def procField(A, year, month, day, RgOb_ne240_x_llOxO, FLD=None ):


    ##########################################
    # Set limits for Latitude slice and make
    # descriptor char
    ###########################################
    lat_lims = [-65.,-35.]  
    lat_lims_A = f'{hemi(lat_lims[0])}_{hemi(lat_lims[1])}'

    #####################################
    # Make date string etc concatenate
    # into filenmae, make directory name
    ####################################
    dayString = f'{str(year).zfill(4)}-{str(month).zfill(2)}-{str(day).zfill(2)}'
    fname=f'{A.exp}.{A.hsPat}_{FLD}_{lat_lims_A}.{dayString}.nc'
    Bdiro=f'/glade/derecho/scratch/juliob/archive/{A.exp}/atm/GWana/'
    ################################
    # Make final file name and
    # create directory if needed
    ###############################
    fname = f'{Bdiro}{fname}'
    os.makedirs( Bdiro , exist_ok=True )

    print(f' Going to make file={fname}', flush=True)


    #########################################
    # Now start calculations 
    #########################################

    Topo=xr.open_dataset(A.X.topography_file)
    topoO=Topo.PHIS.values / 9.8

    print( 'init topo ', np.shape(topoO) )
    
    ##########################################
    # Regrid 'fld' and calculate residual
    # from 2x2 (ne16)
    ##########################################
    aaO=A.X[FLD].values
    long_name = A.X[FLD].long_name
    units = A.X[FLD].units

    print( f'FLD={FLD} {long_name} {units}', flush=True )
    
    topoOxll=RgF.Horz(xfld_Src=topoO , Src='ne240pg3', Dst='latlonOxO' , RegridObj_In=RgOb_ne240_x_llOxO   ) 
    aaOxll=RgF.Horz(xfld_Src=aaO , Src='ne240pg3', Dst='latlonOxO' , RegridObj_In=RgOb_ne240_x_llOxO   ) 

    print( 'Regridded topo', np.shape(topoOxll) )

    latO,lonO = GrU.latlon( grid='latlonOxO' )

    lat_in_slice = np.where( (latO >= lat_lims[0] )&(latO <= lat_lims[1] ))
    latO = latO[ lat_in_slice[0] ]
    topoOxll = topoOxll[lat_in_slice[0], :]
    aaOxll = aaOxll[:,:, lat_in_slice[0], :]
    
    lev=A.X.lev.values
    print( 'Regridded/Subsliced topo', np.shape(topoOxll) )

    
    dims   = ["lon","lat","time","lev"]
    coords = dict( 
        lon  = ( ["lon"],lonO),
        lat  = ( ["lat"],latO ),
        lev  = ( ["lev"],lev),
        time = ( ["time"], A.X.time.values  ), #pd.to_datetime( pdTime_ERA[itim] ) ),
    )
    
    Xw = xr.Dataset( coords=coords  )
    Dar = xr.DataArray( data=topoOxll, 
                        dims=('lat','lon',),
                        attrs=dict( long_name='Topography',units='m',) ,) 
    Xw['topoO'] = Dar

    
    Dar = xr.DataArray( data=aaOxll, 
                        dims=('time','lev','lat','lon',),
                        attrs=dict( long_name=long_name,units=units ,) ,) 
    Xw[FLD] = Dar
    
    Xw.to_netcdf( fname )
        

def driver():

    with open('configure_gwana_genl.yaml', 'r') as file:
        cfg = yaml.safe_load(file)
    
    
    ######################################################
    
    #x='ne240x2_QxQsst'
    x='A_xz_phs1'  #'ne240x2'
    #x='oldCTL' #'waccmL135'
    exp, subd, Src, Hkey, Dst, useri = cfg[x]['name'] , cfg[x]['subdir'] , cfg[x]['SrcGrid'] , cfg[x]['Hkey'] , cfg[x]['DstGrid'] , cfg[x]['user'] 
    ymdPatLs = cfg[x]['ymdPat']
    hsPat = cfg[x]['hsPat']
    print( exp, subd, Src, Hkey, Dst, useri , flush=True )
    print( ymdPatLs , flush=True )

    #created_RegridObjs = False
    #####################################
    # Initialize regrid-object library
    RgObLib={}
    
    for ymdPat in ymdPatLs:
        #########################################################################
        #  Here we will make Lat slices from h{1,2,3}{i,a} files at a time
        #########################################################################
        
        print( f"\n ################## \n Gonna do - {ymdPat}" , flush=True )
    
        # Split the string and convert year and month to integers
        year, month = map(int, ymdPat.split('-')[:2])

        days_in_month = uti.days_in_month(year, month, check_for_leap_year=False )

        for day in np.arange( start=1, stop=days_in_month+1, step=1 ):
            dayPat = f'{str(year).zfill(4)}-{str(month).zfill(2)}-{str(day).zfill(2)}-*'

            A = uti.MakeDict4Exp( exp=exp , user=useri, subd=subd , 
                                 hsPat=hsPat , ymdPat=dayPat,verbose=True, open_dataset=True )
        
            """
            if (created_RegridObjs == False ):
                ##########################
                # Create regridding objects
                #############################
                #RgOb_ne240_x_ne16 = RgF.Horz( Src='ne240pg3', Dst='ne16pg3' ) 
                #RgOb_ne16_x_ne240 = RgF.Horz(  Src='ne16pg3' , Dst='ne240pg3' , RegridMethod='BILINEAR' ) 
                RgOb_ne240_x_llOxO = RgF.Horz( Src='ne240pg3', Dst='latlonOxO' ) 
                created_RegridObjs = True 
            else:
                print( f"Using existing Regrid Objects ", flush=True )
            """
        
            RgOb_ne240_x_llOxO = GrU.regrid_object_lib(RgOb=RgObLib, src='ne240pg3', dst='latlonOxO' ) 
        
            # Call procField with required arguments
            procField(A, year, month, day, RgOb_ne240_x_llOxO, FLD='OMEGA')
            procField(A, year, month, day, RgOb_ne240_x_llOxO, FLD='U')
            procField(A, year, month, day, RgOb_ne240_x_llOxO, FLD='V')

if __name__ == "__main__":
    ################################################

    driver()
    

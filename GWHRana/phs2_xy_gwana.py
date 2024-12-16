#!/usr/bin/env python
#Set up paths
import sys
import os

workdir_ = '/glade/work/juliob'
if ( workdir_ not in sys.path ):
    sys.path.append(workdir_)
    print( f" a path to {workdir_} added in {__name__} ")

from myPythonTools.Utils import utils as uti
from myPythonTools.Utils import numerical_utils as nuti
from myPythonTools.Utils import AveragingUtils as Av
from myPythonTools.Utils import validation_data as Val
from myPythonTools.CASutils import filter_utils as fu 

from PyRegridding.Utils import GridUtils as GrU
from PyRegridding.Utils import MakePressures as MkP
from PyRegridding.Drivers import RegridField as RgF

# The usual
from datetime import date
import numpy as np
import xarray as xr

# Some other useful packages 
import importlib
import copy
import time
import cftime
import yaml
import glob
#from box import Box #???


importlib.reload( uti )
importlib.reload( nuti )

importlib.reload(Av)

#importlib.reload(vAB)
importlib.reload(MkP)
importlib.reload(GrU)


def xpyp_calcs(out_before_wrt=False):
    
    ##########################################################
    # Read YAML file and open xarray Datasets
    #---------------------------------------------------------

    
    with open('configure_gwana_genl.yaml', 'r') as file:
        cfg = yaml.safe_load(file)
    
    x='A_xy_phs2'
    exp, subd, Src, Hkey, Dst, useri = cfg[x]['name'] , cfg[x]['subdir'] , cfg[x]['SrcGrid'] , cfg[x]['Hkey'] , cfg[x]['DstGrid'] , cfg[x]['user'] 
    ymdPatLs = cfg[x]['ymdPat']
    hsPat = cfg[x]['hsPat']
    print( exp, subd, Src, Hkey, Dst, useri , flush=True )
    print( ymdPatLs , flush=True )

    #created_RegridObjs = False
    #####################################
    # Initialize regrid-object library
    RgObLib={}

    
    #------------------------------------------------------------------------------------
    #    c153_topfix_ne240pg3_FMTHIST_QxQsst_xic_x02_GWana_OMEGA_2004-06-01_2004-06-30.nc
    #-------------------------------------------------------------------------------------
    for ymdPat in ymdPatLs:
            
        Bdiro = f"/glade/derecho/scratch/juliob/archive/{exp}/atm/GWana/"

        parts = ymdPat.split('-')
        year,month =int(parts[0]), int(parts[1])
        firstday,lastday = 1, uti.days_in_month( year, month )
        ymdString = f"{str(year)}-{str(month).zfill(2)}-{str(firstday).zfill(2)}_{str(year)}-{str(month).zfill(2)}-{str(lastday).zfill(2)}"

        Uname=f"{exp}_GWana_U_{ymdString}.nc"
        Vname=f"{exp}_GWana_V_{ymdString}.nc"
        Oname=f"{exp}_GWana_OMEGA_{ymdString}.nc"
        TSname=f"{exp}_GWana_TS_{ymdString}.nc"
        
        
        XU = xr.open_mfdataset( Bdiro+Uname )
        XV = xr.open_mfdataset( Bdiro+Vname )
        XO = xr.open_mfdataset( Bdiro+Oname )
        
        
        # Define the pattern
        pattern = Bdiro+TSname #'dir/file_*'
        
        # Use glob to find all matches
        matching_files = glob.glob(pattern)
        
        if matching_files:
            print("We have TS output")
            XTS = xr.open_mfdataset( Bdiro+TSname )
            TS_output=True
        else:
            print("TS output not availbale")
            TS_output=False

        print( "Got ", Bdiro+Uname, flush=True )

        print( "Got all the data ", flush=True )
        ##########################################################
        # Generate regirdding objetcs to be usedf later on
        #---------------------------------------------------------
        RgOb_ne240_x_ne16 = GrU.regrid_object_lib(RgOb=RgObLib, src='ne240pg3', dst='ne16pg3' )
        RgOb_ne16_x_fv1x1 = GrU.regrid_object_lib(RgOb=RgObLib,src='ne16pg3', dst='fv1x1',RegridMethod='BILINEAR') 

       
        ##########################################################
        # Extract data from xarray Datasets
        # Do some calcs
        # Make axis variables 
        #---------------------------------------------------------

        devW = XO.omegaO.values - XO.omegaOx2xO.values
        devU = XU.uO.values - XU.uOx2xO.values
        
        devUW = devU * devW 
        devWW = devW * devW 
        devUU = devU * devU 
        
        devUW_c2=RgF.Horz(xfld_Src=devUW , Src='ne240pg3', Dst='ne16pg3' , RegridObj_In=  RgOb_ne240_x_ne16  ) 
        devUW_x1=RgF.Horz(xfld_Src=devUW_c2 ,  Src='ne16pg3' , Dst='fv1x1', RegridObj_In=  RgOb_ne16_x_fv1x1  ) 
        print( "Finshed with devUW " , flush=True )

        devWW_c2=RgF.Horz(xfld_Src=devWW , Src='ne240pg3', Dst='ne16pg3' , RegridObj_In=  RgOb_ne240_x_ne16  ) 
        devWW_x1=RgF.Horz(xfld_Src=devWW_c2 ,  Src='ne16pg3' , Dst='fv1x1', RegridObj_In=  RgOb_ne16_x_fv1x1  ) 
        print( "Finshed with devWW " , flush=True )
        
        devUU_c2=RgF.Horz(xfld_Src=devUU , Src='ne240pg3', Dst='ne16pg3' , RegridObj_In=  RgOb_ne240_x_ne16  ) 
        devUU_x1=RgF.Horz(xfld_Src=devUU_c2 ,  Src='ne16pg3' , Dst='fv1x1', RegridObj_In=  RgOb_ne16_x_fv1x1  ) 
        print( "Finshed with devUU " , flush=True )
        
        
        if (TS_output==True):
            devTS    = XTS.tsO.values -XTS.tsOx2xO.values 
            varTS    = devTS * devTS 
            varTS_c2=RgF.Horz(xfld_Src=varTS , Src='ne240pg3', Dst='ne16pg3' , RegridObj_In=  RgOb_ne240_x_ne16  ) 
            varTS_x1 =RgF.Horz(xfld_Src=varTS_c2 ,  Src='ne16pg3' , Dst='fv1x1', RegridObj_In=  RgOb_ne16_x_fv1x1  ) 
            print( "Finshed with varTS ", flush=True )
        
        
        clat2,clon2 = GrU.latlon( grid='ne16pg3' )
        lat1,lon1 = GrU.latlon( grid='fv1x1' )
        lev = XO.lev.values

        if out_before_wrt==True:
            dex = { 'exp':exp , 'lon':lon1, 'lat':lat1, 'lev':lev, 'time':XO.time.values ,
                   'devUW':devUW_x1 ,
                   'devWW':devWW_x1 ,
                   'devUU':devUU_x1 ,
                   'varTS':varTS_x1 
                  }

            return dex

        
        days = np.arange( len(XO.time.values) )/4.
    
        dims   = ["lon","lat","time","lev","ilev"]
        coords = dict( 
            lon  = ( ["lon"],lon1 ),
            lat  = ( ["lat"],lat1 ),
            lev  = ( ["lev"],lev),
            time = ( ["time"],  XO.time.values ) )
                    
        Xout = xr.Dataset( coords=coords  )
        
        Dar = xr.DataArray( data=devUW_x1 , dims=('time','lev','lat','lon',),
                                    attrs=dict( description='up_omegap',units='Steradians',) ,) 
    
        Xout['uPomegaP'] = Dar

        Dar = xr.DataArray( data=devWW_x1 , dims=('time','lev','lat','lon',),
                                    attrs=dict( description='omegap+2',units='Pa+2 s-2',) ,) 
    
        Xout['omegaP2'] = Dar

        Dar = xr.DataArray( data=devUU_x1 , dims=('time','lev','lat','lon',),
                                    attrs=dict( description='up_omegap',units='Steradians',) ,) 
    
        Xout['uP2'] = Dar

        
        if (TS_output==True):
            Dar = xr.DataArray( data=varTS_x1 , dims=('time','lat','lon',),
                                        attrs=dict( description='up_omegap',units='Steradians',) ,) 
        
            Xout['tsP2'] = Dar



        outname=f"{exp}_GWana_xPyP_{ymdString}.nc"

        print( "Gonna write ", Bdiro+outname, flush=True )
        Xout.to_netcdf( Bdiro+outname )



if __name__ == "__main__":
    ################################################

    xpyp_calcs()
    

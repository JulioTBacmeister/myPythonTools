# Import packages 
import sys
sys.path.append('../Plotting/')
""" Now you can imprt modules in ../Plotting"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.tri as tri
import matplotlib.colors as colors
from scipy import interpolate as intr

#Models
from sklearn.ensemble import RandomForestRegressor
from sklearn.tree import DecisionTreeRegressor


#Evaluation
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import cross_val_score

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Some useful packages 
import importlib
import copy
import time

# IO packages
import pickle
from scipy.io import FortranFile


# Other modules in myPythonTools
import ana as a
#import xyp_plot as xyp

# Read in Forest

def Get_a_Forest( tag , **kwargs ):
    
    
    if 'DoPrediction' in kwargs:
        doPred=kwargs['DoPrediction']
    else:
        doPred=False
    
    
    print( "  Do a prediction= ",doPred )
    
    filename = "AB_" + tag + ".dat"

    ff=FortranFile( filename , 'r')

    # Integer data seems to be int64 by default
    # when using scipy.io.FortranFile
    ddA = ff.read_record( '<i8'   )
    ddB = ff.read_record( '<i8'   )

    # In current code A is a float64, due to inheriting
    # double precision from ridge data ... . Inheritance must
    # happen when ridge data and model history output are 
    # concatenated with " np.r_ " above
    qA = ff.read_record( '<f8'   ).reshape( ddA[0], ddA[1] )

    # B is a float32 since it is composed of model history 
    # output only.
    qB = ff.read_record( '<f4'   ).reshape( ddB[0], ddB[1] )
    qA_r = ff.read_record( '<f8'   ).reshape( ddA[0], ddA[1] )
    qB_r = ff.read_record( '<f4'   ).reshape( ddB[0], ddB[1] )
    

    ff.close()

    
    filename = "random_forest_full_"+tag+".pkl"
    MLfile = filename

    tic = time.perf_counter()
    # load model
    Model = pickle.load(open(filename, "rb"))
    toc = time.perf_counter()
    LoadingTime = f"Loaded model {filename} in {toc - tic:0.4f} seconds"
    print(LoadingTime)
    
    
    A_test  = qA_r[788201:,:]
    B_test  = qB_r[788201:,:]
    
    if( doPred == True ):
        tic = time.perf_counter()
        B_pred=Model.predict(A_test)
        toc = time.perf_counter()
        PredictionTime = f"Model prediction in {toc - tic:0.4f} seconds"
        print(PredictionTime)
    else:
        B_pred  = - 999999.
        print("no prediction calculated B_pred=-999999.")


    return A_test, B_test, B_pred

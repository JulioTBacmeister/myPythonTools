import xarray as xr
import netCDF4
import cartopy
import numpy as np
import pandas as pd
import metpy
import warnings
warnings.filterwarnings("ignore")
from tqdm import tqdm
from Tracking_Functions_DF import moaap
data_vars = xr.open_dataset('TEST_CAM7_L58_HR_hourly.200001.nc')

dT = 1 # time interval of input files [hours]
Mask = np.copy(data_vars['X']); Mask[:]=1 # tracking is applied globally
DataName = 'CAM7_L58_HR'
OutputFolder = 'CAM7_moaap_output/'
#time_datetime = pd.to_datetime(np.array(data_vars['time'].values, dtype='datetime64'))
time_datetime = pd.to_datetime(data_vars.indexes['time'].to_datetimeindex())
#Track features using  MOAAP
object_split = moaap(
                      data_vars['X'],
                      data_vars['Y'],
                      time_datetime,
                      dT,
                      Mask,
                      v850 = data_vars['V850'],
                      u850 = data_vars['U850'],
                      t850 = data_vars['T850'],
                      q850 = data_vars['Q850'],
                      slp = data_vars['PSL'],
                      ivte = data_vars['uIVT'],
                      ivtn = data_vars['vIVT'],
                      z500 = data_vars['Z500']*9.81,
                      v200 = data_vars['V200'],
                      u200 = data_vars['U200'],
                      pr   = data_vars['PRECT'],
                      tb   = data_vars['Tb'],
                      DataName = DataName,
                      OutputFolder = OutputFolder,
                      js_min_anomaly = 12,
                      MinTimeJS = 12,
                        )

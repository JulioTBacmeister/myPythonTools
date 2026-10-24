#!/usr/bin/env python
# coding: utf-8
# %%

# %%


# !/usr/bin/env python


# # ScaleDiagram.ipynb

# %%


''' 
   ScaleDiagram.ipynb

   This program reads in netcdf output from previously running 
   papers/2021_Object-Tracker-Intro/programs/Global_ProcessIdentification.ipynb
   
   Reads in the size and duration of all features and plots it on a scale diagram

'''


# %%


import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import glob
import os
from pdb import set_trace as stop
#from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage import gaussian_filter
from scipy.ndimage import median_filter
from scipy.ndimage import label
from matplotlib import cm
from scipy import ndimage
import random
import scipy
import pickle
import datetime
import pandas as pd
import subprocess
import matplotlib.path as mplPath
import sys
from calendar import monthrange

#### speed up interpolation
import scipy.interpolate as spint
import scipy.spatial.qhull as qhull
import numpy as np
import h5py
import xarray as xr
import netCDF4

# from mpl_toolkits.basemap import Basemap, cm
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
from pylab import *
import string
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import shapefile

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import shapely.geometry as sgeom

from tqdm import tqdm

# from Tracking_Functions import ObjectCharacteristics


# %%


# https://stackoverflow.com/questions/13542855/algorithm-to-find-the-minimum-area-rectangle-for-given-points-in-order-to-comput/33619018#33619018
import numpy as np
from scipy.spatial import ConvexHull

def minimum_bounding_rectangle(points):
    """
    Find the smallest bounding rectangle for a set of points.
    Returns a set of points representing the corners of the bounding box.

    :param points: an nx2 matrix of coordinates
    :rval: an nx2 matrix of coordinates
    """
    from scipy.ndimage.interpolation import rotate
    pi2 = np.pi/2.

    # get the convex hull for the points
    hull_points = points[ConvexHull(points).vertices]

    # calculate edge angles
    edges = np.zeros((len(hull_points)-1, 2))
    edges = hull_points[1:] - hull_points[:-1]

    angles = np.zeros((len(edges)))
    angles = np.arctan2(edges[:, 1], edges[:, 0])

    angles = np.abs(np.mod(angles, pi2))
    angles = np.unique(angles)

    # find rotation matrices
    # XXX both work
    rotations = np.vstack([
        np.cos(angles),
        np.cos(angles-pi2),
        np.cos(angles+pi2),
        np.cos(angles)]).T
#     rotations = np.vstack([
#         np.cos(angles),
#         -np.sin(angles),
#         np.sin(angles),
#         np.cos(angles)]).T
    rotations = rotations.reshape((-1, 2, 2))

    # apply rotations to the hull
    rot_points = np.dot(rotations, hull_points.T)

    # find the bounding points
    min_x = np.nanmin(rot_points[:, 0], axis=1)
    max_x = np.nanmax(rot_points[:, 0], axis=1)
    min_y = np.nanmin(rot_points[:, 1], axis=1)
    max_y = np.nanmax(rot_points[:, 1], axis=1)

    # find the box with the best area
    areas = (max_x - min_x) * (max_y - min_y)
    best_idx = np.argmin(areas)

    # return the best box
    x1 = max_x[best_idx]
    x2 = min_x[best_idx]
    y1 = max_y[best_idx]
    y2 = min_y[best_idx]
    r = rotations[best_idx]

    rval = np.zeros((4, 2))
    rval[0] = np.dot([x1, y2], r)
    rval[1] = np.dot([x2, y2], r)
    rval[2] = np.dot([x2, y1], r)
    rval[3] = np.dot([x1, y1], r)

    return rval


# %%


def Density(XX,
           YY):
    FIN=((~np.isnan(XX)*~np.isnan(YY) == 1))
    x = XX[FIN]
    y = YY[FIN]
    kde = stats.gaussian_kde([x, y])

    xx, yy = np.mgrid[min(x):max(x):(max(x)-min(x))/100, min(y):max(y):(max(y)-min(y))/100]
    density = kde(np.c_[xx.flat, yy.flat].T).reshape(xx.shape)
    return xx, yy, density

import math
def rotate(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy


# %%





# ###  USER MODIFY SECTION

# %%


DataFolder = '/glade/derecho/scratch/fudan/CAM7/L58_MOAAP/CAM7_moaap_output/'
sPlotFile = '/glade/derecho/scratch/hannay/CAM7/L58_MOAAP/'

StartDay = datetime.datetime(2000, 1, 1, 0) 
StopDay = datetime.datetime(2005, 12, 31, 23)
Time = pd.date_range(StartDay, end=StopDay, freq='h')
TimeMM = pd.date_range(StartDay, end=StopDay, freq='12m')
from calendar import monthrange
DaysPerMonth = np.array([monthrange(TimeMM[mm].year, TimeMM[mm].month)[1] for mm in range(len(TimeMM))])
Years = np.unique(Time.year)


ARM_Lat = 20 # AR centroids have to be poleward of this latitude
dT = 1  # temporal resolution in hours

# precipitation tracking options
SmoothSigmaP = 1  # Gaussion std for precipitation smoothing
Pthreshold = 2.5    # precipitation threshold [mm/h]
# minimum AR volume
MinTimeAR = 3 # in 3-hours
MinVolAR = MinTimeAR*10*10 # hours x lon x lat grid cells
# minimum cyclone volume
MinTimeCY = 6 # minimum livetime in hours
# minimum frontal zone area
MinArea = 10*10 # lon x lat grid cells
# cloud tracking setup
SmoothSigmaC = 1
Cthreshold = 241
MinTimeC = 3 # in 3-hours
MinVolC = MinTimeC*10*10 # hours x lon x lat grid cells
# AR tracking
IVTtrheshold = 500 # [kg m-1 s-1]
MinTimeIVT = 3


# ### Read in the Grid

# %%


ncid=Dataset(DataFolder+'200101_CAM7_L58_HR_ObjectMasks__dt-1h_MOAAP-masks.nc', mode='r')
Lat=np.squeeze(ncid.variables['lat'][:])
Lon=np.squeeze(ncid.variables['lon'][:])
ncid.close()


# %%


EarthCircum = 40075000 #[m]
dLat = np.copy(Lon); dLat[:] = EarthCircum/(Lat.shape[0]*2)
dLon = np.copy(Lon)
for la in range(Lat.shape[0]):
    dLon[la,:] = EarthCircum/(Lat.shape[0]*2)*np.cos(np.deg2rad(Lat[la,0]))
Area = dLat/1000.*dLon/1000.
Area[Area < 0] = 0

# get average distance betwenn grid cells
Lat_mean=np.mean(Lat)
Lon_mean=np.mean(Lon)
dLat_mean = np.mean(Lat[1:,:]-Lat[:-1,:])
dLon_mean = np.mean(Lon[:,1:]-Lon[:,:-1])

# Gridspacing=haversine(Lon_mean, dLat_mean, Lon_mean+dLon_mean, dLat_mean+dLat_mean)*1000. #  horizontal grid spacing in m
Gridspacing = np.mean(np.append(dLat[:,:,None],dLon[:,:,None], axis=2))


# ### Read the netCDF file

# %%


# SAVE_File = DataFolder+'MCS_monthly-Stats.npz'

# if os.path.isfile(SAVE_File) == False:

TC_size_time = np.zeros((999999, 2)); TC_size_time[:] = np.nan; iTC=0
CY_size_time = np.zeros((999999, 2)); CY_size_time[:] = np.nan; iCY=0
AC_size_time = np.zeros((999999, 2)); AC_size_time[:] = np.nan; iAC=0
MCS_size_time = np.zeros((999999, 2)); MCS_size_time[:] = np.nan; iMCS=0
AR_size_time = np.zeros((999999, 2)); AR_size_time[:] = np.nan; iAR=0
CYZ_size_time = np.zeros((999999, 2)); CYZ_size_time[:] = np.nan; iCYZ=0
JET_size_time = np.zeros((999999, 2)); JET_size_time[:] = np.nan; iJET=0
RW_size_time = np.zeros((999999, 2)); RW_size_time[:] = np.nan; iRW=0
MRG_size_time = np.zeros((999999, 2)); MRG_size_time[:] = np.nan; iMRG=0
Kelvin_size_time = np.zeros((999999, 2)); Kelvin_size_time[:] = np.nan; iKelvin=0
IGW_size_time = np.zeros((999999, 2)); IGW_size_time[:] = np.nan; iIGW=0
EIG0_size_time = np.zeros((999999, 2)); EIG0_size_time[:] = np.nan; iEIG0=0

for mm in tqdm(range(len(TimeMM))):
    YYYY = str(TimeMM[mm].year)
    MM = str(TimeMM[mm].month).zfill(2)
#     print('    load '+YYYY+MM)

    InputFile = YYYY+MM+'_CAM7_L58_HR_ObjectMasks__dt-1h_MOAAP-masks.nc'

    ncid=Dataset(DataFolder+InputFile, mode='r')
    Lat=np.squeeze(ncid.variables['lat'][:])
    Lon=np.squeeze(ncid.variables['lon'][:])
    CY_obj = np.squeeze(ncid.variables['CY_Objects'])
    TC_obj = np.squeeze(ncid.variables['TC_Objects'])
    AC_obj = np.squeeze(ncid.variables['ACY_Objects'])
    MCS_obj = np.squeeze(ncid.variables['MCS_Tb_Objects'])
    AR_obj = np.squeeze(ncid.variables['AR_Objects'])
    JET_obj = np.squeeze(ncid.variables['JET_Objects'])
    CYZ_obj = np.squeeze(ncid.variables['CY_z500_Objects'])
    RW_obj = np.squeeze(ncid.variables['ER_Objects'])
    MRG_obj = np.squeeze(ncid.variables['MRG_Objects'])
    Kelvin_obj = np.squeeze(ncid.variables['Kelvin_Objects'])
    IGW_obj = np.squeeze(ncid.variables['IGW_Objects'])
    EIG0_obj = np.squeeze(ncid.variables['EIG0_Objects'])
    ncid.close()

    # Work on TCs
    rgiTCs = np.unique(TC_obj)
    rgiTCs = np.delete(rgiTCs,0)
    Objects=ndimage.find_objects(TC_obj.astype(int))
    for jj in range(0,len(Objects)):
        if Objects[jj] == None:
            continue
        DATA = TC_obj[Objects[jj]]
        rgiTC = DATA == jj+1
        if (Objects[jj][0].start > 0) | (Objects[jj][0].stop != AC_obj.shape[0]-1) == False:
            continue
        TC_size_time[iTC,0] = np.mean(np.sum(rgiTC * Area[None,Objects[jj][1], Objects[jj][2]],axis=(1,2)))
        TC_size_time[iTC,1] = np.sum(np.sum(rgiTC, axis=(1,2)) > 0)
        iTC = iTC+1
        
        
    # Work on CY
    CY_obj[TC_obj ==1] = 0
#     rgiCYs = np.unique(CY_obj)
#     rgiCYs = np.delete(rgiCYs,0)
    Objects=ndimage.find_objects(CY_obj.astype(int))
#     if Objects[0] == None:
#         Objects = Objects[1:]
    for jj in range(len(Objects)):
        if Objects[jj] == None:
            continue
        DATA = CY_obj[Objects[jj]]
        rgiCY = DATA == jj+1
        if (Objects[jj][0].start > 0) | (Objects[jj][0].stop != AC_obj.shape[0]-1) == False:
            continue
        CY_size_time[iCY,0] = np.mean(np.sum(rgiCY * Area[None,Objects[jj][1], Objects[jj][2]],axis=(1,2)))
        CY_size_time[iCY,1] = np.sum(np.sum(rgiCY, axis=(1,2)) > 0)
        iCY = iCY+1

    # Work on AC
    # AC_obj[TC_obj ==1] = 0
    Objects=ndimage.find_objects(AC_obj.astype(int))
    for jj in range(len(Objects)):
        if Objects[jj] == None:
            continue
        DATA = AC_obj[Objects[jj]]
        rgiAC = DATA == jj+1
        if (Objects[jj][0].start > 0) | (Objects[jj][0].stop != AC_obj.shape[0]-1) == False:
            continue
        AC_size_time[iAC,0] = np.mean(np.sum(rgiAC * Area[None,Objects[jj][1], Objects[jj][2]],axis=(1,2)))
        AC_size_time[iAC,1] = np.sum(np.sum(rgiAC, axis=(1,2)) > 0)
        iAC = iAC+1
        
        
    # Work on MCSs
    MCS_obj[TC_obj ==1] = 0
    Objects=ndimage.find_objects(MCS_obj.astype(int))
    for jj in range(len(Objects)):
        if Objects[jj] == None:
            continue
        DATA = MCS_obj[Objects[jj]]
        rgiMCS = DATA == jj+1
        if (Objects[jj][0].start > 0) | (Objects[jj][0].stop != MCS_obj.shape[0]-1) == False:
            continue
        MCS_size_time[iMCS,0] = np.mean(np.sum(rgiMCS * Area[None,Objects[jj][1], Objects[jj][2]],axis=(1,2)))
        MCS_size_time[iMCS,1] = np.sum(np.sum(rgiMCS, axis=(1,2)) > 0)
        iMCS = iMCS+1

    # Work on ARs
    # AR_obj[TC_obj ==1] = 0
    Objects=ndimage.find_objects(AR_obj.astype(int))
    for jj in range(len(Objects)):
        if Objects[jj] == None:
            continue
        DATA = AR_obj[Objects[jj]]
        rgiAR = DATA == jj+1
        if np.sum(rgiAR) < 20: # there are spurious points in the AR field
            continue
        if (Objects[jj][0].start > 0) | (Objects[jj][0].stop != AR_obj.shape[0]-1) == False:
            continue
        AR_size_time[iAR,0] = np.mean(np.sum(rgiAR * Area[None,Objects[jj][1], Objects[jj][2]],axis=(1,2)))
        AR_size_time[iAR,1] = np.sum(np.sum(rgiAR, axis=(1,2)) > 0)
        iAR = iAR+1
        
    # Work on CY-500hPa
    Objects=ndimage.find_objects(CYZ_obj.astype(int))
    for jj in range(len(Objects)):
        if Objects[jj] == None:
            continue
        DATA = CYZ_obj[Objects[jj]]
        rgiCYZ = DATA == jj+1
        if np.sum(rgiCYZ) < 20: # there are spurious points in the AR field
            continue
        if (Objects[jj][0].start > 0) | (Objects[jj][0].stop != CYZ_obj.shape[0]-1) == False:
            continue
        CYZ_size_time[iCYZ,0] = np.mean(np.sum(rgiCYZ * Area[None,Objects[jj][1], Objects[jj][2]],axis=(1,2)))
        CYZ_size_time[iCYZ,1] = np.sum(np.sum(rgiCYZ, axis=(1,2)) > 0)
        iCYZ = iCYZ+1
        
    # Work on Jet
    Objects=ndimage.find_objects(JET_obj.astype(int))
    for jj in range(len(Objects)):
        if Objects[jj] == None:
            continue
        DATA = JET_obj[Objects[jj]]
        rgiJET = DATA == jj+1
        if np.sum(rgiJET) < 20: # there are spurious points in the AR field
            continue
        if (Objects[jj][0].start > 0) | (Objects[jj][0].stop != JET_obj.shape[0]-1) == False:
            continue
        JET_size_time[iJET,0] = np.mean(np.sum(rgiJET * Area[None,Objects[jj][1], Objects[jj][2]],axis=(1,2)))
        JET_size_time[iJET,1] = np.sum(np.sum(rgiJET, axis=(1,2)) > 0)
        iJET = iJET+1
    
    # Work on MRG waves
    Objects=ndimage.find_objects(MRG_obj.astype(int))
    for jj in range(len(Objects)):
        if Objects[jj] == None:
            continue
        DATA = MRG_obj[Objects[jj]]
        rgiMRG = DATA == jj+1
        if np.sum(rgiMRG) < 20: # there are spurious points in the AR field
            continue
        if (Objects[jj][0].start > 0) | (Objects[jj][0].stop != MRG_obj.shape[0]-1) == False:
            continue
        MRG_size_time[iMRG,0] = np.mean(np.sum(rgiMRG * Area[None,Objects[jj][1], Objects[jj][2]],axis=(1,2)))
        MRG_size_time[iMRG,1] = np.sum(np.sum(rgiMRG, axis=(1,2)) > 0)
        iMRG = iMRG+1

    # Work on Rossby waves
    Objects=ndimage.find_objects(RW_obj.astype(int))
    for jj in range(len(Objects)):
        if Objects[jj] == None:
            continue
        DATA = RW_obj[Objects[jj]]
        rgiRW = DATA == jj+1
        if np.sum(rgiRW) < 20: # there are spurious points in the AR field
            continue
        if (Objects[jj][0].start > 0) | (Objects[jj][0].stop != RW_obj.shape[0]-1) == False:
            continue
        RW_size_time[iRW,0] = np.mean(np.sum(rgiRW * Area[None,Objects[jj][1], Objects[jj][2]],axis=(1,2)))
        RW_size_time[iRW,1] = np.sum(np.sum(rgiRW, axis=(1,2)) > 0)
        iRW = iRW+1
        
        
    # Work on Kelvin waves
    Objects=ndimage.find_objects(Kelvin_obj.astype(int))
    for jj in range(len(Objects)):
        if Objects[jj] == None:
            continue
        DATA = Kelvin_obj[Objects[jj]]
        rgiKelvin = DATA == jj+1
        if np.sum(rgiKelvin) < 20: # there are spurious points in the AR field
            continue
        if (Objects[jj][0].start > 0) | (Objects[jj][0].stop != Kelvin_obj.shape[0]-1) == False:
            continue
        Kelvin_size_time[iKelvin,0] = np.mean(np.sum(rgiKelvin * Area[None,Objects[jj][1], Objects[jj][2]],axis=(1,2)))
        Kelvin_size_time[iKelvin,1] = np.sum(np.sum(rgiKelvin, axis=(1,2)) > 0)
        iKelvin = iKelvin+1
        
    # Work on IGW waves
    Objects=ndimage.find_objects(IGW_obj.astype(int))
    for jj in range(len(Objects)):
        if Objects[jj] == None:
            continue
        DATA = IGW_obj[Objects[jj]]
        rgiIGW = DATA == jj+1
        if np.sum(rgiIGW) < 20: # there are spurious points in the AR field
            continue
        if (Objects[jj][0].start > 0) | (Objects[jj][0].stop != IGW_obj.shape[0]-1) == False:
            continue
        IGW_size_time[iIGW,0] = np.mean(np.sum(rgiIGW * Area[None,Objects[jj][1], Objects[jj][2]],axis=(1,2)))
        IGW_size_time[iIGW,1] = np.sum(np.sum(rgiIGW, axis=(1,2)) > 0)
        iIGW = iIGW+1
        
    # Work on EIG0 waves
    Objects=ndimage.find_objects(EIG0_obj.astype(int))
    for jj in range(len(Objects)):
        if Objects[jj] == None:
            continue
        DATA = EIG0_obj[Objects[jj]]
        rgiEIG0 = DATA == jj+1
        if np.sum(rgiEIG0) < 20: # there are spurious points in the AR field
            continue
        if (Objects[jj][0].start > 0) | (Objects[jj][0].stop != EIG0_obj.shape[0]-1) == False:
            continue
        EIG0_size_time[iEIG0,0] = np.mean(np.sum(rgiEIG0 * Area[None,Objects[jj][1], Objects[jj][2]],axis=(1,2)))
        EIG0_size_time[iEIG0,1] = np.sum(np.sum(rgiEIG0, axis=(1,2)) > 0)
        iEIG0 = iEIG0+1
        
#     #     for la in tqdm(range(MCS_obj.shape[1])):
#     #         for lo in range(MCS_obj.shape[2]):
#     #             MonthlyMCSNr[mm,la,lo] = len(np.unique(MCS_obj[:,la,lo]))-1

#     np.savez(SAVE_File,
#             MonthlyMCSFrequency=MonthlyMCSFrequency,
#             MonthlyMCSNr=MonthlyMCSNr,
#             MonthlyPR=MonthlyPR,
#             MonthlyMCS_PR=MonthlyMCS_PR,
#             ExtrPR_mon=ExtrPR_mon,
#             ExtrPRMCS_mon=ExtrPRMCS_mon)
# else:
#     DATA = np.load(SAVE_File, allow_pickle=True)
#     MonthlyMCSFrequency= DATA['MonthlyMCSFrequency']
#     MonthlyMCSNr= DATA['MonthlyMCSNr']
#     MonthlyPR= DATA['MonthlyPR']
#     MonthlyMCS_PR= DATA['MonthlyMCS_PR']
#     ExtrPR_mon= DATA['ExtrPR_mon']
#     ExtrPRMCS_mon= DATA['ExtrPRMCS_mon']


# %%


# np.savez('scale_diagram.npz',
#         TC_size_time=TC_size_time,
#         CY_size_time=CY_size_time,
#         AC_size_time=AC_size_time,
#         MCS_size_time=MCS_size_time,
#         AR_size_time=AR_size_time,
#         CYZ_size_time=CYZ_size_time,
#         )

# data = np.load('scale_diagram.npz',
#         TC_size_time=TC_size_time,
#         CY_size_time=CY_size_time,
#         AC_size_time=AC_size_time,
#         MCS_size_time=MCS_size_time,
#         AR_size_time=AR_size_time,
#         CYZ_size_time=CYZ_size_time,
#         )

# TC_size_time=TC_size_time
# CY_size_time=CY_size_time
# AC_size_time=AC_size_time
# MCS_size_time=MCS_size_time
# AR_size_time=AR_size_time
# CYZ_size_time=CYZ_size_time


# %%


from sklearn.neighbors import KernelDensity

def kde2D(x, y, bandwidth, xbins=100j, ybins=100j, **kwargs): 
    """Build 2D kernel density estimate (KDE)."""

    # create grid of sample locations (default: 100x100)
    xx, yy = np.mgrid[x.min()*0.7:x.max()*1.3:xbins, 
                      y.min()*0.7:y.max()*1.3:ybins]

    xy_sample = np.vstack([yy.ravel(), xx.ravel()]).T
    xy_train  = np.vstack([y, x]).T

    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(xy_train)

    # score_samples() returns the log-likelihood of the samples
    z = np.exp(kde_skl.score_samples(xy_sample))
    return xx, yy, np.reshape(z, xx.shape)


# %%


import seaborn as sns

fig = plt.figure(figsize=(9,7))
plt.rcParams.update({'font.size': 15})

gs1 = gridspec.GridSpec(1,1)
gs1.update(left=0.08, right=0.9,
   bottom=0.20, top=0.85,
   wspace=0.05, hspace=0.35)
ax = plt.subplot(gs1[0,0])

# object_names = [['cold clouds', '#737373', '-', 2],
#                    ['surface cyclones', 'k', '-', 2],
#                    ['mid-level cyclones', 'k', '--', 2],
#                    ['anticyclones', '#ff7f00', '-', 2],
#                    ['MCS', '#33a02c', '-', 2],
#                    ['moisture streams', 'r', '-', 2],
#                    ['jets', '#6a3d9a', '-', 2],
#                    ['Rossby waves', '#8c510a', '-', 3],
#                    ['mixed Rossby gravity waves', '#bf812d', '-', 1.5],
#                    ['inertia gravity waves', '#dfc27d','-', 3],
#                    ['Kelvin waves', '#abd9e9','-', 1.5],
#                    ['eastward inertia gravity waves', '#4575b4', '-', 3],
#                    ['fronts', '#cab2d6', '-', 2]]
object_names = [['cold clouds', '#737373', '-', 2],
                   ['surface cyclones', 'k', '-', 2],
                   ['mid-level cyclones', 'k', '--', 2],
                   ['anticyclones', '#ff7f00', '-', 2],
                   ['MCS', '#33a02c', '-', 2],
                   ['moisture streams', 'r', '-', 2],
                   ['jets', '#6a3d9a', '-', 2],
                   ['Rossby waves', '#8c510a', '-', 3],
                   ['mixed Rossby gravity waves', '#bf812d', '-', 1.5],
                   ['inertia gravity waves', '#dfc27d','-', 3],
                   ['Kelvin waves', '#abd9e9','-', 1.5],
                   ['eastward inertia gravity waves', '#4575b4', '-', 3],
                   ['fronts', '#cab2d6', '-', 2]]

FEATURE = ['surface cyclones', 'tropical cyclones', 'mid-level cyclones', 'anticyclones', 'MCS', 'ARs', 'jets', 
           'Rossby waves','mixed Rossby gravity waves', 'Kelvin waves', 'inertia gravity waves','eastward inertia gravity waves']
COLORS = ['k','k', 'k','#ff7f00','#33a02c','r', '#6a3d9a', '#8c510a', '#bf812d', '#dfc27d', '#abd9e9', '#4575b4']
LINESTYLE = ['-','-','-','-','-','--','-','-','-','-',':','-',':',]
PERCENTILES_YY = np.zeros((len(FEATURE),7))
PERCENTILES_XX = np.zeros((len(FEATURE),7))
Smooth = 0.4
CoreCont = 0.2

# =======================================
# bounds for CYs
FE=0
NAN = (~np.isnan(CY_size_time[:,1])) & (CY_size_time[:,0] != 0) & (CY_size_time[:,1] != 0)
xxO = CY_size_time[NAN,1].flatten()*dT
yyO = 2*(CY_size_time[NAN,0].flatten()/np.pi)**0.5   #CY_size_time[NAN,0].flatten()**0.5
PERCENTILES_YY[FE,:] = np.percentile(yyO, (0,5,25,50,75,95,100))
PERCENTILES_XX[FE,:] = np.percentile(xxO, (0,5,25,50,75,95,100))
xxLog = np.log(xxO)
yyLog = np.log(yyO)
xxS, yyS, zzS = kde2D(xxLog, yyLog, Smooth)
# plt.contour(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=1, linewidths=2)
# plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=0.1)
# plt.contour(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=1, linewidths=2)
# plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=0.1)
plt.plot([],[], color=COLORS[FE], lw=2, label=FEATURE[FE], linestyle = LINESTYLE[FE])

# =======================================
# bounds for TCs
FE=1
NAN = ~np.isnan(TC_size_time[:,1])
xxO = TC_size_time[NAN,1].flatten()*dT
yyO = 2*(TC_size_time[NAN,0].flatten()/np.pi)**0.5   #TC_size_time[NAN,0].flatten()**0.5
xxLog = np.log(xxO)
yyLog = np.log(yyO)
PERCENTILES_YY[FE,:] = np.percentile(yyO, (0,5,25,50,75,95,100))
PERCENTILES_XX[FE,:] = np.percentile(xxO, (0,5,25,50,75,95,100))
xxT, yyT, zzT = kde2D(xxLog, yyLog, Smooth)
# plt.contour(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=1, linewidths=2)
# plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=0.1)
# plt.contour(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=1, linewidths=2)
# plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=0.1)
# plt.scatter(TC_size_time[:,1]*dT,TC_size_time[:,0], c='r')
plt.plot([],[], color=COLORS[FE], lw=2, label=FEATURE[FE], linestyle = LINESTYLE[FE])

#print(yyLog)
#print(yyO)
# =======================================
FE=2
# bounds for CY-500hPa
NAN = ~np.isnan(CYZ_size_time[:,1])
xxO = CYZ_size_time[NAN,1].flatten()*dT
yyO = 2*(CYZ_size_time[NAN,0].flatten()/np.pi)**0.5   #CYZ_size_time[NAN,0].flatten()**0.5
PERCENTILES_YY[FE,:] = np.percentile(yyO, (0,5,25,50,75,95,100))
PERCENTILES_XX[FE,:] = np.percentile(xxO, (0,5,25,50,75,95,100))

print(min(yyO))
yyO[yyO<1]+=2
print(min(yyO))
#xxO = 1
#yyO = 1

xxLog = np.log(xxO)
yyLog = np.log(yyO)
xxZ, yyZ, zzZ = kde2D(xxLog, yyLog, Smooth)
# plt.contour(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=1, linewidths=2, linestyle = LINESTYLE[FE])
# plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=0.1)
# plt.contour(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=1, linewidths=2, linestyle = LINESTYLE[FE])
# plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=0.1)
plt.plot([],[], color=COLORS[FE], lw=2, label=FEATURE[FE], linestyle = LINESTYLE[FE])

ax.set_xscale('log')
ax.set_yscale('log')


# =======================================
### PLOT COMBINED CYCLONE CONTOUR
FE=0
xx = (xxS + xxT + xxZ)/3
yy = (yyS + yyT + yyZ)/3
zz = (zzS + zzT + zzZ)/3
plt.contour(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=1, linewidths=2)
plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=0.1)
#plt.contour(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=1, linewidths=2)
plt.contour(np.exp(yy), np.exp(xx), zz, levels=[0.15,1], colors=COLORS[FE], alpha=1, linewidths=2)
plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=0.1)
plt.plot([],[], color=COLORS[FE], lw=2, label=FEATURE[FE], linestyle = LINESTYLE[FE])

# =======================================
# bounds for ACs
FE=3
NAN = (~np.isnan(AC_size_time[:,1])) & (AC_size_time[:,0] != 0) & (AC_size_time[:,1] != 0)
xxO = AC_size_time[NAN,1].flatten()*dT
yyO = 2*(AC_size_time[NAN,0].flatten()/np.pi)**0.5   # AC_size_time[NAN,0].flatten()**0.5
PERCENTILES_YY[FE,:] = np.percentile(yyO, (0,5,25,50,75,95,100))
PERCENTILES_XX[FE,:] = np.percentile(xxO, (0,5,25,50,75,95,100))
xxLog = np.log(xxO)
yyLog = np.log(yyO)
xx, yy, zz = kde2D(xxLog, yyLog, Smooth)
plt.contour(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=1, linewidths=2)
plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=0.1)
plt.contour(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=1, linewidths=2)
plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=0.1)
# plt.scatter(AC_size_time[:,0],AC_size_time[:,1]*dT, c='r')
plt.plot([],[], color=COLORS[FE], lw=2, label=FEATURE[FE], linestyle = LINESTYLE[FE])

# =======================================
FE=4
# bounds for MCSs
# NAN = (~np.isnan(MCS_size_time[:,1])) & (MCS_size_time[:,0] != 0) & (MCS_size_time[:,1] != 0)
NAN = ~np.isnan(MCS_size_time[:,1])
xxO = MCS_size_time[NAN,1].flatten()*dT
yyO = 2*(MCS_size_time[NAN,0].flatten()/np.pi)**0.5 # MCS_size_time[NAN,0].flatten()**0.5
PERCENTILES_YY[FE,:] = np.percentile(yyO, (0,5,25,50,75,95,100))
PERCENTILES_XX[FE,:] = np.percentile(xxO, (0,5,25,50,75,95,100))
xxLog = np.log(xxO)
yyLog = np.log(yyO)
xx, yy, zz = kde2D(xxLog, yyLog, Smooth)
plt.contour(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=1, linewidths=2)
plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=0.1)
plt.contour(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=1, linewidths=2)
plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=0.1)
# plt.scatter(MCS_size_time[:,0],MCS_size_time[:,1]*dT, c='r', s=3)
plt.plot([],[], color=COLORS[FE], lw=2, label=FEATURE[FE], linestyle = LINESTYLE[FE])

# =======================================
FE=5
# bounds for ARs
NAN = ~np.isnan(AR_size_time[:,1])
xxO = AR_size_time[NAN,1].flatten()*dT
yyO = 2*(AR_size_time[NAN,0].flatten()/np.pi)**0.5   #AR_size_time[NAN,0].flatten()**0.5
PERCENTILES_YY[FE,:] = np.percentile(yyO, (0,5,25,50,75,95,100))
PERCENTILES_XX[FE,:] = np.percentile(xxO, (0,5,25,50,75,95,100))
xxLog = np.log(xxO)
yyLog = np.log(yyO)
xx, yy, zz = kde2D(xxLog, yyLog, Smooth)
plt.contour(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=1, linewidths=2)
plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=0.1)
plt.contour(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=1, linewidths=2)
plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=0.1)
# plt.scatter(AR_size_time[:,0],AR_size_time[:,1]*dT, c='#33a02c', s=3)
plt.plot([],[], color=COLORS[FE], lw=2, label=FEATURE[FE], linestyle = LINESTYLE[FE])

ax.set_xscale('log')
ax.set_yscale('log')



# =======================================
FE=6
# bounds for Jets
NAN = ~np.isnan(JET_size_time[:,1])
xxO = JET_size_time[NAN,1].flatten()*dT
yyO = 2*(JET_size_time[NAN,0].flatten()/np.pi)**0.5   #JET_size_time[NAN,0].flatten()**0.5
PERCENTILES_YY[FE,:] = np.percentile(yyO, (0,5,25,50,75,95,100))
PERCENTILES_XX[FE,:] = np.percentile(xxO, (0,5,25,50,75,95,100))
xxLog = np.log(xxO)
yyO[yyO<1]+=2
yyLog = np.log(yyO)
xx, yy, zz = kde2D(xxLog, yyLog, Smooth)
plt.contour(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=1, linewidths=2)
plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=0.1)
plt.contour(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=1, linewidths=2)
plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=0.1)
# plt.scatter(JET_size_time[:,0],JET_size_time[:,1]*dT, c='#33a02c', s=3)
plt.plot([],[], color=COLORS[FE], lw=2, label=FEATURE[FE], linestyle = LINESTYLE[FE])

ax.set_xscale('log')
ax.set_yscale('log')


# =======================================
FE=7
# bounds for RW
NAN = ~np.isnan(RW_size_time[:,1])
xxO = RW_size_time[NAN,1].flatten()*dT
yyO = 2*(RW_size_time[NAN,0].flatten()/np.pi)**0.5   #EIG0_size_time[NAN,0].flatten()**0.5
PERCENTILES_YY[FE,:] = np.percentile(yyO, (0,5,25,50,75,95,100))
PERCENTILES_XX[FE,:] = np.percentile(xxO, (0,5,25,50,75,95,100))
xxLog = np.log(xxO)
yyO[yyO<1]+=2
yyLog = np.log(yyO)
xxR, yyR, zzR = kde2D(xxLog, yyLog, Smooth)
# plt.contour(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=1, linewidths=2)
# plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=0.1)
# plt.contour(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=1, linewidths=2)
# plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=0.1)
# plt.scatter(EIG0_size_time[:,0],EIG0_size_time[:,1]*dT, c='#33a02c', s=3)
plt.plot([],[], color=COLORS[FE], lw=2, label=FEATURE[FE], linestyle = LINESTYLE[FE])

# =======================================
FE=8
# bounds for MRG
NAN = ~np.isnan(MRG_size_time[:,1])
xxO = MRG_size_time[NAN,1].flatten()*dT
yyO = 2*(MRG_size_time[NAN,0].flatten()/np.pi)**0.5   #MRG_size_time[NAN,0].flatten()**0.5
PERCENTILES_YY[FE,:] = np.percentile(yyO, (0,5,25,50,75,95,100))
PERCENTILES_XX[FE,:] = np.percentile(xxO, (0,5,25,50,75,95,100))
xxLog = np.log(xxO)
yyO[yyO<1]+=2
yyLog = np.log(yyO)
xxM, yyM, zzM = kde2D(xxLog, yyLog, Smooth)
# plt.contour(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=1, linewidths=2)
# plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=0.1)
# plt.contour(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=1, linewidths=2)
# plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=0.1)
# plt.scatter(MRG_size_time[:,0],MRG_size_time[:,1]*dT, c='#33a02c', s=3)
plt.plot([],[], color=COLORS[FE], lw=2, label=FEATURE[FE], linestyle = LINESTYLE[FE])

ax.set_xscale('log')
ax.set_yscale('log')


# =======================================
FE=9
# bounds for Kelvin
NAN = ~np.isnan(Kelvin_size_time[:,1])
xxO = Kelvin_size_time[NAN,1].flatten()*dT
yyO = 2*(Kelvin_size_time[NAN,0].flatten()/np.pi)**0.5   #Kelvin_size_time[NAN,0].flatten()**0.5
PERCENTILES_YY[FE,:] = np.percentile(yyO, (0,5,25,50,75,95,100))
PERCENTILES_XX[FE,:] = np.percentile(xxO, (0,5,25,50,75,95,100))
xxLog = np.log(xxO)
yyO[yyO<1]+=2
yyLog = np.log(yyO)
xxK, yyK, zzK = kde2D(xxLog, yyLog, Smooth)
# plt.contour(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=1, linewidths=2)
# plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=0.1)
# plt.contour(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=1, linewidths=2)
# plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=0.1)
# plt.scatter(Kelvin_size_time[:,0],Kelvin_size_time[:,1]*dT, c='#33a02c', s=3)
plt.plot([],[], color=COLORS[FE], lw=2, label=FEATURE[FE], linestyle = LINESTYLE[FE])


# =======================================
FE=10
# bounds for IGW
NAN = ~np.isnan(IGW_size_time[:,1])
xxO = IGW_size_time[NAN,1].flatten()*dT
yyO = 2*(IGW_size_time[NAN,0].flatten()/np.pi)**0.5   #IGW_size_time[NAN,0].flatten()**0.5
PERCENTILES_YY[FE,:] = np.percentile(yyO, (0,5,25,50,75,95,100))
PERCENTILES_XX[FE,:] = np.percentile(xxO, (0,5,25,50,75,95,100))
xxLog = np.log(xxO)
yyO[yyO<1]+=2
yyLog = np.log(yyO)
xxI, yyI, zzI = kde2D(xxLog, yyLog, Smooth)
# plt.contour(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=1, linewidths=2)
# plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=0.1)
# plt.contour(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=1, linewidths=2)
# plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=0.1)
# plt.scatter(IGW_size_time[:,0],IGW_size_time[:,1]*dT, c='#33a02c', s=3)
plt.plot([],[], color=COLORS[FE], lw=2, label=FEATURE[FE], linestyle = LINESTYLE[FE])


# =======================================
FE=11
# bounds for EIG0
NAN = ~np.isnan(EIG0_size_time[:,1])
xxO = EIG0_size_time[NAN,1].flatten()*dT
yyO = 2*(EIG0_size_time[NAN,0].flatten()/np.pi)**0.5   #EIG0_size_time[NAN,0].flatten()**0.5
PERCENTILES_YY[FE,:] = np.percentile(yyO, (0,5,25,50,75,95,100))
PERCENTILES_XX[FE,:] = np.percentile(xxO, (0,5,25,50,75,95,100))
xxLog = np.log(xxO)
yyO[yyO<1]+=2
yyLog = np.log(yyO)
xxE, yyE, zzE = kde2D(xxLog, yyLog, Smooth)
# plt.contour(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=1, linewidths=2)
# plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=0.1)
# plt.contour(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=1, linewidths=2)
# plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=0.1)
# plt.scatter(EIG0_size_time[:,0],EIG0_size_time[:,1]*dT, c='#33a02c', s=3)
plt.plot([],[], color=COLORS[FE], lw=2, label=FEATURE[FE], linestyle = LINESTYLE[FE])



# # =======================================
# # PLOT EASTWARD COMBINED TROPICAL WAVES
# FE=7
# xx = (xxE + xxK + xxI + xxM)/4
# yy = (yyE + yyK + yyI + yyM)/4
# zz = (zzE + zzK + zzI + zzM)/4
# plt.contour(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=1, linewidths=2)
# plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=0.1)
# plt.contour(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=1, linewidths=2)
# plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=0.1)


# ax.set_xscale('log')
# ax.set_yscale('log')

# =======================================
# PLOT EASTWARD COMBINED TROPICAL WAVES
FE=7
xx = (xxR + xxM)/2
yy = (yyR + yyM)/2
zz = (zzR + zzM)/2
plt.contour(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=1, linewidths=2)
plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=0.1)
plt.contour(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=1, linewidths=2)
plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=0.1)

ax.set_xscale('log')
ax.set_yscale('log')

# =======================================
# PLOT WESTWARD COMBINED TROPICAL WAVES
FE=11
xx = (xxE + xxK + xxI)/3
yy = (yyE + yyK + yyI)/3
zz = (zzE + zzK + zzI)/3
plt.contour(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=1, linewidths=2)
plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[0.10,1], colors=COLORS[FE], alpha=0.1)
plt.contour(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=1, linewidths=2)
plt.contourf(np.exp(yy), np.exp(xx), zz, levels=[CoreCont,1], colors=COLORS[FE], alpha=0.1)


ax.set_xscale('log')
ax.set_yscale('log')

    
# plt.legend(loc='center left', bbox_to_anchor=(1.3, 0.5))

# https://www.researchgate.net/publication/331404250_Three-Dimensional_Observation_of_Atmospheric_Processes_in_Cities/figures?lo=1
# plot sound wave speed
plt.plot([1188,1188*24*100],[1,24*100], c='k')
plt.text(4000,2.8, 'sound waves', rotation=32)

plt.plot([(10**1)/1000,(10**6)/1000],[(10**2)/3600,(10**7)/3600], c='k')
plt.text(80,5*24, 'wind speed detection limit', rotation=32, ha='center')

ax.set_ylabel('Time [hours]')
# ax.set_xticks(range(3,1000))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.ylim(2,1000)
ax.set_xlabel('horizontal length scale [km]')
plt.xlim(20,10000)

# add x labels
Ylabels = ['6 hours', '1 day', '10 days']
YlabelLoc = [6,24,240]
ax.set_yticks(YlabelLoc)
ax.set_yticklabels(Ylabels)

for yy in range(len(Ylabels)):
    plt.axhline(y=YlabelLoc[yy], color='k', linestyle='-', lw=0.5)
#     tt = plt.text(4300,YlabelLoc[yy], Ylabels[yy],ha='left',va='center',weight="bold")
#     tt.set_bbox(dict(facecolor='w', alpha=1, edgecolor='w'))


### Add labels for length scale
pos1 = ax.get_position()
ax2 = axes([pos1.x0,0.01,pos1.width, 0.08])

# https://www.sciencedirect.com/topics/earth-and-planetary-sciences/meteorological-phenomenon
rgsSpatScales = ['Meso\ngamma','Meso\nbeta','Meso\nalpha','Macro\nbeta','Macro\nalpha']
rgrSpatScales = [0.2, 2,20,200,2000,10000]
Ceter_SpatScales = [0.5,5,80,800,5000]
for ss in range(len(rgsSpatScales))[2:]:
    plt.text(Ceter_SpatScales[ss],0.5, rgsSpatScales[ss],ha='right', va='center')
    plt.axvline(x=rgrSpatScales[ss], c='k', lw=2)
    ax.axvline(x=rgrSpatScales[ss], c='k', lw=1)
plt.xlim(20,10000)
ax2.set_xscale('log')
plt.axis('off')



# ====================
# Add box-whisker plots for horizongal scale
ax3 = axes([pos1.x0,0.86,pos1.width, 0.3])
plt.xlim(20,10000)
ax3.set_xscale('log')

for FE in range(len(FEATURE)):
    plt.plot([PERCENTILES_YY[FE,2], PERCENTILES_YY[FE,4]],[FE,FE], lw=8, color=COLORS[FE], solid_capstyle='butt')
    plt.plot([PERCENTILES_YY[FE,1], PERCENTILES_YY[FE,5]],[FE,FE], lw=2, color=COLORS[FE], solid_capstyle='butt')
    plt.plot(PERCENTILES_YY[FE,3],FE,'wo', markersize=3)
    print(' '+FEATURE[FE]+' spatial scale ->  5th prct: '+str(PERCENTILES_YY[FE,1])+'; 25th prct: '+str(PERCENTILES_YY[FE,2])+'; 50th prct: '+str(PERCENTILES_YY[FE,3])+'; 75th prct: '+str(PERCENTILES_YY[FE,4])+'; 95th prct: '+str(PERCENTILES_YY[FE,5]))
    plt.text(20 ,FE, FEATURE[FE], 
             color=COLORS[FE], fontsize = 11, ha='left', va = 'center')
    # plt.plot(PERCENTILES_YY[FE,0],FE,marker='o', color = COLORS[FE], markersize=3)
    # plt.plot(PERCENTILES_YY[FE,6],FE,marker='o', color = COLORS[FE], markersize=3)
plt.ylim(-1,len(FEATURE))
plt.axis('off')

# ====================
# Add box-whisker plots for vertical scale
ax3 = axes([pos1.x0+pos1.width+0.01,pos1.y0,0.2, pos1.height])
plt.ylim(2,1000)
ax3.set_yscale('log')

for FE in range(len(FEATURE)):
    plt.plot([FE,FE],[PERCENTILES_XX[FE,2], PERCENTILES_XX[FE,4]], lw=8,color=COLORS[FE], solid_capstyle='butt')
    print(' '+FEATURE[FE]+' temporal scale ->  5th prct: '+str(PERCENTILES_XX[FE,1])+'; 25th prct: '+str(PERCENTILES_XX[FE,2])+'; 50th prct: '+str(PERCENTILES_XX[FE,3])+'; 75th prct: '+str(PERCENTILES_XX[FE,4])+'; 95th prct: '+str(PERCENTILES_XX[FE,5]))
    #plt.plot([FE,FE],[PERCENTILES_XX[FE,1], PERCENTILES_XX[FE,5]], lw=2,color=COLORS[FE], solid_capstyle='butt')
    #print(' '+FEATURE[FE]+' temporal  5th prct: '+str(PERCENTILES_XX[FE,1])+';  95th prct: '+str(PERCENTILES_XX[FE,5]))
    #print(' '+FEATURE[FE]+' temporal   50th prct: '+str(PERCENTILES_XX[FE,3]))
    plt.plot(FE,PERCENTILES_XX[FE,3],'wo', markersize=3)
    plt.text(FE, 40*24 , FEATURE[FE], 
             color=COLORS[FE], fontsize = 11, ha='center', va = 'bottom', rotation = 90)
    # plt.plot(FE,PERCENTILES_XX[FE,0],marker='o', color = COLORS[FE], markersize=3)
    # plt.plot(FE,PERCENTILES_XX[FE,6],marker='o', color = COLORS[FE], markersize=3)
plt.xlim(-1,len(FEATURE))
plt.axis('off')
    
sPlotFile=sPlotFile
# sPlotName= 'BottomUp-'+str(rgrClustersFin[1].max()+1)+'WT_Centroids.pdf'
sPlotFile=''
sPlotName= 'Atmospheric-Feature_Scale-Diagram_CAM7_HR_L58.pdf'
if os.path.isdir(sPlotFile) != 1:
    subprocess.call(["mkdir","-p",sPlotFile])
print('        Plot map to: '+sPlotFile+sPlotName)
fig.savefig(sPlotFile+sPlotName, bbox_inches='tight')


# %%





# %%





# %%





# Import packages 
import sys
sys.path.append('../Plotting/')
""" Now you can imprt modules in ../Plotting"""
import xyp_plot as xyp

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.tri as tri
import matplotlib.colors as colors

from scipy.io import FortranFile

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
#import ana as a

import ESMF as E

import importlib

import copy
import glob


def Regrid( srcScrip , dstScrip , srcType , dstType ,  RegridMethod="CONSERVE" , **kwargs ):

    if(RegridMethod.upper()=='CONSERVE'):
        regrid_method=E.RegridMethod.CONSERVE 
    if(RegridMethod.upper()=='BILINEAR'):
        regrid_method=E.RegridMethod.BILINEAR 
    
    if(srcType.lower()=='mesh'):
        srcDesc=E.Mesh( filename=srcScrip,
            filetype=E.FileFormat.SCRIP )
        
    if(srcType.lower()=='grid'):
        srcDesc=E.Grid( filename=srcScrip ,
            filetype=E.FileFormat.SCRIP ,
            add_corner_stagger=True   )
        
    if(dstType.lower()=='mesh'):
        dstDesc=E.Mesh( filename=dstScrip,
            filetype=E.FileFormat.SCRIP )
        
    if(dstType.lower()=='grid'):
        dstDesc=E.Grid( filename=dstScrip ,
            filetype=E.FileFormat.SCRIP ,
            add_corner_stagger=True   )
        
        
    if(srcType.lower()=='mesh'):
        srcField = E.Field(srcDesc, meshloc=E.MeshLoc.ELEMENT )
    if(srcType.lower()=='grid'):
        srcField = E.Field(srcDesc )

    if(dstType.lower()=='mesh'):
        dstField = E.Field(dstDesc, meshloc=E.MeshLoc.ELEMENT )
    if(dstType.lower()=='grid'):
        dstField = E.Field(dstDesc )

    
    dstField.data[:]=1e20
    srcField.data[:]=1e20


    Regrd = E.Regrid( srcField , dstField , 
                      regrid_method=regrid_method,
                      unmapped_action=E.UnmappedAction.IGNORE)


    return Regrd, srcField , dstField 

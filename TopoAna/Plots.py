import os
import sys

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gsp
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import importlib
import glob
import copy

"""
Some routines for plotting. Mostly topo.
"""

def sixpanel(aa,clev=21,alpha=1.0, goofy=False, cmap='viridis' ):
    fig=plt.figure( figsize=(21,15))
    gs = gsp.GridSpec( 3, 4, figure=fig ,wspace=0.0, hspace=0.0 )
    
    ax0=fig.add_subplot(gs[1, 0])
    ax1=fig.add_subplot(gs[1, 1])
    ax2=fig.add_subplot(gs[1, 2])
    ax3=fig.add_subplot(gs[1, 3])
    if (goofy==False):
        ax4=fig.add_subplot(gs[2, 0])
    else:
        ax4=fig.add_subplot(gs[2, 3])
    ax5=fig.add_subplot(gs[0, 0])
    for i, ax in enumerate(fig.axes):
        if ( i==4 and goofy==True ):
            aax = np.flip( np.transpose( aa[i,:,:] ) , axis=0 )
        else:
            aax = aa[ i, :,:]
        ax.contourf(  aax , levels=clev , alpha=alpha , cmap=cmap )

    
    
def format_axes(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        ax.tick_params(labelbottom=False, labelleft=False)

def CAMridgelet( lat, lon, angll, hwdth=0., clngt=100. ):
    
    RR=6371.0
    rlat, rlon, angr = lat*np.pi/180. , lon*np.pi/180. , angll*np.pi/180.
    x0 = RR * np.cos( rlat ) * rlon
    y0 = RR * rlat
    dx = 0.5*clngt * np.sin( angr )
    dy = 0.5*clngt * np.cos( angr )
    x1 = x0 + dx
    y1 = y0 + dy
    x2 = x0 - dx
    y2 = y0 - dy
    
    lon1 = (180./np.pi) * x1/(RR*np.cos(rlat) )
    lat1 = (180./np.pi) * y1/RR
    lon2 = (180./np.pi) * x2/(RR*np.cos(rlat) )
    lat2 = (180./np.pi) * y2/RR
    
    return lat1,lon1,lat2,lon2

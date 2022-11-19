#!/usr/bin/env python
import sys
sys.path.append('../Plotting/')

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import argparse as arg

import xyp_plot as xyp
from datawork import camdata as cm

def test_xp_sxn():

    f1='c6_3_59.ne30pg3_L32_latlon025_SAMwrf.ndg05.cam.h1.2010-06-15-00000.nc'

    a=xr.open_dataset(f1)

    U=a['U']
    lon=a['lon']
    plev=a['lev']

    """
    hyam=a['hyam']
    hybm=a['hybm']
    PS=a['PS']
    plev=a['lev']
    nlev=plev.size

    P3=U*0.

    for L in np.arange( nlev  ):
        P3[0,L,:,:]=hyam[L]*100000. + hybm[L]*PS[0,:,:]
    """
    cc=cm()
    P3 = cc.make_pmid(a)

    ulv=np.linspace(-40,40,num=31)


    plt.ion()

    xyp.pltxp(a=U,p=P3,x=lon,j0=200,plev=plev,xlim=[270,340],ylim=[100000.,1000.],clevs=ulv)

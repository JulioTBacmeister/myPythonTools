#!/usr/bin/env python

import hfdata as h
import xarray as xr
import numpy as np
import getopt as go
import subprocess as sp
import sys
import os

#----------------------------------
# Run from python prompt like this:
#  exec(open("./testh0.py").read())


argv=sys.argv
lonfill_=False
point_=False

try:
   opts, args = go.getopt( argv[1:], "y:m:X:", 
                           ["year=","month=","case="] )
except:
    print( "something is wrong")
    exit()

for opt, arg in opts:
    if opt in ("-y","--year"):
        year = int(arg)
    elif opt in ("-m","--month"):
        month = int(arg)
    elif opt in ("-X","--case"):
        xp = arg

idir='/project/amp/juliob/CAM/'+xp+'/f09_omega/L58/2010/'
odir='/project/amp/juliob/CAM/'+xp+'/f09_omega_phis/L58/2010/'

# create output directory
sp.run("mkdir -p "+ odir,shell=True)

days_in_month =[31 , 28, 31, 30, 31, 30, 31, 31, 30, 31,30, 31 ]

f="/fs/cgd/csm/inputdata/atm/cam/topo/fv_0.9x1.25_nc3000_Nsw042_Nrs008_Co060_Fi001_ZR_sgh30_24km_GRNL_c170103.nc"
topo=xr.open_dataset( f )
phis=topo['PHIS']

x = h.hfdata(xp=xp,dir=idir)
y = h.hfdata(xp=xp,dir=odir)

lfilo=[]
m=month-1
nd = days_in_month[m]
for d in np.arange(nd):
    dd=d+1
    for h in np.arange(start=0,stop=23,step=6):
        ss=h*3600
        fili=x.filename( year=year,month=month,day=dd,second=ss,moniker='cam.h1')
        #print(fili)
        filo=y.filename( year=year,month=month,day=dd,second=ss,moniker='cam.h1')
        y.add_phis( fili, filo, phis )
        lfilo.append( filo ) 
      

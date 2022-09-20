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

#xp='c6_3_59.f09_L58.CTL01'
#year=2010
#month=1

argv=sys.argv
lonfill_=False
point_=False

try:
   opts, args = go.getopt( argv[1:], "y:m:X:FP", 
                           ["year=","month=","case=","LonFill","Point"] )
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
    elif opt in ("-F","--LonFill"):
        lonfill_=True
    elif opt in ("-P","--Point"):
        point_=True

if (point_ == True ):
   idir='/project/amp/juliob/CAM/'+xp+'/f09_omega/L58/2010/'
   odir='/project/amp/juliob/CAM/'+xp+'/f09_omega_point/L58/2010/'
elif (lonfill_ == True ):
   idir='/project/amp/juliob/CAM/'+xp+'/f09_omega/L58/2010/'
   odir='/project/amp/juliob/CAM/'+xp+'/f09_omega_zonav/L58/2010/'
elif (lonfill_ == False ):
   idir='/project/amp/juliob/CAM/'+xp+'/f09_omega/L58/2010/'
   odir='/project/amp/juliob/CAM/'+xp+'/f09_omega_yz/L58/2010/'

# create output directory
sp.run("mkdir -p "+ odir,shell=True)

days_in_month =[31 , 28, 31, 30, 31, 30, 31, 31, 30, 31,30, 31 ]


x = h.hfdata(xp=xp,dir=idir)
y = h.hfdata(xp=xp,dir=odir)

lfilo=[]
for m in np.arange(12):
    mm=m+1
    dd=1
    ss=0
    if (point_ == True):
       filo=y.filename( year=year,month=mm,day=dd,second=ss,moniker='cam_point.h1')
    elif (lonfill_ == True):
       filo=y.filename( year=year,month=mm,day=dd,second=ss,moniker='cam_zonav.h1')
    elif (lonfill_ == False):
       filo=y.filename( year=year,month=mm,day=dd,second=ss,moniker='cam_yz.h1')
       
    lfilo.append( y.monthly_file ) 
      
print(lfilo)

ds=xr.open_mfdataset( lfilo )
ds.to_netcdf( y.yearly_file )

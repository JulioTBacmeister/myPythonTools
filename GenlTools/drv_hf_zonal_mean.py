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
m=month-1
nd = days_in_month[m]
for d in np.arange(nd):
    dd=d+1
    for h in np.arange(start=0,stop=23,step=6):
        ss=h*3600
        fili=x.filename( year=year,month=month,day=dd,second=ss,moniker='cam.h1')
        print(fili)

        if (point_ == True):
           filo=y.filename( year=year,month=month,day=dd,second=ss,moniker='cam_point.h1')
           y.point_output( fili, filo )
        elif (lonfill_ == True):
           filo=y.filename( year=year,month=month,day=dd,second=ss,moniker='cam_zonav.h1')
           y.zonal_mean_lonfill( fili, filo )
        elif (lonfill_ == False):
           filo=y.filename( year=year,month=month,day=dd,second=ss,moniker='cam_yz.h1')
           y.zonal_mean( fili, filo )

        lfilo.append( filo ) 
      
print(lfilo)
print(y.monthly_file)
print(y.yearly_file)

if (lonfill_ == False):
   ds=xr.open_mfdataset( lfilo )
   ds.to_netcdf( y.monthly_file )

#!/usr/bin/env python

import hfdata as h
import numpy as np
import getopt as go
import sys
import os

#----------------------------------
# Run from python prompt like this:
#  exec(open("./testh0.py").read())

idir='/project/amp/juliob/CAM/f09_omega/L58/2010/'
odir='/project/amp/juliob/CAM/f09_omega_zonav/L58/2010/'
xp='c6_3_59.f09_L58.ndg01'

year=2010
month=1
argv=sys.argv

try:
   opts, args = go.getopt( argv[1:], "y:m:", 
                           ["year=","month="] )
except:
    print( "something is wrong")
    exit()



for opt, arg in opts:
    if opt in ("-y","--year"):
        year = int(arg)
    elif opt in ("-m","--month"):
        month = int(arg)


days_in_month =[31 , 28, 31, 30, 31, 30, 31, 31, 30, 31,30, 31 ]


x = h.hfdata(xp=xp,dir=idir)
y = h.hfdata(xp=xp,dir=odir)

m=month-1
nd = days_in_month[m]
for d in np.arange(nd):
    dd=d+1
    for h in np.arange(start=0,stop=23,step=6):
        ss=h*3600
        fili=x.filename( year=year,month=month,day=dd,second=ss,moniker='cam.h1')
        filo=y.filename( year=year,month=month,day=dd,second=ss,moniker='cam_zonav.h1')
        #print(fili)
        print(filo)
        
        y.zonal_mean_lonfill( fili, filo )

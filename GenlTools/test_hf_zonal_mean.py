import hfdata as h

#----------------------------------
# Run from python prompt like this:
#  exec(open("./testh0.py").read())

idir='/project/amp/juliob/CAM/f09_omega/L58/2010/'
odir='/project/amp/juliob/CAM/f09_omega_zonav/L58/2010/'
xp='c6_3_59.f09_L58.ndg01'

x = h.hfdata(xp=xp,dir=idir)
y = h.hfdata(xp=xp,dir=odir)

fili=x.filename( year=2010,month=1,day=15,second=0,moniker='cam.h1')
filo=y.filename( year=2010,month=1,day=15,second=0,moniker='cam_zonav.h1')

print(fili)
print(filo)


y.zonal_mean_lonfill( fili, filo )

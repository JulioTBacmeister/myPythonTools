import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.tri as tri

#  exec(open("./testplot.py").read())

input_dir =  "/glade/campaign/cgd/projects/NCGD0051/ENSO_2010/L32/f.e22r.SAMwrf01.ne30x16.L32.NODEEP_2010_01/atm/hist/"
output_dir = "/glade/p/cesm/amwg_dev/juliob/SAMwrf/ne30x16/i/"
work_dir = "/glade/work/juliob/SAMwrf_grids/"


ifile_root = "f.e22r.SAMwrf01.ne30x16.L32.NODEEP_2010_01.cam.h1."
#ofile_root = "f.e22r.SAMwrf01.f09.L32.NODEEP_2010_01.cam.h1."
ofile_root = "f.e22r.SAMwrf01.ne30.L32.NODEEP_2010_01.cam.h1."

iy=2010
imm=6
idd=1
iss=3600

yy=str(iy).zfill(4)
mm=str(imm).zfill(2)
dd=str(idd).zfill(2)
ss=str(iss).zfill(5)
tstamp=yy+'-'+mm+'-'+dd+'-'+ss+'.nc'

ifile=input_dir+ifile_root+tstamp

#plt.ion()

a=xr.open_dataset( ifile )
h=a['PHIS']
lat=a['lat'].values
lon=a['lon'].values
hh=h[0,:].values
plt.xlim(250.,350)
plt.ylim(-60,20)
plt.tricontourf(lon, lat, hh , levels=30, cmap="RdBu_r" )

plt.show()

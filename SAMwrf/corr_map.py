import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.tri as tri


tag='ndg04'
fil1='/glade/scratch/juliob/SAMwrf_'+tag+'_ML_super_v3.nc'
tag='ndg05'
fil2='/glade/scratch/juliob/SAMwrf_'+tag+'_ML_super_v3.nc'


d1=xr.open_dataset( fil1 )
d2=xr.open_dataset( fil2 )

utn1=d1['UTEND_NDG'].values
utgw2=d2['UTEND_GWDTOT'].values
lon=d1['lon'].values
lat=d1['lat'].values

mos= np.array( [  [2010,6,1,0]
             ] )

"""
tag='ndg05'
xp='c6_3_59.ne30pg3_L32_SAMwrf.'+tag
dir='/glade/scratch/juliob/archive/'+xp+'/atm/hist/'
i=0
fili=dir+xp+'.cam.h1.'+str(mos[i,0]).zfill(4) + '-' + str(mos[i,1]).zfill(2) + '-' + str(mos[i,2]).zfill(2) + '-' +str(mos[i,3]).zfill(5) + '.nc'
a=xr.open_mfdataset( fili )
topo=xr.open_dataset( a.attrs['topography_file'] )
"""



sh=np.asarray(np.shape(utn1))
ntime=sh[0]
nlev=sh[1]
ncol=sh[2]


clevs=np.linspace(-1.,1.,num=21)
vlevs=np.arange(nlev)



corro=np.zeros( (sh[1],sh[2]) )



for L in np.arange( sh[1]-1,0,-1 ):
    #for L in np.arange( 31,30,-1 ):
    print("Level= ",L)
    for i in np.arange( sh[2] ):
        poo = np.corrcoef( x=utn1[:,L,i], y=utgw2[:,L,i] )
        corro[L,i]=poo[0,1]
        
corro=np.nan_to_num(corro,nan=0.0)

#plt.ion()

#plt.figure(0)

fig,axs=plt.subplots(2,2)
plt.xlim(270.,340)
plt.ylim(-60,20)

Lplots = [31,29,25,20]

iplot=0
for col in np.arange(2):
    for row in np.arange(2):
        ax=axs[row,col]
        Lplot=Lplots[iplot]
        ax.set_title('Level ='+str(Lplot) )
        p1=ax.tricontourf(lon, lat, corro[Lplot,:] , levels=clevs , cmap="RdBu_r" )
        pbar=plt.colorbar( p1 ,ax=ax)
        iplot=iplot+1

plt.figure(1)


fig,axs=plt.subplots(2,2)
plt.xlim(270.,340)
plt.ylim(-60,20)

Lplots = [15,10,5,2]

iplot=0
for col in np.arange(2):
    for row in np.arange(2):
        ax=axs[row,col]
        Lplot=Lplots[iplot]
        ax.set_title('Level ='+str(Lplot) )
        p1=ax.tricontourf(lon, lat, corro[Lplot,:] , levels=clevs , cmap="RdBu_r" )
        pbar=plt.colorbar( p1 ,ax=ax)
        iplot=iplot+1

plt.show()


######################
# Pretty latlon plot
######################

#MapProj = ccrs.Robinson(central_longitude=180.)
#MapProj = ccrs.PlateCarree(central_longitude=180.)
MapProj = ccrs.Orthographic(central_longitude=295.,central_latitude=-60.)
DataProj = ccrs.PlateCarree()

"""
Add axes method 
ax = fig.add_axes([xmin , ymin , X , dy , ])
"""

fig = plt.figure(figsize=(20, 12))
levs=[29,26,20,10]
   
monAs=[ 'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec' ]    
mon=10-1
cmap='terrain' #gist_ncar'
npo=0

season=monAs[mon%12]

kl0,kl1=78,82
kl0,kl1=88,92

zA = f" ( )"
clevs=[1,100,200,500,600,700,800,900,1000,1100,1200]  #2000,3000,5000,6000,7000]
scale=(1./9.81)

npo=npo+1
X=Dx
Axes = Pu.axes_def(n=npo,nxplo=2,nyplo=2) 

ax1 = fig.add_axes( Axes , projection=MapProj)
ax1.set_global()
ax1.coastlines(resolution='110m',color='black',linewidth=2)

AAxy = X['PHIS'][0,:,:]

co1=ax1.contourf(X.lon,X.lat, scale*AAxy ,transform=DataProj,levels=clevs,cmap=cmap)
ax1.set_extent([280., 310, -75., -40.], crs= ccrs.PlateCarree() )

cbar = plt.colorbar(co1, shrink=.6)
#ax1.set_title( CLUBBparm + zA, fontsize=16)
ax1.set_title( f" Rougher " , fontsize=16)

npo=npo+1
X=Dc
Axes = Pu.axes_def(n=npo,nxplo=2,nyplo=2) 

ax1 = fig.add_axes( Axes , projection=MapProj)
ax1.set_global()
ax1.coastlines(resolution='110m',color='black',linewidth=2)

AAxy = X['PHIS'][0,:,:]

co1=ax1.contourf(X.lon,X.lat, scale*AAxy ,transform=DataProj,levels=clevs,cmap=cmap)
#ax1.set_xlim((200,330))
ax1.set_extent([280., 310, -75., -40.], crs= ccrs.PlateCarree() )

cbar = plt.colorbar(co1, shrink=.6)
#ax1.set_title( CLUBBparm + zA, fontsize=16)
ax1.set_title( f" Smoother " , fontsize=16)

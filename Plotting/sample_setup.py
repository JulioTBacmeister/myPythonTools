fig = plt.figure(figsize=(20, 5))
clev=np.linspace( -60,140,num=21) 
dlev=np.linspace( -20,20,num=21) 
cmap='gist_ncar'

Axes1 = Pu.axes_def(n=1,nxplo=3,nyplo=1 ) 
ax1 = fig.add_axes( Axes1 )
co1=ax1.contourf(lat,zlev,UUc_sz ,levels=clev, cmap=cmap )
co2=ax1.contour(lat,zlev,UUc_sz ,levels=clev, colors='black')
ax1.set_title(f"Control <{exp_C}> {season.upper()}" )
ax1.set_ylim(0,82)
cb=plt.colorbar(co1)

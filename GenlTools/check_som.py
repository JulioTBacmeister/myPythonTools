>>> d3=xr.open_dataset(f3)

>>> mask=d3['mask']


>>> zz=np.where( abs(mask.values.flatten())>0 )
>>> len(zz[0])
86212
>>> len(mask)
384
>>> mask.shape
(384, 320)
>>> 384*320
122880
>>> hh=h1[0,:,:]
>>> hh=h1[0,:,:].values.flatten()
>>> oo=zz[0]
>>> fini=np.isfinite( h1[0,:,:])
>>> fini.shape
(384, 320)

>>> h3=d3['hblt']
>>> hh3=h3[0,:,:].values.flatten()
>>> plt.plot(hh3[oo])

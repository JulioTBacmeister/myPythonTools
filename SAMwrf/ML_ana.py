#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Import packages 
import sys
sys.path.append('../Plotting/')
""" Now you can imprt modules in ../Plotting"""
import xyp_plot as xyp

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.tri as tri
import matplotlib.colors as colors
from scipy import interpolate as intr

#Models
from sklearn.ensemble import RandomForestRegressor
from sklearn.tree import DecisionTreeRegressor


#Evaluation
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import cross_val_score

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import ana as a

# Some useful packages 
import importlib
import copy
import time
import pickle

# In[2]:


# Constants
Pi=3.141592653589793
R_d = 287.0 # J K-1 kg-1
C_p = 1005.7 # J K-1 kg-1
grav = 9.8 # m s-2
kappa = R_d/C_p
print(kappa)


# In[3]:


# Specify directory for data
drct='/glade/p/cesm/amwg_dev/juliob/SAMwrf/Curtains/'
# Get data (v7 has corrected ANGLL)
tag='ndg04'
fi=drct+'SAMwrf_'+tag+'_ML_super_v7.nc'
ds=xr.open_dataset( fi )


# In[4]:


##plt.ion()

# Get some ridge parameters and other unresolved topo params

MXD=ds['MXDIS']
CLN=ds['CLNGT']
ANG=ds['ANGLL']
SGH=ds['SGH']
ANX=ds['ANGLX']
print(np.shape(ANG.values.flatten()))
print(np.shape(ANX.values.flatten()))
plt.xlim(-360,360)
plt.ylim(-360,360)
#plt.plot(ANX.values.flatten(),'.')
#plt.plot(ANG.values.flatten(),'.')
plt.scatter(ANG.values.flatten(), ANX.values.flatten() )
plt.show()
# Better look at ANGLL in SE


# In[5]:


mxd=MXD.values
cln=CLN.values
ang=ANG.values
anx=ANX.values
sgh=SGH.values
print(mxd.shape)


# In[6]:


# Get nugding tendencies, horz winds and GW tendencies
# Note for intial training, results from run w/out GWD param is used - 'ndg04'

print(list(ds.variables))

Utn=ds['UTEND_NDG']*86400.
Vtn=ds['VTEND_NDG']*86400.
Utc=ds['UTEND_CORE']*86400.
Vtc=ds['VTEND_CORE']*86400.
Utgw=ds['UTEND_GWDTOT']*86400.
Vtgw=ds['VTEND_GWDTOT']*86400.
U=ds['U']
V=ds['V']
T=ds['T']

# Save-off dimensions of met-data
nT,nL,nS_0=np.shape(U)
print( nL,nT,nS_0 )


# In[7]:


#plt.plot(Utn[31,0,:])


# In[8]:


# Calculate 3D pressure field for full [nT,nL,nS_0] range
ps=ds['PS']
# These have dims=[nT,nL] ...!!!
hyam=ds['hyam']
hybm=ds['hybm']

#plev=d4['lev']
#gps=np.average( ps, axis=0 )
#ghya=np.average( hyam , axis=0 )
#ghyb=np.average( hybm , axis=0 )

p3=a.press(PS=ps,hybm=hybm,hyam=hyam )
print(np.shape(p3),np.shape(T))


# In[9]:


hybi=ds['hybi']
hyai=ds['hyai']

p3e=a.press(PS=ps,hybm=hybi,hyam=hyai )


# In[10]:


# Testing random number stuff (local cell w/ no side effects)
rng=np.random.default_rng(seed=42)
uuu=rng.uniform(0,2,1000)

huuu,bins=np.histogram(uuu)
#plt.stairs( huuu, bins  )


# # WARNING !!!!!

# In[11]:


# DO NOT DO THIS !!!!! :
#   > theta=T.values
#   > te = T.values
#   > theta = calculations .... 
#
# '=' in Python is not really a copy. Acts like a pointer.
# So, in above 'te' is modified by calculations 
# even though it is never on the LHS
#
# Need to use, e.g.,
#   > theta =copy.deepcopy(T.values)


# In[12]:


# Calculate dry Theta at mid-levels
# 

theta=np.zeros( (nT, nL, nS_0 )    ) # T.values
te=T.values
#plt.plot(te[1,:,2000])

print( id(te) )
print( id(theta) )

print(np.shape(hybm))
for iT in np.arange(nT):
    theta[iT,:,:] = (( 100_000. / p3[iT,:,:]  )**kappa ) * te[iT,:,:]
    
#plt.plot(te[1,:,2000])


# In[ ]:


# Calculate dry air density at mid-levels (kg/m+3)

rho_d = p3 / (R_d * te )


# In[ ]:


print( np.arange(nL) )
print(' .. ')
print( np.arange( nL-1,-1,-1))


# In[ ]:


# Integrate hydrostatic relation to get heights (i.e. geoph/grav) at edges

zhgte=np.zeros( (nT,nL+1,nS_0) )
for iL in np.arange(nL-1,-1,-1):
    zhgte[:,iL,:] = zhgte[:,iL+1,:] + \
                    ( p3e[:,iL+1,:]-p3e[:,iL,:]  )  \
                    /( grav * rho_d[:,iL,:] )
    print(iL,end=',')


# In[ ]:


# Simple averaging to get mid-level heights from edge heights

zhgtm = 0.5 * (zhgte[:,0:nL,:]+zhgte[:,1:nL+1,:])
print(np.shape(zhgtm))

#plt.plot( theta[1,:,2000] , zhgtm[1,:,2000] )


# In[17]:


# Calculate Theta values at edges by using 1D linear interpolation in height (w/
# extrapolation) from mid-level heights and Theat 

thetae=np.zeros( (nT,nL+1,nS_0) )
for iT in np.arange(nT):
    print(iT,end=',')
    for iS in np.arange(nS_0):
        fint=intr.interp1d( x=zhgtm[iT, :, iS],y=theta[iT,:,iS] , fill_value='extrapolate'  )
        thetae[iT, :, iS] = fint( zhgte[iT, :, iS ] )


# In[18]:

"""
plt.plot( thetae[1,:,2000] , zhgte[1,:,2000] )

plt.plot( theta[1,:,2000] , zhgtm[1,:,2000], 'x' )

plt.ylim(0,1000)
plt.xlim(280,320)
#plt.plot(rho_d[1,:,2000] )
"""

# In[19]:


# Calculate "N-squared" from Theta on edges
# Simple centered differences
nsq = np.zeros( (nT, nL, nS_0 ) )
for iL in np.arange( start=1,stop=nL-1):
    nsq[:,iL,:] = grav * ( ( thetae[:,iL-1,:]-thetae[:,iL+1,:] ) / \
                           ( zhgte[:,iL-1,:] -zhgte[:,iL+1,:]  ) ) / \
                                theta[:,iL,:] 


# Take boundary values from adjacent levels
nsq[:,nL-1,:]=nsq[:,nL-2,:]
nsq[:,0,:]=nsq[:,1,:]


# In[20]:


# Calculate "N" (Brunt-Vaisalla freq) from N**2.
# Here account for negative (unstable) stratification using np.where
#----------------------------------------
# Once again, np.where works like this:
#  B = np.where( {condition on A} , {value where condition=True}, {value where condition=False} )
#  B will have the same shape as A.
#  2nd and 3rd arguments can be scalars or shaped like A.

n_bv=np.where( nsq>=0., nsq , -nsq )
stab=np.where( nsq>=0., 1.0 , -1.0 )
n_bv = stab * np.sqrt( n_bv )

"""
plt.ylim(0,30000.)
#plt.xlim(-0.0001,0.0001)
plt.xlim(-0.01,0.04)
#plt.plot(  nsq[1,:,2000] , zhgtm[1,:,2000] )
plt.plot(  n_bv[1,:,2000] , zhgtm[1,:,2000], '+' )
"""

# ## Begin setting up data for ML
# 

# ## Here we isolate attention to regions with topography

# In[21]:


# Pick out grid cells with MXD[0,:] bigger than 50m.
#
# Note here np.where is working more like IDL where.
#
# With no 2nd and 3rd args np.where documentaion says:
#   When only condition is provided, this function is a shorthand for np.asarray(condition).nonzero(). 
#   Using nonzero directly should be preferred, as it behaves correctly for subclasses. The rest of 
#   this documentation covers only the case where all three arguments are provided.


oo=np.where(mxd[0,:]>50.)
print(np.shape(oo))



mxd=mxd[0,oo[0][:]]
cln=cln[0,oo[0][:]]
ang=ang[0,oo[0][:]]
angrad=ang*Pi/180.

cosrdg=np.cos( angrad )
sinrdg=np.sin( angrad )

"""
fig=plt.figure( figsize =(20,9))
ax=fig.add_subplot(2,2,1)
po=ax.plot( mxd, '.')
po=ax.set_title( 'MXDIS' ,fontsize=12, loc='center')
ax=fig.add_subplot(2,2,2)
po=ax.plot( cln, '.')
po=ax.set_title( 'CLNGT' ,fontsize=12, loc='center')
ax=fig.add_subplot(2,2,3)
po=ax.plot( ang,cosrdg, '.')
po=ax.set_title( 'Cos(ANGLL)' ,fontsize=12, loc='center')
ax=fig.add_subplot(2,2,4)
po=ax.plot( ang,sinrdg, '.')
po=ax.set_title( 'Sin (ANGLL)' ,fontsize=12, loc='center')
"""

# In[22]:


# Numpy arrays are [nT,nL,ns]
# Select columns with MXD>h_crit nS_0=>nS
# and restrict to bottom 10-layers
# Transpose to [nL,nt,ns_crit]
# Then reshape to [nL,nT*ns_crit]
u=U.values
u=u[:,21:,oo[0][:]]
u=np.transpose(u, (1,0,2) )
nL,nT,nS=np.shape(u)
u=np.reshape( u, (nL,nT*nS ) )/10.

v=V.values
v=v[:,21:,oo[0][:]]
v=np.transpose(v, (1,0,2) )
v=np.reshape( v, (nL,nT*nS ) )/10.

te=T.values
te=te[:,21:,oo[0][:]]
te=np.transpose(te, (1,0,2) )
te=np.reshape( te, (nL,nT*nS ) )/10.


# In[23]:


utn=Utn.values
print(np.shape(utn))



#plt.plot(utn[100:500,30,:].flatten())


# In[24]:


# Numpy arrays are [nT,nL,ns]
# Select columns with MXD>h_crit ns=>ns_crit
# and restrict to bottom 10-layers
# Transpose to [nL,nt,ns_crit]
# Then reshape to [nL,nT*ns_crit]
utn=Utn.values
utn=utn[:,21:,oo[0][:]]
utn=np.transpose(utn, (1,0,2) )
utn=np.reshape( utn, (nL,nT*nS ) )/10.

vtn=Vtn.values
vtn=vtn[:,21:,oo[0][:]]
vtn=np.transpose(vtn, (1,0,2) )
vtn=np.reshape( vtn, (nL,nT*nS ) )/10.




# In[25]:


mxdx=np.zeros([1,nT,nS])
clnx=np.zeros([1,nT,nS])
cosrx=np.zeros([1,nT,nS])
sinrx=np.zeros([1,nT,nS])
for iT in np.arange(nT):
    mxdx[0,iT,:]=mxd
    clnx[0,iT,:]=cln
    cosrx[0,iT,:]=cosrdg
    sinrx[0,iT,:]=sinrdg
mxdx = np.reshape( mxdx, (1,nT*nS ) )/1000.
clnx = np.reshape( clnx, (1,nT*nS ) )/100.


# In[26]:


print(np.shape(mxdx))
#plt.plot(cln.view())


# ## Creating training data (A) and target ('label') data (B)

# In[27]:


A=np.r_[u,v,mxdx,clnx]
B=np.r_[utn,vtn]


# In[28]:


print(np.shape(A))


# In[29]:


#A=A[:,0:20000]
#B=B[:,0:20000]
A=np.transpose(A)
B=np.transpose(B)

Ashp = np.shape(A)


# In[30]:


idxs=np.arange( Ashp[0] )
print( np.shape(idxs) )
Ridxs = copy.deepcopy(idxs) 


tic = time.perf_counter()
np.random.shuffle(Ridxs)
toc = time.perf_counter()
ShuffleTime = f"Shuffled indices in {toc - tic:0.4f} seconds"

print(ShuffleTime)
print( np.shape(Ridxs) )

#plt.scatter( idxs, Ridxs)


# In[31]:


A_r = A[ Ridxs, :]
B_r = B[ Ridxs, :]
#print( Ridxs[10] )
#print( A_r[10,:])
#print( A[Ridxs[10],:])


# In[32]:


A_train = A_r[0:788200,:]
B_train = B_r[0:788200,:]
A_test  = A_r[788201:,:]
B_test  = B_r[788201:,:]

A_mini = A_r[0:155_000,:]
B_mini = B_r[0:155_000,:]


# In[33]:


forest_reg = RandomForestRegressor(n_estimators=100, random_state=42)


# In[34]:


tic = time.perf_counter()
#forest_reg.fit(A_mini, B_mini)
forest_reg.fit(A_train, B_train)
toc = time.perf_counter()
TrainingTime = f"Trained model in {toc - tic:0.4f} seconds"
print(TrainingTime)


# In[35]:


B_pred=forest_reg.predict(A_test)


# In[36]:


plt.scatter(B_pred.flatten(),B_test.flatten())
plt.show()

# In[37]:


poo = np.corrcoef( x=  B_pred.flatten(), y = B_test.flatten() )
print(poo)


# In[38]:


# 155_000 rows : Trained model in 562.1953 seconds
# 55_000  rows: r=0.59
# 155_000 rows: r=0.636



#How to save a RandomForest model


filename = "random_forest_full.pkl"

# save model
pickle.dump(forest_reg , open(filename, "wb"))

#
# load model
#loaded_model = pickle.load(open(filename, "rb"))

# you can use loaded model to compute predictions
#y_predicted = loaded_model.predict(X)



# In[ ]:





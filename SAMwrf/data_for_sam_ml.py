import numpy as np
import xarray as xr



#/glade/scratch/juliob/archive/c6_3_59.ne30pg3_L32_SAMwrf.ndg04/atm/hist/c6_3_59.ne30pg3_L32_SAMwrf.ndg04.cam.h1.*.nc

tag='ndg04'
xp='c6_3_59.ne30pg3_L32_SAMwrf.'+tag
dir='/glade/scratch/juliob/archive/'+xp+'/atm/hist/'
odir='/glade/p/cesm/amwg_dev/juliob/SAMwrf/Curtains/'

topofile='/glade/p/cgd/amp/juliob/bndtopo/latest/ne30pg3_gmted2010_modis_bedmachine_nc3000_Laplace0100_20230105.nc'

mos= np.array( [  [2010,6]
                 ,[2010,7]
                 ,[2010,8] 
                 ,[2010,9] 
                 ,[2010,10] 
                 ,[2010,11] 
             ] )

print(mos.shape)

smos=mos.shape

lfilo=[]

for i in np.arange( smos[0] ):
    fili=dir+xp+'.cam.h1.'+str(mos[i,0]).zfill(4) + '-' + str(mos[i,1]).zfill(2) +'*.nc'
    a=xr.open_mfdataset( fili )
    lons=a['lon'].values[0,:]
    lats=a['lat'].values[0,:]
    oo=np.where( (lats>-60)&(lats<20.)&(lons>270.)&(lons<340.0) )
    ooo=np.asarray(oo)[0,:]
    print(" Opened .. ",fili)

    #topo=xr.open_dataset( a.attrs['topography_file'] ) #Overridden above by corrected topo file
    topo=xr.open_dataset( topofile ) #Overridden above by corrected topo file

    """
    Set up 'dict' for xarray dataset.  This can then
    turned into basis for dataset, which can then be 
    added to.  'dict' sets up dimensions and puts 
    lat and lon vectors
    """
    d = { 
        'lon':{'dims':('ncol'), 'data':lons[ooo] },
        'lat':{'dims':('ncol'), 'data':lats[ooo] },
        'time':{'dims':('time'), 'data':a['time'].values },
        'lev':{'dims':('lev'), 'data':a['lev'].values },
        'ilev':{'dims':('ilev'), 'data':a['ilev'].values },
        'nbnd':{'dims':('nbnd'), 'data':a['nbnd'].values },
        'irdg':{'dims':('nrdg') , 'data':np.arange(16) },
        }
    b = xr.Dataset.from_dict(d)
    """
    Add some attributes
    """
    b.attrs={'topography_file':a.attrs['topography_file']}
    """
    NOw add variables to base dataset one at a time
    """
    b['hyam']=a['hyam']
    b['hybm']=a['hybm']
    b['hyai']=a['hyai']
    b['hybi']=a['hybi']


    b['PHIS']=topo['PHIS'][ooo]
    b['SGH']=topo['SGH'][ooo]
    b['MXDIS']=topo['MXDIS'][:,ooo]
    b['ANGLL']=topo['ANGLL'][:,ooo]
    b['ANGLX']=topo['ANGLX'][:,ooo]
    b['CLNGT']=topo['CLNGT'][:,ooo]
    b['ANISO']=topo['ANISO'][:,ooo]
    b['ISOVAR']=topo['ISOVAR'][ooo]

    b['PS']=a['PS'][:,ooo]

    b['Target_U']=a['Target_U'][:,:,ooo]
    b['U']=a['U'][:,:,ooo]
    b['UTEND_NDG']=a['UTEND_NDG'][:,:,ooo]
    b['UTEND_GWDTOT']=a['UTEND_GWDTOT'][:,:,ooo]
    b['UTEND_CORE']=a['UTEND_CORE'][:,:,ooo]
    b['Target_V']=a['Target_V'][:,:,ooo]
    b['V']=a['V'][:,:,ooo]
    b['VTEND_NDG']=a['VTEND_NDG'][:,:,ooo]
    b['VTEND_GWDTOT']=a['VTEND_GWDTOT'][:,:,ooo]
    b['VTEND_CORE']=a['VTEND_CORE'][:,:,ooo]

    b['T']=a['T'][:,:,ooo]

    
    filo=odir+'SAMwrf_'+tag+'_ML_'+str(mos[i,0]).zfill(4) + '-' + str(mos[i,1]).zfill(2) +'.nc'
    print(" will try to write curtain to ",filo)
    b.to_netcdf( filo )
    print(" wrote curtain to ",filo)
    lfilo.append( filo ) 

print(lfilo)
"""
 'different' keyword prevents variables that don't change 
 from one file to the next from being replicated
"""
ds=xr.open_mfdataset( lfilo, data_vars='different' )

filo=odir+'SAMwrf_'+tag+'_ML_super_v7.nc'
ds.to_netcdf( filo )


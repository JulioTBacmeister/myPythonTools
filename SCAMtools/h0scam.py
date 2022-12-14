class h0scam:

    def __init__(self, xp, dir, base='NA', machine='NA'):
        self.case    = xp
        self.archdir = dir
        self.basecase = base
        self.machine = machine

    def curtain(self,fld):
        import numpy as np
        import xarray as xr
        import glob
        import txtutil as tx
        import os
        import subprocess as sp

        host=os.environ['HOST']
        user=os.environ['USER']
 
        if ('izumi' in host):
            self.machine="izumi"
        elif ('cheyenne' in host):
            self.machine="cheyenne"        
        elif ('thorodin' in host):
            self.machine="thorodin"        

        xp=self.case
        dir=self.archdir
        basedir=self.basecase
        

        fili=dir+'/atm_in'
        nhtfrq = tx.nmlread( fili, 'nhtfrq' )
        nhtfrq = nhtfrq.split(',')
        h0frq  = int( nhtfrq[0] )
        
        fili=dir+'/nuopc.runconfig'
        freqc = int( tx.nmlread( fili, 'atm_cpl_dt' ) )

        if (h0frq > 0):
            freqw=freqc*h0frq

        freqs = str( freqw )
        freqs = freqs.strip()+'S'
        print('write interval='+freqs)
 
        fl = sorted( glob.glob( dir +'/*cam.h0*') )
        nf = len( fl )

        #fl=fl[0:4]
        ird=0
        for f in fl:
            print(f)
            try:
                a=xr.open_dataset( f )
                print("Successfully opened"+f)
                nx=a.lon.size
                ny=a.lat.size
                nl=a.lev.size
                nli=a.ilev.size
                nt=a.time.size
  
                # Start looping over multiple fields here .... 
                aa=a[fld] #.isel(time=0)
                print(aa.dims)

                if ('lon' in aa.dims) and ('lat' in aa.dims) and ('lev' not in aa.dims) and ('ilev' not in aa.dims):
                    print(fld+' is a surface var ' )
                    varType = 'surface'
                if ('lon' in aa.dims) and ('lat' in aa.dims) and ('lev' in aa.dims):
                    print(fld+' is a profile var ' )
                    varType = 'profile'
                if ('lon' in aa.dims) and ('lat' in aa.dims) and ('ilev' in aa.dims):
                    print(fld+' is an Iprofile var ' )
                    varType = 'iprofile'

                if (ird == 0):
                    timeData=a.time.data
                    bigTime = xr.cftime_range( timeData[0]  , periods=nt*nf, freq=freqs , calendar="noleap" )
                    if varType == 'surface':
                        dummy=np.zeros( [nt*nf] )
                        #dummy=dummy.reshape( nt*nf ,1,1)
                        #cu=xr.DataArray( dummy , coords=[bigTime,a.lat,a.lon], dims=['time','lat','lon'] , name=fld )                    
                    if varType == 'profile':
                        dummy=np.zeros( [nt*nf,nl] )
                        #dummy=dummy.reshape( nt*nf,nl ,1,1)
                        #cu=xr.DataArray( dummy , coords=[bigTime,a.lev,a.lat,a.lon], dims=['time','lev','lat','lon'], name=fld  )                    
                    if varType == 'iprofile':
                        dummy=np.zeros( [nt*nf,nli] )
                        #dummy=dummy.reshape( nt*nf,nli ,1,1)
                        #cu=xr.DataArray( dummy , coords=[bigTime,a.ilev,a.lat,a.lon], dims=['time','ilev','lat','lon'], name=fld  )                    
                        
                    print('Prepped numpy ARRAY for data with time')

                if varType == 'surface':
                    #cu.values[  ird*nt :(ird+1)*nt , 0, 0] = aa[:,0,0]  
                    dummy[  ird*nt :(ird+1)*nt ] = aa[:,0,0]  
 
                if varType == 'profile' or varType == 'iprofile':
                    #cu.values[  ird*nt :(ird+1)*nt , : , 0, 0] = aa[:,:,0,0]  
                    dummy[  ird*nt :(ird+1)*nt , : ] = aa[:,:,0,0]  

                ird=ird+1

            except ValueError:
                print('******** VALUE ERROR *************')
                print('File \n'+f+'\n probably no good')
                if (ird == 0):
                    exit('First dataset not valid')
                else:
                    ird=ird+1


        if varType == 'surface':                
            #dummy=dummy.reshape( nt*nf ,1,1)
            #cu=xr.DataArray( dummy , coords=[bigTime,a.lat,a.lon], dims=['time','lat','lon'] , name=fld )                    
            cu=xr.DataArray( dummy , coords=[bigTime], dims=['time'] , name=fld )                    
        if varType == 'profile':                
            #dummy=dummy.reshape( nt*nf,nl ,1,1)
            #cu=xr.DataArray( dummy , coords=[bigTime,a.lev,a.lat,a.lon], dims=['time','lev','lat','lon'] , name=fld )                    
            cu=xr.DataArray( dummy , coords=[bigTime,a.lev], dims=['time','lev'] , name=fld )                    
        if varType == 'iprofile':                
            #dummy=dummy.reshape( nt*nf,nli ,1,1)
            #cu=xr.DataArray( dummy , coords=[bigTime,a.ilev,a.lat,a.lon], dims=['time','ilev','lat','lon'] , name=fld )                    
            cu=xr.DataArray( dummy , coords=[bigTime,a.ilev], dims=['time','ilev'] , name=fld )                    

        if (self.machine=='cheyenne'):
            wrtdir='/glade/scratch/'+user+'/'
        elif (self.machine=='thorodin') or (self.machine=='izumi'):
            if (self.basecase != 'NA'):
                wrtdir='/project/amp/juliob/scam/archive/pyProc/'  #'+self.basecase+'/'+xp+'/atm/pyProc/'
            else:
                wrtdir='/project/amp/juliob/scam/archive/pyProc/'
        else:
            wrtdir='./'
  
        cmd1="mkdir -p "+wrtdir
        print('File write: '+wrtdir+xp+'_'+fld+'.nc')
        sp.run(cmd1,shell=True)
        cu.to_netcdf(wrtdir+xp+'_'+fld+'.nc')

        return cu


    def ncmerge(self):
        import numpy as np
        import xarray as xr
        import glob
        import txtutil as tx

        

        xp=self.case
        dir=self.archdir
             
        fl = sorted( glob.glob( dir +'/*cam.h0*') )
        nf = len( fl )

        flt = fl[0:10]
 
        #ds=xr.open_mfdataset( flt , concat_dim='time')
        #print("merged nc files")
        #ds.to_netcdf(path = dir +'/merged_h0.nc')
        #exit()

        ird=0
        for f in flt:
            print(f)
            try:
                a=xr.open_dataset( f )
                print("Successfully opened"+f)

                if (ird == 0):
                    b=a
                elif (ird>0):
                    b.merge(a,join='inner' )
                    b.to_netcdf( 'testmerge.nc' )

                ird=ird+1

            except ValueError:
                print('******** VALUE ERROR *************')
                print('File \n'+f+'\n probably no good')
                if (ird == 0):
                    exit('First dataset not valid')
                else:
                    ird=ird+1
 
        return

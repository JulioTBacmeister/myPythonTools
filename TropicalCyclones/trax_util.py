import sys
# import modules in other directories
sys.path.append('../Utils/')

import numpy as np
import xarray as xr

#import MyConstants as Con
import constants as Con



#===========================================================
# Class to allow things to be accessed via dict.thing syntax
# as well as dict['thing'] syntax. They are equivalent.
#===========================================================
class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


#===========================================================
#               functions
#===========================================================

# ---------------------------
# --- Basic reader
# ---------------------------
def readtrx( fileN ,power_wind=1.0):
    
    dS = xr.open_dataset( fileN )
    #print( list(dS.variables) )
    
    date=dS.date.values
    datesec=dS.datesec.values
    year = (date/10000).astype(int)
    month = ((date-10000*year)/100).astype(int)
    day = (date - 10000*year-100*month).astype(int)
    hour = (datesec/3600).astype(int)
    
    wind = dS.wind.values * power_wind

    speed   = translation_speed( dS.lon.values , dS.lat.values ) 
    intfW3,intfP3 = intensification( wind, dS.pres.values , dTime=3. ) 
    intfW24,intfP24 = intensification( wind, dS.pres.values , dTime=24. ) 
    
    struc = {'date':date,'year':year,'month':month,'day':day,'hour':hour,
            'basin':dS.basin.values, 'lon':dS.lon.values , 'lat':dS.lat.values ,
            'wind':wind , 'speed':speed , 
            'intfW3':intfW3 , 'intfP3':intfP3 ,
            'intfW24':intfW24 , 'intfP24':intfP24 ,
            'pres':dS.pres.values, 'sub_basin': 0*dS.basin.values - 99 }
    
    Atruc = AttrDict( struc )
    
    return Atruc


# ---------------------------
# --- filename constructors
# ---------------------------
def pdfname(ens=1,sst='',beginyear=1979,endyear=2012,special='',justBaseName=False):
    dirA='/glade/campaign/cgd/amp/juliob/TC-cesm1/10mWinds_tracfiles/pd/'
    bname='f.e13.FAMIPC5.ne120_ne120.1979_2012'
    yr0A = str(beginyear).zfill(4)
    yr1A = str(endyear).zfill(4)
    yrsA = yr0A+'-'+yr1A
    if (ens>0):
        ensA=str(ens).zfill(3)
    else:
        ensA=''
    
    if (special != ''):
        spclA = special+'.'
    else:
        spclA = ''
        
       

    
    if (justBaseName==True):
        fname=bname+'.'+spclA+ensA+'.'
    else:
        fname = bname+'.'+spclA+ensA+'.'+yrsA+'_tracfile.nc'
        fname=dirA+fname
    
    return fname
    
def rcp85fname(ens=1,sst='',beginyear=2070,endyear=2099,special='',justBaseName=False):
    dirA='/glade/campaign/cgd/amp/juliob/TC-cesm1/10mWinds_tracfiles/rcp85/'
    bname='f.e13.FAMIPC5.ne120_ne120.RCP85_2070_2099'
    yr0A = str(beginyear).zfill(4)
    yr1A = str(endyear).zfill(4)
    yrsA = yr0A+'-'+yr1A
    if (ens>0):
        ensA=str(ens).zfill(3)
    else:
        ensA=''
        
    if ((sst != '') and (sst != 'sst1')):
        sstA = '_'+sst
    else:
        sstA = ''
        
    if (special != ''):
        spclA = '.'+special
    else:
        spclA = ''
        
       
    if (justBaseName==True):
        fname = bname+sstA+'.'+ensA+spclA+'.'
    else:
        fname = bname+sstA+'.'+ensA+spclA+'.'+yrsA+'_tracfile.nc'
        fname=dirA+fname
    
    return fname

#----------------------------
# --- Utilities
# ---------------------------
def basinACE(wind,basin,power_wind=1.0):
    
    nstorm,nt=np.shape( wind )
    Nb = ( np.nanmax(basin)-np.nanmin(basin) + 1 ).astype(int)
    bace = np.zeros( (Nb , nt) )
    awind = np.nan_to_num( wind )
    
    w10m_kts = power_wind*1.944*awind
    aace = 1e-4 * w10m_kts**2 
    
    bmin=np.nanmin(basin).astype(int)
    bmax=np.nanmax(basin).astype(int)
    
        
    # Yes, this is an awkward way to do this 
    for b in np.arange( start=bmin , stop=bmax+1 ):
        ib = b - bmin
        for s in np.arange(nstorm):
            oo=np.where( basin[s,:]==b )
            if (len(oo[0])>0):
                #print(oo[0])
                indx=oo[0][:]
                #print(ib,np.max(aace[s,indx]))
                bace[ib,indx] = bace[ib,indx]+aace[s,indx]

    Bace = np.zeros(Nb)
    # This is a clearer way to do it
    for b in np.arange( start=bmin , stop=bmax+1 ):
        ib = b - bmin
        oo=np.where( basin.flatten()==b )
        if (len(oo[0])>0):
                indx=oo[0][:]
                Bace[ib] = np.sum(aace.flatten()[indx])
 

    return Bace

def basinACEyear(wind,basin,year,power_wind=1.0):

    nstorm,nt=np.shape( wind )
    Nb = ( np.nanmax(basin)-np.nanmin(basin) + 1 ).astype(int)
    bmin=np.nanmin(basin).astype(int)
    bmax=np.nanmax(basin).astype(int)

    awind = np.nan_to_num( wind )
    
    w10m_kts = power_wind*1.944*awind
    aace = 1e-4 * w10m_kts**2 
    
    oo=np.where( year.flatten() > 0)
    yrmin=np.nanmin( year.flatten()[oo[0]] )
    yrmax=np.nanmax( year.flatten()[oo[0]] )

    print( yrmin,yrmax )
    Nyr = (yrmax - yrmin + 1).astype(int)
    
    BaceYr = np.zeros((Nb,Nyr))
    Years = np.zeros( Nyr )
    print(np.shape(BaceYr))
    # This is a clearer way to do it
    for y in np.arange( start=yrmin , stop=yrmax+1 ):
        iy = y - yrmin
        Years[ iy ] = y
        for b in np.arange( start=bmin , stop=bmax+1 ):
            ib = b - bmin
            oo=np.where( ( basin.flatten()==b ) & ( year.flatten()==y ) )
            if (len(oo[0])>0):
                indx=oo[0][:]
                BaceYr[ib,iy] = np.sum(aace.flatten()[indx])
 

    return BaceYr,Years


def basinACEmonth(wind,basin,year,month,power_wind=1.0):

    nstorm,nt=np.shape( wind )
    #Nb = ( np.nanmax(basin)-np.nanmin(basin) + 1 ).astype(int)
    bmin=np.nanmin(basin).astype(int)
    bmax=np.nanmax(basin).astype(int)
    # Always have -1 as minimum basin id number
    #-----------------------------------------
    if (bmin<-1):
        bmin = -1
    Nb = ( bmax - bmin + 1 ).astype(int)
    
    
    awind = np.nan_to_num( wind )
    
    w10m_kts = power_wind*1.944*awind
    aace = 1e-4 * w10m_kts**2 
    
    oo=np.where( year.flatten() > 0)
    yrmin=np.nanmin( year.flatten()[oo[0]] )
    yrmax=np.nanmax( year.flatten()[oo[0]] )

    print( yrmin,yrmax )
    Nyr = (yrmax - yrmin + 1).astype(int)
    
    BaceMon = np.zeros((Nb,Nyr,12))
    Years = np.zeros( Nyr )
    Months= np.zeros( 12 )
    print(np.shape(BaceMon))
    # This is a clearer way to do it
    for y in np.arange( start=yrmin , stop=yrmax+1 ):
        iy = y - yrmin
        Years[ iy ] = y
        for m in np.arange(start=1,stop=13):
            im = m - 1
            Months[ im ] = m
            for b in np.arange( start=bmin , stop=bmax+1 ):
                ib = b - bmin
                oo=np.where( ( basin.flatten()==b ) & ( year.flatten()==y ) & (month.flatten()==m) )
                if (len(oo[0])>0):
                    indx=oo[0][:]
                    BaceMon[ib,iy,im] = np.sum(aace.flatten()[indx])
 

    return BaceMon,Years,Months

def basinACEyearHemi(wind,basin,year,month,power_wind=1.0):

    #         -1.       0.        1.         2            3             4              5            6 
    bnames=['None','N \nAtl','S \n Atl' , 'NW \n Pac', 'NE \n Pac', 'SW \n Pac','N \n Indian','S \nIndian']

    nstorm,nt=np.shape( wind )
    Nb = ( np.nanmax(basin)-np.nanmin(basin) + 1 ).astype(int)
    bmin=np.nanmin(basin).astype(int)
    bmax=np.nanmax(basin).astype(int)

    awind = np.nan_to_num( wind )
    
    w10m_kts = power_wind*1.944*awind
    aace = 1e-4 * w10m_kts**2 
    
    oo=np.where( year.flatten() > 0)
    yrmin=np.nanmin( year.flatten()[oo[0]] )
    yrmax=np.nanmax( year.flatten()[oo[0]] )

    fyear = year + month/12.
    
    print( "Years  ",yrmin,yrmax )
    print( "Basins ",bmin,bmax )
    Nyr = (yrmax - yrmin + 1).astype(int)
    
    BaceYr = np.zeros((Nb,Nyr))
    Years = np.zeros( Nyr )
    print(np.shape(BaceYr))
    # This is a clearer way to do it
    for y in np.arange( start=yrmin , stop=yrmax+1 ):
        iy = y - yrmin
        Years[ iy ] = y
        for b in np.arange( start=bmin , stop=bmax+1 ):
            indx=None
            ib = b - bmin
            if ((b<=0) or (b==2) or (b==3) or (b==5) ):
                oo=np.where( ( basin.flatten()==b ) & ( year.flatten()==y ) )
                if (len(oo[0])>0):
                    indx=oo[0][:]
            elif ((b==1) or (b==4) or (b==6) ):  #Southern Hemisphere
                oo=np.where( ( basin.flatten()==b ) & ( fyear.flatten()>y-0.5) & ( fyear.flatten()<=y+0.5 ) )
                if (len(oo[0])>0):
                    indx=oo[0][:]
            if (indx is not None):
                BaceYr[ib,iy] = np.sum(aace.flatten()[indx])
                BaceYr[ib,iy] = np.sum(aace.flatten()[indx])
 

    return BaceYr,Years


def basinNUMyearHemi(wind,basin,year,month,power_wind=1.0,category='TS+'):


    if (category=='TD+'):
        thresh0=0.
        thresh1=10000.
    if (category=='TS+'):
        thresh0=18.
        thresh1=10000.
    if (category=='Cat1+'):
        thresh0=33.
        thresh1=10000.
    if (category=='Cat2+'):
        thresh0=43.
        thresh1=10000.
    if (category=='Cat3+'):
        thresh0=50.
        thresh1=10000.
    if (category=='Cat4+'):
        thresh0=58.
        thresh1=10000.
    if (category=='Cat5+'):
        thresh0=70.
        thresh1=10000.
    CatThreshhold  = thresh0*1.944

    #         -1.       0.        1.         2            3             4              5            6 
    bnames=['None','N \nAtl','S \n Atl' , 'NW \n Pac', 'NE \n Pac', 'SW \n Pac','N \n Indian','S \nIndian']

    nstorm,nt=np.shape( wind )
    Nb = ( np.nanmax(basin)-np.nanmin(basin) + 1 ).astype(int)
    bmin=np.nanmin(basin).astype(int)
    bmax=np.nanmax(basin).astype(int)

    awind = np.nan_to_num( wind )
    
    w10m_kts = power_wind*1.944*awind
    aace = 1e-4 * w10m_kts**2 
    
    oo=np.where( year.flatten() > 0)
    yrmin=np.nanmin( year.flatten()[oo[0]] )
    yrmax=np.nanmax( year.flatten()[oo[0]] )

    fyear = year + month/12.
    
    print( "Years  ",yrmin,yrmax )
    print( "Basins ",bmin,bmax )
    Nyr = (yrmax - yrmin + 1).astype(int)
    
    NumYr = np.zeros((Nb,Nyr))
    Years = np.linspace(yrmin,yrmax, num=Nyr )

    # This is a clearer way to do it
    for s in np.arange(nstorm):
        if max( w10m_kts[s,:] >= CatThreshhold ):
            b  = basin[s,0].astype(int)
            y  = year[s,0].astype(int)
            ib = b - bmin
            iy = y - yrmin
            NumYr[ib,iy]=NumYr[ib,iy]+1

    return NumYr,Years

def add_sub_basin(trx,X,Y,sub_basin_number):
    
    trx_o = trx
    

    nstorm,nt=np.shape( trx.wind )
    
    for s in np.arange( nstorm ):
        for t in np.arange( nt ):
            lon0=trx.lon[s,t]
            lat0=trx.lat[s,t]
        
            In_sub_basin = winding_number(lon0, lat0, X, Y)
            if (In_sub_basin==True):
                trx_o.sub_basin[s,t]=sub_basin_number
                     
    
    return trx_o

def translation_speed(lon,lat):
    
    Rer  = Con.Rearth()
    pi  = Con.pi()
    d2r = pi/180.
    
    nstorm,nt=np.shape( lon )
    speed=np.zeros( (nstorm,nt) )
    delxy=np.zeros( (nstorm,nt) )
    dTime = 3. * 3600. # SHould make this a calculation not an assumption
    print(f"Rer {Rer}")
    
    for s in np.arange( nstorm ):
        for t in np.arange( start=1,stop=nt-1 ):
            lat0=d2r*lat[s,t]
            lonI=d2r*lon[s,t+1]
            latI=d2r*lat[s,t+1]
            lonD=d2r*lon[s,t-1]
            latD=d2r*lat[s,t-1]
        
            dx = np.cos(lat0) * (lonI - lonD)
            dy = (latI - latD )
            delxy[s,t] = np.sqrt( dx**2 + dy**2 ) #np.sqrt( dx**2 + dy**2 ) / dTime
            speed[s,t] = Rer*delxy[s,t] / dTime #np.sqrt( dx**2 + dy**2 ) / dTime
        
    
    return speed

def intf_by_24hr(wind,pres,year,month,day,hour):
    #Returns rate in m s-1 h-1, hPa h-1
    
    nstorm,nt=np.shape( hour )
    intfW=np.empty( (nstorm,nt) )
    intfW.fill(np.nan)
    intfP=np.empty( (nstorm,nt) )
    intfP.fill(np.nan)

    dTime = 1. # days # 24. # Hours. SHould make this a calculation not an assumption

    for s in np.arange( nstorm ):
        oo=np.where( ( hour[s,:]==0 ) & (wind[s,:]>0.) )
        W_oo = wind[ s, oo[0] ]
        P_oo = pres[ s, oo[0] ]
        nt_oo = np.size( oo[0] )
        if (nt_oo > 2):
            for t in np.arange( start=1,stop=nt_oo-1 ):
                wind0=W_oo[t]
                pres0=P_oo[t]
                windD=W_oo[t-1]
                presD=P_oo[t-1]

                intfW[s,t-1] = (wind0 - windD)/dTime       
                intfP[s,t-1] = (pres0 - presD)/dTime       
    
    return intfW,intfP

def intensification(wind,pres,dTime=3.):
    #Returns rate in m s-1 h-1, hPa h-1
    
    nstorm,nt=np.shape( wind )
    intfW=np.zeros( (nstorm,nt) )
    intfP=np.zeros( (nstorm,nt) )
    #dTime = 3. # Hours. SHould make this a calculation not an assumption
    tback = int(dTime/3.)
    
    for s in np.arange( nstorm ):
        for t in np.arange( start=1,stop=nt-tback ):
            wind0=wind[s,t]
            pres0=pres[s,t]
            windD=wind[s,t-tback]
            presD=pres[s,t-tback]
        
            intfW[s,t] = (wind0 - windD)/dTime       
            intfP[s,t] = (pres0 - presD)/dTime       
    
    return intfW,intfP

def density(trx,x,y,category='TD+',**kwargs):
    
    if (category=='TD-only'):
        thresh0=0.
        thresh1=18.
    if (category=='TS-only'):
        thresh0=18.
        thresh1=33.

    if (category=='TD+'):
        thresh0=1.
        thresh1=10000.
    if (category=='TS+'):
        thresh0=18.
        thresh1=10000.
    if (category=='Cat1+'):
        thresh0=33.
        thresh1=10000.
    if (category=='Cat2+'):
        thresh0=43.
        thresh1=10000.
    if (category=='Cat3+'):
        thresh0=50.
        thresh1=10000.
    if (category=='Cat4+'):
        thresh0=58.
        thresh1=10000.
    if (category=='Cat5+'):
        thresh0=70.
        thresh1=10000.
    
    if 'limits' in kwargs:
        limits = kwargs['limits']
    if 'genesis' in kwargs:
        genesis = kwargs['genesis']
    else:
        genesis = False
    
    oo=np.where( trx.year.flatten() > 0)
    yrmin=np.nanmin( trx.year.flatten()[oo[0]] )
    yrmax=np.nanmax( trx.year.flatten()[oo[0]] )
    Nyrs = yrmax-yrmin+1
    
    nx=np.shape( x )[0]
    ny=np.shape( y )[0]
    
    dens = np.zeros( (ny-1, nx-1) , dtype=float )
    
    nstorm,nt=np.shape( trx.wind )
    
    if (genesis==True):
        nt=1
        thresh0=-100
    
    for s in np.arange( nstorm ):
        for t in np.arange( nt ):
            wind0=trx.wind[s,t]
            if ((wind0 >= thresh0) and (wind0 < thresh1)):
                lon0=trx.lon[s,t]
                lat0=trx.lat[s,t]
                if (lon0 < 0):
                    lon0=lon0+360.
                # west edge index
                iw=np.count_nonzero( x<=lon0) - 1
                iw=np.max((0,iw))
                iw=np.min((nx-2,iw))
                # south edge index
                js=np.count_nonzero( y<=lat0) - 1
                js=np.max((0,js))
                js=np.min((ny-2,js))
                dens[ js, iw ] = dens[ js, iw ] + 1.
    
    dens = dens/Nyrs
    
    return dens
            
            
def adensity(trx,x,y,alpha=0,**kwargs):
    
    if ('byYear' in kwargs):
        doYears=kwargs['byYear']
    else:
        doYears=False

    if ('byMonth' in kwargs):
        doMonths=kwargs['byMonth']
    else:
        doMonths=False

    oo=np.where( trx.year.flatten() > 0)
    yrmin=np.nanmin( trx.year.flatten()[oo[0]] )
    yrmax=np.nanmax( trx.year.flatten()[oo[0]] )
    Nyrs = yrmax-yrmin+1
    
    nx=np.shape( x )[0]
    ny=np.shape( y )[0]
    
    nstorm,nt=np.shape( trx.wind )

    if (doYears==True):
        dens = np.zeros( (Nyrs, ny-1, nx-1) , dtype=float )
        for s in np.arange( nstorm ):
            for t in np.arange( nt ):
                wind0=trx.wind[s,t]
                year0=trx.year[s,t]
                if ((wind0 >= 1.) and (year0==Y)):
                    lon0=trx.lon[s,t]
                    lat0=trx.lat[s,t]
                    if (lon0 < 0):
                        lon0=lon0+360.
                    # west edge index
                    iw=np.count_nonzero( x<=lon0) - 1
                    iw=np.max((0,iw))
                    iw=np.min((nx-2,iw))
                    # south edge index
                    js=np.count_nonzero( y<=lat0) - 1
                    js=np.max((0,js))
                    js=np.min((ny-2,js))
                    Y  = int(year0)
                    Yx = Y-yrmin
                    dens[ Yx, js, iw ] = dens[ Yx, js, iw ] + wind0**alpha
        
    elif (doMonths==True):
        dens = np.zeros( (Nyrs, 12, ny-1, nx-1) , dtype=float )
        for s in np.arange( nstorm ):
            for t in np.arange( nt ):
                wind0=trx.wind[s,t]
                year0=trx.year[s,t]
                month0=trx.month[s,t]
                if ((wind0 >= 1.)):
                    lon0=trx.lon[s,t]
                    lat0=trx.lat[s,t]
                    if (lon0 < 0):
                        lon0=lon0+360.
                    # west edge index
                    iw=np.count_nonzero( x<=lon0) - 1
                    iw=np.max((0,iw))
                    iw=np.min((nx-2,iw))
                    # south edge index
                    js=np.count_nonzero( y<=lat0) - 1
                    js=np.max((0,js))
                    js=np.min((ny-2,js))
                    Y  = int(year0)
                    M  = int(month0)
                    Yx = Y-yrmin
                    Mx = M-1
                    dens[ Yx, Mx, js, iw ] = dens[ Yx, Mx, js, iw ] + wind0**alpha

    else:
        dens = np.zeros( (ny-1, nx-1) , dtype=float )
        for s in np.arange( nstorm ):
            for t in np.arange( nt ):
                wind0=trx.wind[s,t]
                if (wind0 >= 1.):
                    lon0=trx.lon[s,t]
                    lat0=trx.lat[s,t]
                    if (lon0 < 0):
                        lon0=lon0+360.
                    # west edge index
                    iw=np.count_nonzero( x<=lon0) - 1
                    iw=np.max((0,iw))
                    iw=np.min((nx-2,iw))
                    # south edge index
                    js=np.count_nonzero( y<=lat0) - 1
                    js=np.max((0,js))
                    js=np.min((ny-2,js))

                    dens[ js, iw ] = dens[ js, iw ] + (wind0**alpha)
                    
    
    
    return dens
            
            
def Prec500( lon0,lat0, x,y,prec ):

    
    Ret=6378.1
    jwin2=20
    iwin2=15

    ny,nx = np.shape( prec )
    if (lon0 < 0):
        lon0=lon0+360.
    
    # center I-index
    ic=np.count_nonzero( x<=lon0) - 1
    # center J-index
    jc=np.count_nonzero( y<=lat0) - 1
    iw,ie = np.max( (0,ic-iwin2) ),np.min( (nx-1,ic+iwin2) )
    js,jn = np.max( (0,jc-jwin2) ),np.min( (ny-1,jc+jwin2) )
    
    precsub = prec[ js:jn, iw:ie ] 

    X=x[iw:ie]
    Y=y[js:jn]
    XX,YY=np.meshgrid(X,Y)
    Xkm=Ret * np.cos( lat0 * np.pi/180.) * (XX-lon0) * np.pi/180.
    Ykm=Ret * (YY-lat0) * np.pi/180.

    weight=1. - np.sqrt( Xkm**2 + Ykm**2 )//500
    precsub=precsub*weight
    prec0 = np.sum( precsub ) / np.sum( weight )

    return prec0    #,precsub,Ykm,Xkm

def Prec500grid( lon0,lat0, x,y,prec ):

    
    Ret=6378.1
    jwin2=20
    iwin2=15

    ny,nx = np.shape( prec )
    if (lon0 < 0):
        lon0=lon0+360.
    
    # center I-index
    ic=np.count_nonzero( x<=lon0) - 1
    # center J-index
    jc=np.count_nonzero( y<=lat0) - 1
    iw,ie = np.max( (0,ic-iwin2) ),np.min( (nx-1,ic+iwin2) )
    js,jn = np.max( (0,jc-jwin2) ),np.min( (ny-1,jc+jwin2) )
    
    precsub = prec[ js:jn, iw:ie ] 

    X=x[iw:ie]
    Y=y[js:jn]
    XX,YY=np.meshgrid(X,Y)
    Xkm=Ret * np.cos( lat0 * np.pi/180.) * (XX-lon0) * np.pi/180.
    Ykm=Ret * (YY-lat0) * np.pi/180.

    weight=1. - np.sqrt( Xkm**2 + Ykm**2 )//500
    precsub=precsub*weight

    return precsub,js,jn,iw,ie  

def Heal_rc4_rc6( rc4, rc6 ):
    ## Heal rc4 and rc6
    for s in np.arange( start=826,stop=864 ):
        oo=np.where( rc4.year[s,:] > 0)
        rc4.year[s,oo[0]]=2083
    for s in np.arange( start=1158,stop=1175 ):
        oo=np.where( rc4.year[s,:] > 0)
        rc4.year[s,oo[0]]=2088
    for s in np.arange( start=1654,stop=1676 ):
        oo=np.where( rc4.year[s,:] > 0)
        rc4.year[s,oo[0]]=2096
    for s in np.arange( start=1266,stop=1288 ):
        oo=np.where( rc6.year[s,:] > 0)
        rc6.year[s,oo[0]]=2089

    return rc4,rc6

def Heal_rc4( rc4 ):
    ## Heal rc4
    for s in np.arange( start=826,stop=864 ):
        oo=np.where( rc4.year[s,:] > 0)
        rc4.year[s,oo[0]]=2083
    for s in np.arange( start=1158,stop=1175 ):
        oo=np.where( rc4.year[s,:] > 0)
        rc4.year[s,oo[0]]=2088
    for s in np.arange( start=1654,stop=1676 ):
        oo=np.where( rc4.year[s,:] > 0)
        rc4.year[s,oo[0]]=2096

    return rc4

def Heal_rc6( rc6 ):
    ## Heal rc4 and rc6
    for s in np.arange( start=1266,stop=1288 ):
        oo=np.where( rc6.year[s,:] > 0)
        rc6.year[s,oo[0]]=2089

    return rc6

def HealZeroWind_rc4( rc4 ):
    ## Heal rc4
    for s in np.arange( start=826,stop=864 ):
        oo=np.where( rc4.year[s,:] > 0)
        rc4.wind[s,oo[0]]=0.
        rc4.year[s,oo[0]]=-99999.
    for s in np.arange( start=1158,stop=1175 ):
        oo=np.where( rc4.year[s,:] > 0)
        rc4.wind[s,oo[0]]=0.
        rc4.year[s,oo[0]]=-99999.
    for s in np.arange( start=1654,stop=1676 ):
        oo=np.where( rc4.year[s,:] > 0)
        rc4.wind[s,oo[0]]=0.
        rc4.year[s,oo[0]]=-99999.

    return rc4

def HealZeroWind_rc6( rc6 ):
    ## Heal rc4 and rc6
    for s in np.arange( start=1266,stop=1288 ):
        oo=np.where( rc6.year[s,:] > 0)
        rc6.wind[s,oo[0]]=0.
        rc6.year[s,oo[0]]=-99999.

    return rc6


def is_left(P0, P1, P2):
    """
    Test if point P2 lies to the left of the line formed by P0 and P1.
    >0 for P2 left of the line through P0 and P1
    =0 for P2 on the line
    <0 for P2 right of the line
    """
    return (P1[0] - P0[0]) * (P2[1] - P0[1]) - (P2[0] - P0[0]) * (P1[1] - P0[1])

def winding_number(x, y, X, Y):
    """
    Determine if the point (x, y) is inside the polygon defined by X and Y.
    :param x, y: Coordinates of the point to test
    :param X, Y: Lists of x and y coordinates of polygon's vertices
    :return: Boolean. True if point is inside polygon, else False
    """
    wn = 0  # Initialize the winding number to be 0

    # Loop through all edges of the polygon
    for i in range(len(X)):
        # If the point is exactly on a vertex, consider it inside
        if x == X[i] and y == Y[i]:
            return True

        # Get start and end points of the edge
        x1, y1 = X[i], Y[i]
        x2, y2 = X[(i + 1) % len(X)], Y[(i + 1) % len(X)]

        if y1 <= y:
            if y2 > y:
                # An upward crossing
                if is_left((x1, y1), (x2, y2), (x, y)) > 0:
                    # Point is to the left of edge, so increment winding number
                    wn += 1
        else:
            if y2 <= y:
                # A downward crossing
                if is_left((x1, y1), (x2, y2), (x, y)) < 0:
                    # Point is to the right of edge, so decrement winding number
                    wn -= 1

    # Point is inside the polygon if the winding number is nonzero
    return wn != 0
"""
# Example usage
X = [0, 2, 2, 0]
Y = [0, 0, 2, 2]
print(winding_number(1, 1, X, Y))  # Should return True (inside)
print(winding_number(3, 3, X, Y))  # Should return False (outside)
"""
# This allow both dict.key and dict['key'] syntax
class AttrDict(dict):
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(f"'AttrDict' object has no attribute '{key}'")

    def __setattr__(self, key, value):
        self[key] = value

    def __delattr__(self, key):
        try:
            del self[key]
        except KeyError:
            raise AttributeError(f"'AttrDict' object has no attribute '{key}'")

import sys
import os
This_module_path = os.getcwd()  #os.path.dirname(os.path.abspath(__file__))
sys.path.append(This_module_path )

import constants as co

def IwIe( i , nx, wrap=True ):
    iw = i - 1
    if (i==0):
        if (wrap==True):
            iw=nx-1
        else:
            iw=0
    ie = i + 1
    if (i==nx-1):
        if (wrap==True):
            ie=0
        else:
            ie=nx-1
    return iw,ie

def JsJn( j , ny ):
    js = j - 1
    if (j==0):
           js=0
    jn = j + 1
    if (j==ny-1):
        jn=ny-1
    return js,jn


def Sphere_Grad2( f, lat, lon, wrap=True ):
    import numpy as np

    Pi  = co.pi()
    R_e = co.Rearth()

    nx,ny = len( lon ), len( lat )
    
    rlat = (Pi/180.)*lat
    rlon = (Pi/180.)*lon
    dlat=rlat[2]-rlat[0]
    dlon=rlon[2]-rlon[0]

    f_x = np.zeros( (ny,nx) )
    for j in np.arange( ny ):
        coslat = np.cos( rlat[j] )
        for i in np.arange( nx ):
            iw,ie = IwIe(i,nx,wrap=wrap)
            f_x[j,i] = (1./(R_e * coslat))*( f[j,ie]-f[j,iw] ) / dlon
    f_y = np.zeros( (ny,nx) )
    for j in np.arange( ny ):
        for i in np.arange( nx ):
            js,jn = JsJn(j,ny)
            f_y[j,i] = (1./ R_e)*( f[jn,i]-f[js,i] ) / dlat

    return f_x, f_y

def Sphere_Div2( f_x, f_y, lat, lon, wrap=True ):
    import numpy as np

    Pi  = co.pi()
    R_e = co.Rearth()

    nx,ny = len( lon ), len( lat )
    
    rlat = (Pi/180.)*lat
    rlon = (Pi/180.)*lon
    dlat=rlat[2]-rlat[0]
    dlon=rlon[2]-rlon[0]
    coslat = np.cos( rlat )

    divf=np.zeros( (ny, nx ))
    for j in np.arange( ny ):
        for i in np.arange( nx ):
            iw,ie = IwIe(i,nx,wrap=wrap)
            js,jn = JsJn(j,ny)
            divf[j,i] = (1./(R_e * coslat[j]))*( 
                        ( f_x[j,ie]-f_x[j,iw] ) / dlon
                    +   ( coslat[jn]*f_y[jn,i] - coslat[js]*f_y[js,i] )/dlat )

    return divf


def Sphere_Curl2( f_x, f_y, lat, lon, wrap=True ):
    ##################################################
    # Curl of a vector field f=[f_x,f_y] in latlon
    # coords:
    #        Zeta ~ d(f_y)/dx - d(f_x)/dy 
    ##################################################
    import numpy as np

    Pi  = co.pi()
    R_e = co.Rearth()

    nx,ny = len( lon ), len( lat )
    
    rlat = (Pi/180.)*lat
    rlon = (Pi/180.)*lon
    dlat=rlat[2]-rlat[0]
    dlon=rlon[2]-rlon[0]
    coslat = np.cos( rlat )

    curlf_z = np.zeros( (ny, nx ))
    for j in np.arange( ny ):
        for i in np.arange( nx ):
            iw,ie = IwIe(i,nx,wrap=wrap)
            js,jn = JsJn(j,ny)
            curlf_z[j,i] = (1./(R_e * coslat[j]))*( 
                     1.0* ( f_y[j,ie]-f_y[j,iw] ) / dlon
                    -   ( coslat[jn]*f_x[jn,i] - coslat[js]*f_x[js,i] )/dlat )

    return curlf_z



def Sphere_Lap2( f , lat, lon, wrap=True ):

    f_x,f_y = Sphere_Grad2( f=f, lon=lon,lat=lat,wrap=wrap )
    lapf    = Sphere_Div2( f_x=f_x, f_y=f_y, lon=lon,lat=lat,wrap=wrap )
    
    return lapf




    
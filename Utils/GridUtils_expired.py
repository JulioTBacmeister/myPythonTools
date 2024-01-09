import numpy as np
import MyConstants as C

pi = C.pi()


def area2d(lon,lat):
    #inputs 
    #   lon 1D lon vector
    #   lat 1D lat vector
    
    nx = np.size( lon )
    ny = np.size( lat )
    area = np.zeros( (ny, nx) )
    latr = (pi/180.) * lat
    
    for j in np.arange( ny ):
        area[j,:] = np.cos(latr[j])
        
    area = ( 4*pi / np.sum(area) ) * area
    
    return area
    
def gridKey( Var ):
    
    VarDims = Var.dims
    ndim = len(VarDims)
    
    print( "why am I here " )
    
    if (VarDims[0]=='time'):
        gridKey = 't'
    elif (VarDims[0]=='lev'):
        gridKey = 'z'
    elif (VarDims[0]=='lat'):
        gridKey = 'y'
    elif (VarDims[0]=='ncol'):
        gridKey = 'c'
    else:
        gridKey = 'not found'
        return gridKey
    
    if (ndim>1):
        if (VarDims[1]=='time'):
            gridKey = gridKey + 't'
        elif (VarDims[1]=='lev'):
            gridKey = gridKey + 'z'
        elif (VarDims[1]=='lat'):
            gridKey = gridKey + 'y'
        elif (VarDims[1]=='ncol'):
            gridKey = gridKey + 'c'
        else:
            gridKey = gridKey + ' not found'
            return gridKey
    
    if (ndim>2):
        if (VarDims[2]=='time'):
            gridKey = gridKey + 't'
        elif (VarDims[2]=='lev'):
            gridKey = gridKey + 'z'
        elif (VarDims[2]=='lat'):
            gridKey = gridKey + 'y'
        elif (VarDims[2]=='ncol'):
            gridKey = gridKey + 'c'
        else:
            gridKey = gridKey + ' not found'
            return gridKey
    
    if (ndim>3):
        if (VarDims[3]=='time'):
            gridKey = gridKey + 't'
        elif (VarDims[3]=='lev'):
            gridKey = gridKey + 'z'
        elif (VarDims[3]=='lat'):
            gridKey = gridKey + 'y'
        elif (VarDims[3]=='lon'):
            gridKey = gridKey + 'x'
        elif (VarDims[3]=='ncol'):
            gridKey = gridKey + 'c'
        else:
            gridKey = gridKey + ' not found'
            return gridKey
    
    
    return gridKey
#workdir_ = '/glade/work/juliob/'
import sys
#######################################
# Leave this for now. But it should change to better
# method as here:
import os
This_module_path = os.path.dirname(os.path.abspath(__file__))
workdir_ = os.path.join(This_module_path, '../../' )
# sys.path.append(utils_path)
# print( f" a path added in {__name__} {utils_path} ")

print( f" In {__name__} we have This_module_path={This_module_path} " )
print( f" In {__name__} we have workdir_={workdir_} " )
########################################
sys.path.append(workdir_ + 'myPythonTools/GenlTools/')
sys.path.append(workdir_ + 'myPythonTools/Utils/')
sys.path.append(workdir_ + 'PyRegridding/Regridder/')
sys.path.append(workdir_ + 'PyRegridding/Utils/')

# Own local packages
import AveragingUtils as Av
#import VertRegridFlexLL as Vrg  # This is toxic for some reason
import PlotUtil as Pu
import utils as uti
import validation_data as Val
import var_A_x_B as vAB
import MakePressures as MkP
import GridUtils as GU
import esmfRegrid as erg

# The usual
from datetime import date
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.tri as tri

# Cartopy for pretty maps
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Some other useful packages 
import math
import importlib
import copy
import time
import cftime

importlib.reload( uti )
importlib.reload( Pu )
importlib.reload(Av)
importlib.reload(Val)
importlib.reload(vAB)
importlib.reload(MkP)


def Maps( fields , lons, lats, **kwargs ):
    ##################################################
    #           "Pretty" latlon plots
    #---------------------------------------------
    # Only mandatory input should be a fields list
    # containing 2D (lat,lon) variables to be plotted, 
    # along with lons and lats.  Currently, lons and lats
    # must also be list with same length as fields
    ##################################################
    
    MapProj = ccrs.PlateCarree(central_longitude=180.)
    DataProj = ccrs.PlateCarree()

    nplots = len( fields )
    
    fig = plt.figure(figsize=(20, 12))
    
    if ( 'CoastColors' in kwargs ):
        CoastColors = kwargs[ 'CoastColors' ]
    else:
        value = 'black'
        CoastColors = [value for _ in range(nplots)]
    if ( 'clevs' in kwargs ):
        clevs = kwargs[ 'clevs' ]
    else:
        value = 21
        clevs = [value for _ in range(nplots)]
    if ( 'cmaps' in kwargs ):
        cmaps = kwargs[ 'cmaps' ]
    else:
        value = 'gist_ncar'
        cmaps = [value for _ in range(nplots)]
    if ( 'scale' in kwargs ):
        scale = kwargs[ 'scale' ]
    else:
        value = 1.0
        scale = [value for _ in range(nplots)]
    

    if (nplots==1):
        nx,ny=1,1
    elif (nplots==2):
        nx,ny=2,1
    elif (nplots==3):
        nx,ny=3,1
    elif (nplots==4):
        nx,ny=2,2
    elif ((nplots>4)and(nplots<=9)):
        nx,ny=3,3
    else:
        ny = 3
        nx = math.ceil( nplots / ny)

    
    npo=0
    for field in fields:
        ipo=npo
        npo=npo+1
        Axes = Pu.axes_def(n=npo,nxplo=nx,nyplo=ny ) 
    
        ax1 = fig.add_axes( Axes , projection=MapProj)
        ax1.set_global()
        ax1.coastlines(resolution='110m',color='black',linewidth=2)
    
        AAxy = scale[ipo]*field 
        lon_ = lons[ipo]
        lat_ = lats[ipo]
        
        co1=ax1.contourf(lon_ ,lat_ , AAxy ,transform=DataProj,cmap=cmaps[ipo] ,levels=clevs[ipo])
    
        cbar = plt.colorbar(co1, ax=ax1, shrink=.6)
        #ax1.set_title( CLUBBparm + zA, fontsize=16)
        #ax1.set_title( f"{season.upper()} {flds[ipo]} {exps[ipo]}" , fontsize=16)


def Regrid( fields , lon, lat, **kwargs ):

    try:
        import ESMF as E
    except ImportError:
        import esmpy as E

    
    # Initialize ESMF
    E.Manager()
    
    # Define the unstructured source grid
    # Assuming you have lists or arrays of latitudes and longitudes for your unstructured grid
    src_lons = lon 
    src_lats = lat 
    num_nodes = len(src_lons)
    
    """
    # Create an E.Mesh for the unstructured grid
    src_mesh = E.Mesh(parametric_dim=2, spatial_dim=2)

    src_mesh.add_nodes(np.arange(num_nodes), np.stack((src_lons, src_lats), axis=1))
    """


    # Create an E.Mesh for the unstructured grid
    src_mesh = E.Mesh(parametric_dim=2, spatial_dim=2)
    
    # Define node coordinates and node owners
    #node_coords = np.stack((src_lons, src_lats), axis=0)
    node_coords = np.concatenate((src_lats, src_lons)).reshape( 2*num_nodes, 1 )
    node_owners = np.zeros(num_nodes, dtype=int)
    node_ids = np.arange(num_nodes)  # Unique identifiers for each node

    print( "Node coords ",np.shape(node_coords))

    
    # Add nodes to the mesh
    #src_mesh.add_nodes(np.arange(num_nodes), node_coords = node_coords , node_ids=node_ids ) #, node_owners)
    src_mesh.add_nodes(node_count=num_nodes, node_coords = node_coords , node_ids=node_ids , node_owners=node_owners)

    
    
    # Define the destination grid (structured lat-lon grid as before)
    dst_lon = np.linspace(0, 360, 360, endpoint=False)
    dst_lat = np.linspace(-80, 80, 181, endpoint=True)
    dst_grid = E.Grid(np.array([len(dst_lon), len(dst_lat)]), staggerloc=[E.StaggerLoc.CENTER])
    dst_lon_par, dst_lat_par = [dst_grid.get_coords(i) for i in range(2)]
    dst_lon_par[...] = dst_lon.reshape((dst_lon.shape[0], 1))
    dst_lat_par[...] = dst_lat.reshape((1, dst_lat.shape[0]))
    
    # Create field objects on the source mesh and destination grid
    src_field = E.Field(src_mesh, name='Source Data', meshloc=E.MeshLoc.NODE ) # meshloc=E.MeshLoc.NODE)
    dst_field = E.Field(dst_grid, name='Destination Data')
    
    # Initialize the source field with some data, for example, a simple pattern
    # This will depend on how your data is structured; here's a placeholder:
    print( np.shape(src_field.data ))
    src_field.data[...] = field[...]
    
    # Create a regridder object from source to destination
    regridder = E.Regrid(src_field, dst_field, regrid_method=E.RegridMethod.BILINEAR,
                            unmapped_action=E.UnmappedAction.IGNORE)
    
    # Perform the regridding
    dst_field = regridder(src_field, dst_field)
    
    # dst_field.data now contains the regridded data

    return dst_field.data, dst_lat, dst_lon
        
def ScripRegrid( fields , Src, Dst, RegridMethod = 'BILINEAR' ):

    if ( isinstance(fields, np.ndarray) == True ):
        fields = [fields]
        InputIsNdarray = True
    else:
        InputIsNdarray = False

    
    
    Src_gridInfo = GU.gridInfo( grid=Src )
    src_scrip, src_Hkey, src_type = Src_gridInfo['scrip'],Src_gridInfo['Hkey'],Src_gridInfo['type']
    
    Dst_gridInfo = GU.gridInfo( grid=Dst )
    dst_scrip, dst_Hkey, dst_type = Dst_gridInfo['scrip'],Dst_gridInfo['Hkey'],Dst_gridInfo['type']

    print(f" Src grid {Src} {src_scrip} , {src_Hkey} , {src_type}")
    print(f" Dst grid {Dst} {dst_scrip} , {dst_Hkey} , {dst_type}")

    lat_dst,lon_dst = GU.latlon( scrip=dst_scrip , Hkey = dst_Hkey )

    regrd, srcf, dstf = erg.Regrid( srcScrip = src_scrip , 
                                    srcType  = src_type  ,
                                    dstScrip = dst_scrip ,
                                    dstType  = dst_type  ,
                                    RegridMethod = RegridMethod )

    dst_fields = []
    for field in fields:
        field_dst = erg.HorzRG( aSrc = field , 
                     regrd = regrd , 
                     srcField= srcf , 
                     dstField= dstf , 
                     srcGridkey= src_Hkey ,
                     dstGridkey= dst_Hkey )
        dst_fields.append( field_dst )

    if InputIsNdarray == True:
        dst_fields=np.array(dst_fields[0] ) 
        
    return dst_fields, lat_dst, lon_dst
        

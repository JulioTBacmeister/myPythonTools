import sys
import os

workdir_ = '/glade/work/juliob'
if ( workdir_ not in sys.path ):
    sys.path.append(workdir_)
    print( f" a path to {workdir_} added in {__name__} ")



# Own local packages
from myPythonTools.Utils import AveragingUtils as Av
from myPythonTools.Utils import PlotUtil as Pu
from myPythonTools.Utils import utils as uti
#from myPythonTools.Utils import validation_data as Val
from PyRegridding.Utils import MakePressures as MkP
from PyRegridding.Utils import GridUtils as GU
from PyRegridding.Regridder import esmfRegrid as erg

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
#importlib.reload(Val)
importlib.reload(MkP)

def titleGen( exp,fld,season,years ):
    title=rf"{fld} <{exp}> {season.upper()} Years:{str(years[0]).zfill(4)}-{str(years[-1]).zfill(4)}" 
    return title



def Maps( fields , lons, lats, overplot=None, **kwargs ):
    ##################################################
    #           "Pretty" latlon plots
    #---------------------------------------------
    # Only mandatory input should be a fields list
    # containing 2D (lat,lon) variables to be plotted, 
    # along with lons and lats.  Currently, lons and lats
    # must also be list with same length as fields
    ##################################################
    

    
    print(" LL top 0.0.0 " )
    nplots = len( fields )
    
    # Set up defaults for kwargs that are needed:
    if ( 'Projection' in kwargs ):
        MapProj = kwargs[ 'Projection' ]
    else:
        MapProj = ccrs.PlateCarree(central_longitude=180.)

    DataProj = ccrs.PlateCarree()
    # Get the name of the projection
    proj_name = MapProj.__class__.__name__
    print(" LL 1.0.0 " )

    
    if ( 'map_extent' in kwargs ):
        map_extent_ = kwargs[ 'map_extent' ]
    else:
        map_extent_ = None

    if ( 'verbose' in kwargs ):
        verbose_ = kwargs[ 'verbose' ]
    else:
        verbose_ = False
    
    if ( 'CoastColors' in kwargs ):
        CoastColors = kwargs[ 'CoastColors' ]
    else:
        value = 'black'
        CoastColors = [value for _ in range(nplots)]

    if ( 'CoastWidth' in kwargs ):
        CoastWidth = kwargs[ 'CoastWidth' ]
    else:
        value = 2
        CoastWidth = [value for _ in range(nplots)]

    if ( 'coastline_res' in kwargs ):
        if ( kwargs[ 'coastline_res' ].lower()=='low'):
            coastline_res = '110m'
        elif ( kwargs[ 'coastline_res' ].lower()=='medium'):
            coastline_res = '50m'
        elif ( kwargs[ 'coastline_res' ].lower()=='high'):
            coastline_res = '10m'
        else:
            print( "No valid input for coastline res. Using 110m " )
            coastline_res = '110m'
    else:
        coastline_res = '110m'

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

    if ( 'titles' in kwargs ):
        titles = kwargs[ 'titles' ]
    else:
        value = 'Title Placeholder'
        titles = [value for _ in range(nplots)]

    if ( 'cbarlabels' in kwargs ):
        cbarlabels = kwargs[ 'cbarlabels' ]
        #print( "Set cbarlabels ")
    else:
        value = '.'
        cbarlabels = [value for _ in range(nplots)]

    if ( 'insets' in kwargs ):
        insets = kwargs[ 'insets' ]
    else:
        insets = None 

    if ( 'Landmask' in kwargs ):
        add_landmask = kwargs[ 'Landmask' ]
    else:
        add_landmask = False 
    
    if ( 'fig_xsize_per_panel' in kwargs ):
        xsize = kwargs[ 'fig_xsize_per_panel' ]
    else:
        xsize = 10.
    
    print( f"xsize per panel ={xsize}" )
    
    if ( 'Arrangement' in kwargs ):
        nx,ny = kwargs['Arrangement']
        print( f"Overiding arrange with {nx}x{ny} " )
    else:
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

    print(" LL 2.0.0 " )

    print(f'proj_name = {proj_name}')
    if (proj_name=='PlateCarree'):
        ysize=xsize*0.5
    else:
        ysize=xsize
    fig = plt.figure(figsize=( nx*xsize, ny*ysize ))
    
    
    npo=0
    for field in fields:
        ipo=npo
        npo=npo+1
        Axes = Pu.axes_def(n=npo,nxplo=nx,nyplo=ny ) 
    
        ax1 = fig.add_axes( Axes , projection=MapProj)

        if (map_extent_ is None):
            ax1.set_global()
        else:
            # Set extent: [lon_min, lon_max, lat_min, lat_max]
            ax1.set_extent( map_extent_ , crs=DataProj)

        ax1.coastlines(resolution=coastline_res,color=CoastColors,linewidth=CoastWidth )

        # Add state boundaries
        ax1.add_feature(cfeature.STATES,
                       edgecolor='black',
                       linewidth=0.6)

        
        print(f" LL 1.0.{ipo} " )
    
        AAxy = scale[ipo]*field 
        if ( len(lons) >= len(fields) ):
            lon_ = lons[ipo]
            lat_ = lats[ipo]
        else:
            lon_ = lons[0]
            lat_ = lats[0]
        if ( len(clevs) >= len(fields) ):
            clev_ = clevs[ipo]
        else:
            clev_ = clevs[0]
        if ( len(cmaps) >= len(fields) ):
            cmap_ = cmaps[ipo]
        else:
            cmap_ = cmaps[0]
        if ( len(titles) >= len(fields) ):
            title_ = titles[ipo]
        else:
            title_ = titles[0]
        if ( len(cbarlabels) >= len(fields) ):
            cbarlabel_ = cbarlabels[ipo]
        else:
            cbarlabel_ = None
        if insets is not None:
            if ( ipo < len(insets)  ):
                inset=insets[ipo]
            else:
                inset = None
        else:
            inset = None

        if( verbose_==True ):
            print( f' color map={cmap_} clevels={clev_} ' )

        co1=ax1.contourf(lon_ ,lat_ , AAxy ,transform=DataProj,cmap=cmap_ ,levels=clev_, extend='both' )

        if (overplot is not None):        
            ovfields,lat2,lon2 = overplot['field'],overplot['lat'],overplot['lon']

            ovfield = ovfields[ipo]
            
            if ('clevels' in overplot):
                clevels = overplot['clevels']
            else:
                clevels=21
            if ('alpha' in overplot):
                alpha = overplot['alpha']
            else:
                alpha = 0.5
            if ('cmap' in overplot):
                ovcmap = overplot['cmap']
            else:
                ovcmap = 'Greys'
                
            co2 = ax1.contourf( lon2 ,lat2 , ovfield ,transform=DataProj,cmap=ovcmap ,levels=clevels, extend='both' ,alpha=alpha)
            li2 = ax1.contour( lon2 ,lat2 , ovfield ,transform=DataProj, cmap=ovcmap ,levels=clevels, extend='both' ,alpha=0.8 )
            
        
        gls=ax1.gridlines(draw_labels=True)
    
        cbar = plt.colorbar(co1, ax=ax1, fraction=0.046, pad=0.04 ,aspect=10)

        if cbarlabel_ is not None:
            cbar.set_label( cbarlabel_ ,fontsize=18 )

        ax1.set_title( title_ , fontsize=16)
        #ax1.set_title( f"{season.upper()} {flds[ipo]} {exps[ipo]}" , fontsize=16)
        if (add_landmask==True):
            print('Adding Landmask' )
            ax1.add_feature(cfeature.LAND, facecolor='tan', zorder=10) 
        else:
            print('no Landmask')


        ax1.text(-0.08, 1.05, f"{chr(97 +npo-1+1)})", transform=ax1.transAxes,
            fontsize=16, fontweight='bold', va='top')

        if inset is not None:
            ax1.text(
                inset['x'], inset['y'],  # X, Y position in axes coordinates
                inset['text'],  # Text content
                transform=ax1.transAxes,  # Use axes coordinates
                fontsize=10,
                verticalalignment='top',
                bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='lightgrey')
                )
        #ax1.set_xlims(50,100)
        #ax1.set_ylims(15,50)

def Maps_NoProj( fields , lons, lats, **kwargs ):
    ##################################################
    #           "Simple" latlon plots for sanity
    #---------------------------------------------
    # Only mandatory input should be a fields list
    # containing 2D (lat,lon) variables to be plotted, 
    # along with lons and lats.  Currently, lons and lats
    # must also be list with same length as fields
    ##################################################
    

    nplots = len( fields )
    
    fig = plt.figure(figsize=(20, 12))
    #fig = plt.figure(figsize=(30, 18))
    
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

    if ( 'titles' in kwargs ):
        titles = kwargs[ 'titles' ]
    else:
        value = 'Title Placeholder'
        titles = [value for _ in range(nplots)]
    

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
    
        ax1 = fig.add_axes( Axes ) #, projection=MapProj)
    
        AAxy = scale[ipo]*field 
        if ( len(lons) >= len(fields) ):
            lon_ = lons[ipo]
            lat_ = lats[ipo]
        else:
            lon_ = lons[0]
            lat_ = lats[0]
        if ( len(clevs) >= len(fields) ):
            clev_ = clevs[ipo]
        else:
            clev_ = clevs[0]
        if ( len(cmaps) >= len(fields) ):
            cmap_ = cmaps[ipo]
        else:
            cmap_ = cmaps[0]
        if ( len(titles) >= len(fields) ):
            title_ = titles[ipo]
        else:
            title_ = titles[0]

        co1=ax1.contourf(lon_ ,lat_ , AAxy , cmap=cmap_ ,levels=clev_ , extend='both' )
    
        cbar = plt.colorbar(co1, ax=ax1, fraction=0.046, pad=0.04 ,aspect=10)
        ax1.set_title( title_ , fontsize=16)
        #ax1.set_title( f"{season.upper()} {flds[ipo]} {exps[ipo]}" , fontsize=16)

def Maps_Unstruc( fields , lons, lats, **kwargs ):
    ##################################################
    #           "Simple" latlon plots for sanity
    #---------------------------------------------
    # Only mandatory input should be a fields list
    # containing 2D (lat,lon) variables to be plotted, 
    # along with lons and lats.  Currently, lons and lats
    # must also be list with same length as fields
    ##################################################
    

    nplots = len( fields )
    
    #fig = plt.figure(figsize=(30, 18))
    
    if ( 'figsize' in kwargs ):
        figsize = kwargs[ 'figsize' ]
    else:
        figsize=(20,12)
    
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

    if ( 'titles' in kwargs ):
        titles = kwargs[ 'titles' ]
    else:
        value = 'Title Placeholder'
        titles = [value for _ in range(nplots)]
    

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

    fig = plt.figure(figsize=figsize )
    
    npo=0
    for field in fields:
        ipo=npo
        npo=npo+1
        Axes = Pu.axes_def(n=npo,nxplo=nx,nyplo=ny ) 
    
        ax1 = fig.add_axes( Axes ) #, projection=MapProj)
    
        AAxy = scale[ipo]*field 
        if ( len(lons) >= len(fields) ):
            lon_ = lons[ipo]
            lat_ = lats[ipo]
        else:
            lon_ = lons[0]
            lat_ = lats[0]
        if ( len(clevs) >= len(fields) ):
            clev_ = clevs[ipo]
        else:
            clev_ = clevs[0]
        if ( len(cmaps) >= len(fields) ):
            cmap_ = cmaps[ipo]
        else:
            cmap_ = cmaps[0]
        if ( len(titles) >= len(fields) ):
            title_ = titles[ipo]
        else:
            title_ = titles[0]

        co1=ax1.tricontourf(lon_ ,lat_ , AAxy , cmap=cmap_ ,levels=clev_  , extend='both' )
        if ( 'topo' in kwargs ):
            topo=kwargs['topo']
            ltop=ax1.tricontour(lon_ ,lat_ , topo , colors='white' ,levels= [5,1000,2000,5000] ) #[1,10,100,1000] )
            
    
        cbar = plt.colorbar(co1, ax=ax1, fraction=0.046, pad=0.04 ,aspect=10)
        ax1.set_title( title_ , fontsize=16)
        #ax1.set_title( f"{season.upper()} {flds[ipo]} {exps[ipo]}" , fontsize=16)

def ZonalMeans( fields , lats, **kwargs ):
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
    # Get the name of the projection
    proj_name = MapProj.__class__.__name__

    nplots = len( fields )
    
    # Set up defaults for kwargs that are needed:
    if ( 'verbose' in kwargs ):
        verbose_ = kwargs[ 'verbose' ]
    else:
        verbose_ = False
    
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

    if ( 'titles' in kwargs ):
        titles = kwargs[ 'titles' ]
    else:
        value = 'Title Placeholder'
        titles = [value for _ in range(nplots)]
    

    if ( 'Arrangement' in kwargs ):
        nx,ny = kwargs['Arrangement']
    else:
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

    xsize=10.
    print(f'proj_name = {proj_name}')
    if (proj_name=='PlateCarree'):
        ysize=xsize*0.5
    else:
        ysize=xsize
    fig = plt.figure(figsize=( nx*xsize, ny*ysize ))
    
    
    npo=0
    for field in fields:
        ipo=npo
        npo=npo+1
        Axes = Pu.axes_def(n=npo,nxplo=nx,nyplo=ny ) 
    
        ax1 = fig.add_axes( Axes , projection=MapProj)
        ax1.set_global()
        ax1.coastlines(resolution='110m',color='black',linewidth=2)
    
        AAxy = scale[ipo]*field 
        if ( len(lons) >= len(fields) ):
            lon_ = lons[ipo]
            lat_ = lats[ipo]
        else:
            lon_ = lons[0]
            lat_ = lats[0]
        if ( len(clevs) >= len(fields) ):
            clev_ = clevs[ipo]
        else:
            clev_ = clevs[0]
        if ( len(cmaps) >= len(fields) ):
            cmap_ = cmaps[ipo]
        else:
            cmap_ = cmaps[0]
        if ( len(titles) >= len(fields) ):
            title_ = titles[ipo]
        else:
            title_ = titles[0]

        if( verbose_==True ):
            print( f' color map={cmap_} clevels={clev_} ' )

        co1=ax1.contourf(lon_ ,lat_ , AAxy ,transform=DataProj,cmap=cmap_ ,levels=clev_, extend='both' )
    
        cbar = plt.colorbar(co1, ax=ax1, fraction=0.046, pad=0.04 ,aspect=10)
        ax1.set_title( title_ , fontsize=16)
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
        

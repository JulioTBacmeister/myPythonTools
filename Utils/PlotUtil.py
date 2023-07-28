import numpy as np


def axes_def(n,nxplo=1,nyplo=1,xpan=0.9,ypan=0.9,xcor=1.,ycor=1.):
    dxplo=xpan/nxplo 
    dyplo=ypan/nyplo
    dxcor=xcor/nxplo
    dycor=ycor/nyplo

    jplo = (n-1)//nxplo #+ 1
    iplo = n-1 - jplo*nxplo
    Axes = [ iplo*dxcor , 1.-(jplo+1)*dycor , dxplo , dyplo ] 
    
    return Axes

def right_cbar_axes(Axes,xshift=0,yshift=0, xsize=0.01,ysize=0.3 ):
    CrAxes = [ Axes[0]+Axes[2]+xshift ,Axes[1]+yshift , xsize , ysize  ] 
    return CrAxes

def step_cmap(cmapN,levels):
    import matplotlib.colors as colors
    import matplotlib.pyplot as plt


    # Normalized contour levels
    norm_levels = (levels - levels.min()) / (levels.max() - levels.min())
    # Creating a colormap from the "viridis" colormap
    colomap = plt.get_cmap(cmapN)
    colors_array = colomap(np.linspace(0, 1, len(norm_levels)))
    # Creating a new colormap
    cmap = colors.LinearSegmentedColormap.from_list("nonuniform", list(zip(norm_levels, colors_array)))

    return cmap
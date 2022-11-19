from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay

def examplot(npoints=100):

    rng = np.random.default_rng()
    x = rng.random(npoints) - 0.5
    y = rng.random(npoints) - 0.5
    z = np.hypot(x, y)
    z = rng.random(npoints) 
    X = np.linspace(min(x), max(x), num=100) # linspace defaults to 50 samples
    Y = np.linspace(min(y), max(y), num=100)

    XX, YY = np.meshgrid(X, Y)
    
    ooo=np.transpose( np.asarray([x,y]) )
    poo=Delaunay( ooo )
    interp = LinearNDInterpolator( poo , z)
    Z1=interp(XX,YY)

    clevels=np.linspace(0.,1.,num=20)
    
    interp2 = LinearNDInterpolator( list(zip(x,y)), z)
    Z2=interp2(XX,YY)
    
    ZOG = np.hypot(XX, YY)

    plt.ion()

    plt.figure(1).clear
    plt.pcolormesh(XX, YY, ZOG ,vmin=0, vmax=1.)
    plt.colorbar()
    plt.contour(XX, YY, ZOG, levels=clevels, colors='black' )
    #plt.plot(x, y, "ok", label="input point")
    #plt.legend()
    #plt.title(" zip(x,y) input " )
    plt.title(" Gridded input " )
    plt.axis("equal") # gives aspect ratio=1

    plt.figure(2)

    plt.pcolormesh(XX, YY, Z1 ,vmin=0, vmax=1. )
    plt.colorbar()
    #plt.contour(XX, YY, Z1, levels=clevels, colors='gray' )
    plt.plot(x, y, "ok", label="input point")
    #plt.legend()
    plt.title(" Precomputed Delaunay input " )
    plt.axis("equal")

    diff=Z1-Z2
    #print(diff[:,24])

    plt.figure(3)

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    ax.plot_surface(XX, YY, Z1 ,vmin=0, vmax=1. )

    
    return(Z1,XX,YY)

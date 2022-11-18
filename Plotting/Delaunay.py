from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay

def examplot():

    rng = np.random.default_rng()
    x = rng.random(10) - 0.5
    y = rng.random(10) - 0.5
    z = np.hypot(x, y)
    X = np.linspace(min(x), max(x)) # linspace defaults to 50 samples
    Y = np.linspace(min(y), max(y))

    XX, YY = np.meshgrid(X, Y)
    
    ooo=np.transpose( np.asarray([x,y]) )
    poo=Delaunay( ooo )
    interp = LinearNDInterpolator( poo , z)
    Z1=interp(XX,YY)
    
    
    interp2 = LinearNDInterpolator( list(zip(x,y)), z)
    Z2=interp2(XX,YY)
    
    
    plt.ion()

    plt.pcolormesh(XX, YY, Z2 )
    plt.colorbar()
    plt.plot(x, y, "ok", label="input point")
    plt.legend()
    plt.title(" zip(x,y) input " )
    plt.axis("equal") # gives aspect ratio=1

    plt.figure(2)

    plt.pcolormesh(XX, YY, Z1 )
    plt.colorbar()
    plt.plot(x, y, "ok", label="input point")
    plt.legend()
    plt.title(" Precomputed Delaunay input " )
    plt.axis("equal")

    #plt.show()
    diff=Z1-Z2
    #print(diff[:,24])

    

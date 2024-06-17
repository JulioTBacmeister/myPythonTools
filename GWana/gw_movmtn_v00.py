workdir_ = "/glade/work/juliob/"
import sys

import numpy as np
import xarray as xr

"""
  !+++jtb
  boti=pver
  topi=desc%k-5

  do i=1,ncol
     udiff(i) = u(i,topi(i)) - usteer(i)
     vdiff(i) = v(i,topi(i)) - vsteer(i)
     do k=1,pver
        uwavef(i, k ) = u(i, k ) - usteer(i)
        vwavef(i, k ) = v(i, k ) - vsteer(i)
     end do
    ! xv(i) = -1._r8*udiff(i)/abs(udiff(i))
    ! yv(i) = -1._r8*vdiff(i)/abs(vdiff(i))
  end do
  call get_unit_vector(udiff , vdiff , xv, yv, ubisrc )
"""

gravit = 9.81

def get_unit_vector(u, v ):  #, u_n, v_n, mag)

    """
    real(r8), intent(in) :: u(:)
    real(r8), intent(in) :: v(:)
    real(r8), intent(out) :: u_n(:)
    real(r8), intent(out) :: v_n(:)
    real(r8), intent(out) :: mag(:)
    """
    u_n = np.zeros( len(u) )
    v_n = np.zeros( len(u) )
    mag = np.sqrt(u*u + v*v)

    # Has to be a loop/if instead of a where, because floating point
    # exceptions can trigger even on a masked divide-by-zero operation
    # (especially on Intel).
 
    for i in np.arange( len(mag) ):
        if (mag[i] > 0.):
            u_n[i] = u[i]/mag[i]
            v_n[i] = v[i]/mag[i]

        else:
            u_n[i] = 0.
            v_n[i] = 0.

    return u_n,v_n,mag

# Vectorized version of a 2D dot product (since the intrinsic dot_product
# is more suitable for arrays of contiguous vectors).
def dot_2d (u1, v1, u2, v2):

    dot_2d = u1*u2 + v1*v2

    return dot_2d

def midpoint_interp(arr): 
    ncol,nz = np.shape(arr)
    for i in np.arange(nz):
        interp[:,i] = 0.5 * ( arr[:,i]+arr[:,i+1] )
   
    return interp

def tends_from_tau(tau,pint,ubm):
    ncol,pverp = np.shape(tau)
    pver=pverp-1
    delp=np.zeros((ncol,pver))
    gwut=np.zeros((ncol,pver))
    c=np.zeros(ncol)
    for k in np.arange( pver ):
        delp[:,k] = pint[:,k+1] - pint[:,k] 

    for k in np.arange( pver ):
        ubtl = gravit * (tau[:,k+1]-tau[:,k]) / delp[:,k]
        # Assuming gwut, ubtl, c, and ubm are already defined NumPy arrays
        gwut[:, k ] = np.abs(ubtl) * np.sign(c[:] - ubm[:, k])

    return gwut

def movmtn_profiles(u,v,steer_level):

    ncol,pver = np.shape(u)
    # initialize vars to wirk with
    usteer = np.zeros(ncol)
    vsteer = np.zeros(ncol)
    topi = np.zeros(ncol).astype(int) 
    boti = np.zeros(ncol).astype(int)  
    # Initialize the Cell_Retro_Speed array
    Cell_Retro_Speed = np.zeros(ncol)

    
    isteer_level = np.round(steer_level).astype(int) 
    #Winds at 'steering level' 
    for i in np.arange( ncol ):
        usteer[i] = u[i, isteer_level[i] ]  
        vsteer[i] = v[i, isteer_level[i] ]

    xv_steer, yv_steer, umag_steer = get_unit_vector(usteer, vsteer ) 
    
    # Calculate the Cell_Retrograde_Speed for each column
    for i in range(ncol):
        Cell_Retro_Speed[i] = min(np.sqrt(usteer[i]**2 + vsteer[i]**2), 10.0)
    
    # Modify usteer and vsteer
    for i in range(ncol):
        usteer[i] = usteer[i] - xv_steer[i] * Cell_Retro_Speed[i]
        vsteer[i] = vsteer[i] - yv_steer[i] * Cell_Retro_Speed[i]

    boti[:]=pver-1
    topi=isteer_level-5

    
    # Assuming necessary arrays and variables are already defined
    # usteer, vsteer, xv_steer, yv_steer, u, v, topi, pver
    
    # Calculate CS1 and CS
    CS1 = np.sqrt(usteer**2 + vsteer**2)
    CS = CS1 * xv_steer + CS1 * yv_steer
    
    # Initialize arrays for udiff, vdiff, uwavef, vwavef
    udiff = np.zeros(ncol)
    vdiff = np.zeros(ncol)
    uwavef = np.zeros((ncol, pver))
    vwavef = np.zeros((ncol, pver))
    
    # Loop over ncol
    for i in range(ncol):
        udiff[i] = u[i, topi[i]] - usteer[i]
        vdiff[i] = v[i, topi[i]] - vsteer[i]
        for k in range(pver):
            uwavef[i, k] = u[i, k] - usteer[i]
            vwavef[i, k] = v[i, k] - vsteer[i]
    
    # Call get_unit_vector (assuming it's already defined in Python)
    # get_unit_vector(udiff, vdiff, xv, yv, ubisrc)
    xv, yv, ubisrc = get_unit_vector(udiff, vdiff )  #, xv, yv, ubisrc)


    # Assuming necessary arrays and variables are already defined
    # u, v, uwavef, vwavef, xv, yv, pver, ncol, topi
    
    # Initialize ubm array
    ubm = np.zeros((ncol, pver))
    
    # Looping through altitudes to calculate ubm
    for k in range(pver):
        ubm[:, k] = dot_2d(uwavef[:, k], vwavef[:, k], xv, yv)
    
    # Source level on-crest wind
    ubmsrc = np.zeros(ncol)
    for i in range(ncol):
        ubmsrc[i] = ubm[i, topi[i]]
    
    # Adjusting ubm, xv, yv based on the sign of ubmsrc
    for k in range(pver):
        for i in range(ncol):
            sign_factor = np.sign(ubmsrc[i])
            ubm[i, k] *= sign_factor
    
    # Adjust xv and yv
    xv *= np.sign(ubmsrc)
    yv *= np.sign(ubmsrc)
    
    # Set ubmsrc to its absolute value
    ubmsrc = np.abs(ubmsrc)

    return {
        'ubm': ubm,
        'uwavef': uwavef,
        'vwavef': vwavef,
        'vsteer': vsteer,
        'usteer': usteer,
        'vdiff' : vdiff,
        'udiff' : udiff,
        'xv': xv,
        'yv': yv,
        'xv_steer': xv_steer,
        'yv_steer': yv_steer
        }


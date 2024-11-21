############################################
# Python reproduction of gravity wave schemes in CAM.
# For simplicity we will stick with the Fortran-style dimensions,
# i.e., U(ncol,pver) etc ...
# ###########################################
import numpy as np
import sys
sys.path.append('../Utils/')
from shr_const import ShrConst as Cs

rair   = Cs.RDAIR
gravit = Cs.G

class BandType:
    def __init__(self, ngwv=None, dc=None, effkwv=None, kwv=None, fcrit2=1.0 ):
        self.ngwv = ngwv
        # Delta between nearest phase speeds [m/s].
        self.dc = dc
        if (effkwv==None):
            self.effkwv = kwv*fcrit2 #effkwv
        else:
            self.effkwv = effkwv
        self.kwv = kwv
        self.fcrit2 = fcrit2

class PType:
    def __init__(self, ifc=None, mid=None, rdel=None, delp=None, dst=None, rdst=None ):
        if ifc is not None :
            self.ifc = ifc
        else:
            print('You need to supply at least ifc')
            return
        if mid is None:
            self.mid = midpoint_interp( ifc )
        else:
            self.mid = mid
        if delp is None:
            self.delp = edgediff( ifc )
        else:
            self.delp = delp
        if dst is None:
            self.dst = edgediff( self.mid )
        else:
            self.dst = dst
        if rdel is None:
            self.rdel = 1./self.delp
        else:
            self.rdel = rdel
        if rdst is None:
            self.rdst = 1./self.dst
        else:
            self.rdst = rdst
        """
        self.delp = delp
        self.mid  = mid
        self.ifc  = ifc
        self.rdst = rdst
        """

#from gw_utils import midpoint_interp
def midpoint_interp(arr):
    """
    Interpolates the values of the input array along dimension 2.
    
    Parameters:
    arr (numpy.ndarray): Input 2D array of shape, e.g. (ncol, pver+1)
    
    Returns:
    numpy.ndarray: Interpolated array of shape (ncol, pver-1)
    """
    # Ensure input is a numpy array
    arr = np.asarray(arr)

    if ( len(arr.shape)==2 ):
        # Initialize the result array
        interp = np.zeros((arr.shape[0], arr.shape[1] - 1), dtype=arr.dtype)
        
        # Perform the interpolation
        interp = 0.5 * (arr[:, :-1] + arr[:, 1:])
    elif ( len(arr.shape)==2 ):
            # Initialize the result array
        interp = np.zeros((arr.shape[0] - 1), dtype=arr.dtype)
        
        # Perform the interpolation
        interp = 0.5 * (arr[:-1] + arr[1:])

    return interp
def edgediff(arr):
    """
    arr (numpy.ndarray): Input 2D array of shape, e.g. (ncol, pver+1)
    
    Returns:
    numpy.ndarray: diffs array of shape (ncol, pver-1)
    """
    # Ensure input is a numpy array
    arr = np.asarray(arr)

    if ( len(arr.shape)==2 ):
        # Initialize the result array
        diffs = np.zeros((arr.shape[0], arr.shape[1] - 1), dtype=arr.dtype)
        
        # Perform the interpolation
        diffs = (arr[:, :-1] - arr[:, 1:])
    elif ( len(arr.shape)==2 ):
            # Initialize the result array
        diffs = np.zeros((arr.shape[0] - 1), dtype=arr.dtype)
        
        # Perform the interpolation
        diffs = (arr[:-1] - arr[1:])

    return diffs
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


def gw_prof(ncol, p, cpair, t):
    # Initialize arrays for interface temperature, density, and Brunt-Vaisala frequencies
    pver=t.shape[1]
    ti = np.zeros((ncol, pver + 1))
    rhoi = np.zeros((ncol, pver + 1))
    ni = np.zeros((ncol, pver + 1))

    # Minimum value of Brunt-Vaisala frequency squared
    n2min = 5e-5

    # The top interface values are calculated assuming an isothermal atmosphere above the top level
    k = 0
    for i in range(ncol):
        ti[i, k] = t[i, k]
        rhoi[i, k] = p.ifc[i, k] / (rair * ti[i, k])
        ni[i, k] = np.sqrt(gravit * gravit / (cpair * ti[i, k]))

    # Interior points use centered differences
    ti[:, 1:pver] = midpoint_interp(t)
    for k in range(1,pver):
        for i in range(ncol):
            rhoi[i, k] = p.ifc[i, k] / (rair * ti[i, k])
            dtdp = (t[i, k] - t[i, k - 1]) * p.rdst[i, k - 1]
            n2 = gravit * gravit / ti[i, k] * (1 / cpair - rhoi[i, k] * dtdp)
            ni[i, k] = np.sqrt(max(n2min, n2))

    # Bottom interface uses bottom level temperature, density; next interface B-V frequency
    k = pver
    for i in range(ncol):
        ti[i, k] = t[i, k - 1]
        rhoi[i, k] = p.ifc[i, k] / (rair * ti[i, k])
        ni[i, k] = ni[i, k - 1]

    # Determine the midpoint Brunt-Vaisala frequencies
    nm = midpoint_interp(ni)

    return rhoi, nm, ni

def gw_drag_prof(ncol, band, p, src_level, tend_level, dt, 
                t, vramp,   
                piln, rhoi,    nm,   ni,  ubm,  ubi,  xv,    yv,   
                effgw,      c, kvtt, q,   dse,  tau,  utgw,  vtgw, 
                ttgw, qtgw, egwdffi,   gwut, dttdf, dttke, ro_adjust, 
                kwvrdg, satfac_in, lapply_effgw_in, lapply_vdiff, tau_diag ,
                perform_second_half=False, tau_0_ubc=False, do_vertical_diffusion=False ):

    pver =t.shape[1]

    """
    # Assuming all input variables are defined and have the correct dimensions
    # Initialize output and local arrays
    tau = np.zeros((ncol, 2 * ngwv + 1, pver + 1))
    utgw = np.zeros((ncol, pver))
    vtgw = np.zeros((ncol, pver))
    ttgw = np.zeros((ncol, pver))
    qtgw = np.zeros_like(q)  # Assuming q is already defined
    egwdffi = np.zeros((ncol, pver + 1))
    gwut = np.zeros((ncol, pver, 2 * ngwv + 1))
    dttdf = np.zeros((ncol, pver))
    dttke = np.zeros((ncol, pver))
    tau_diag = np.zeros((ncol, pver + 1))  # Optional
    """

    kbot_src = np.max( src_level )
    kbot_tend = np.max( tend_level )
    ktop = 0
    
    # Local storage
    d = np.zeros(ncol)
    mi = np.zeros(ncol)
    taudmp = np.zeros(ncol)
    tausat = np.zeros(ncol)
    ubmc = np.zeros(ncol)
    ubmc2 = np.zeros(ncol)
    ubt = np.zeros((ncol, pver))
    ubtl = np.zeros(ncol)
    wrk = np.zeros(ncol)
    ubt_lim_ratio = np.zeros(ncol)
    satfac = 2.0  # Default value

    dback=0.05
        
    # Boolean variables
    lapply_effgw = False
    do_vertical_diffusion = False    
    print( "Hey from gw_drag_prof " )
    # Loop from bottom to top to get stress profiles
    for k in range(kbot_src, ktop - 1,-1 ):
    
        # Determine the diffusivity for each column
        d = dback + kvtt[:, k]

        print( "K ",k)
        for l_ in range(-band.ngwv, band.ngwv + 1):
            print( "L" ,l_ )
            #return
            l= l_ + band.ngwv
            # Determine the absolute value of the saturation stress
            # Define critical levels where the sign of (u-c) changes between interfaces
            ubmc = ubi[:, k] - c[:, l]
            tausat = np.zeros_like(ubmc)

            # The lines below should simply ensure that if a critical exists between
            # levels k and k+1, tausat should be set= 0.
            if kwvrdg is not None:
                mask = src_level >= k
                sign_change_mask = np.logical_xor(ubmc > 0, ubi[:, k + 1] > c[:, l])
                tausat[mask & sign_change_mask] = np.abs(
                    kwvrdg * rhoi[mask, k] * ubmc[mask & sign_change_mask] ** 3 /
                    (satfac * ni[mask, k])
                )
            else:
                mask = src_level >= k
                sign_change_mask = np.logical_xor(ubmc > 0, ubi[:, k + 1] > c[:, l])
                tausat[mask & sign_change_mask] = np.abs(
                    band.effkwv * rhoi[mask, k] * ubmc[mask & sign_change_mask] ** 3 /
                    (satfac * ni[mask, k])
                )
    
            if ro_adjust is not None:
                mask = src_level >= k
                tausat[mask] *= np.sqrt(ro_adjust[:, l, k])
    
            if kwvrdg is not None:
                mask = src_level >= k
                ubmc2 = np.maximum(ubmc ** 2, ubmc2mn)
                mi = ni[:, k] / (2 * kwvrdg * ubmc2) * (alpha[k] + ni[:, k] ** 2 / ubmc2 * d)
                wrk = -2 * mi * rog * t[:, k] * (piln[:, k + 1] - piln[:, k])
                taudmp = tau[:, l, k + 1]
    
                tausat[tausat <= taumin] = 0
                taudmp[taudmp <= taumin] = 0
    
                tau[:, l, k] = np.minimum(taudmp, tausat)
    
            else:
                mask = src_level >= k
                ubmc2 = np.maximum(ubmc ** 2, ubmc2mn)
                mi = ni[:, k] / (2 * band.kwv * ubmc2) * (alpha[k] + ni[:, k] ** 2 / ubmc2 * d)
                wrk = -2 * mi * rog * t[:, k] * (piln[:, k + 1] - piln[:, k])
                taudmp = tau[:, l, k + 1] * np.exp(wrk)
    
                tausat[tausat <= taumin] = 0
                taudmp[taudmp <= taumin] = 0
    
                tau[:, l, k] = np.minimum(taudmp, tausat)





    if not perform_second_half:
        return

    # Force tau at the top of the model to zero, if requested
    if tau_0_ubc:
        tau[:, :, ktop] = 0.0

    # Write out pre-adjustment tau profile for diagnostic purposes
    if tau_diag is not None:
        tau_diag[:, :] = tau[:, 0, :]

    # Apply efficiency to completed stress profile
    if lapply_effgw_in:
        for k in range(ktop, kbot_tend + 1):
            for l_ in range(-band.ngwv, band.ngwv + 1):
                l= l_ + band.ngwv
                mask = (k - 1 <= tend_level)
                tau[mask, l, k] = tau[mask, l, k] * effgw

    # Compute the tendencies from the stress divergence
    for k in range(ktop, kbot_tend):
        ubt = np.zeros_like(ubm[:, k])
        for l_ in range(-band.ngwv, band.ngwv + 1):
            l= l_ + band.ngwv
            ubtl = gravit * (tau[:, l, k + 1] - tau[:, l, k]) * p.rdel[:, k]
            ubtl = np.minimum(ubtl, umcfac * np.abs(c[:, l] - ubm[:, k]) / dt)
            if not lapply_effgw_in:
                ubtl = np.minimum(ubtl, tndmax)
            mask = (k <= tend_level)
            gwut[mask, k, l] = np.sign(ubtl, c[:, l] - ubm[:, k])
            ubt[mask, k] = ubt[mask, k] + gwut[mask, k, l]

        if lapply_effgw_in:
            mask = np.abs(ubt[:, k]) > tndmax
            ubt_lim_ratio = np.where(mask, tndmax / np.abs(ubt[:, k]), 1.0)
            ubt[:, k] = ubt_lim_ratio * ubt[:, k]
        else:
            ubt_lim_ratio = 1.0

        for l_ in range(-band.ngwv, band.ngwv + 1):
            l= l_ + band.ngwv
            gwut[:, k, l] = ubt_lim_ratio * gwut[:, k, l]
            mask = np.abs(gwut[:, k, l]) < 1.e-15
            gwut[mask, k, l] = 0.0
            mask = (k <= tend_level)
            tau[mask, l, k + 1] = tau[mask, l, k] + np.abs(gwut[mask, k, l]) * p.del_field[mask, k] / gravit

        mask = (k <= tend_level)
        utgw[mask, k] = ubt[mask, k] * xv
        vtgw[mask, k] = ubt[mask, k] * yv

        if vramp is not None:
            utgw[:, k] *= vramp[k]
            vtgw[:, k] *= vramp[k]

    if not lapply_effgw_in:
        for k in range(ktop, kbot_tend + 1):
            for l_ in range(-band.ngwv, band.ngwv + 1):
                l= l_ + band.ngwv
                mask = (k - 1 <= tend_level)
                tau[mask, l, k] = tau[mask, l, k] * effgw
        for k in range(ktop, kbot_tend):
            for l_ in range(-band.ngwv, band.ngwv + 1):
                l= l_ + band.ngwv
                gwut[:, k, l] *= effgw

            utgw[:, k] *= effgw
            vtgw[:, k] *= effgw

    if do_vertical_diffusion:
        gw_ediff(ncol, pver, band.ngwv, kbot_tend, ktop, tend_level, gwut, ubm, nm, rhoi, dt, prndl, gravit, p, c, vramp, egwdffi, decomp, ro_adjust=ro_adjust)
        for m in range(q.shape[2]):
            gw_diff_tend(ncol, pver, kbot_tend, ktop, q[:, :, m], dt, decomp, qtgw[:, :, m])
        gw_diff_tend(ncol, pver, kbot_tend, ktop, dse, dt, decomp, dttdf)

    for l_ in range(-band.ngwv, band.ngwv + 1):
        for k in range(ktop, kbot_tend):
            dttke[:, k] = dttke[:, k] - (ubm[:, k] - c[:, l]) * gwut[:, k, l]

    ttgw = dttke + dttdf

    if vramp is not None:
        for k in range(ktop, kbot_tend):
            ttgw[:, k] *= vramp[k]


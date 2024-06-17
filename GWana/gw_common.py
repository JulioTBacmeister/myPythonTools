############################################
# Python reproduction of gravity wave schemes in CAM.
# For simplicity we will stick with the Fortran-style dimensions,
# i.e., U(ncol,pver) etc ...
############################################
import numpy as np
import sys
sys.path.append('../Utils/')
from shr_const import ShrConst as Cs

rair   = Cs.RDAIR
gravit = Cs.G
rog    = rair/gravit

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

def rad_alphas(pref):
    from scipy.interpolate import interp1d

    # Define the number of alpha values
    nalph = 71
    
    # Define alpha0 array
    alpha0 = np.array([
        0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.10133333, 0.104,
        0.108, 0.112, 0.116, 0.12066667,
        0.126, 0.132, 0.138, 0.144,
        0.15133333, 0.16, 0.17, 0.18,
        0.19, 0.19933333, 0.208, 0.216,
        0.224, 0.232, 0.23466667, 0.232,
        0.224, 0.216, 0.208, 0.20133333,
        0.196, 0.192, 0.188, 0.184,
        0.18266667, 0.184, 0.188, 0.192,
        0.196, 0.19333333, 0.184, 0.168,
        0.152, 0.136, 0.12133333, 0.108,
        0.096, 0.084, 0.072, 0.061,
        0.051, 0.042, 0.033, 0.024,
        0.017666667, 0.014, 0.013, 0.012,
        0.011, 0.010333333, 0.01, 0.01,
        0.01, 0.01, 0.01
    ], dtype=np.float64)
    
    # Define palph array (pressure levels used to calculate alpha0 in hPa)
    palph = np.array([
        2.06115E-06, 2.74280E-06, 3.64988E-06, 4.85694E-06,
        6.46319E-06, 8.60065E-06, 1.14450E-05, 1.52300E-05,
        2.02667E-05, 2.69692E-05, 3.58882E-05, 4.77568E-05,
        6.35507E-05, 8.45676E-05, 0.000112535, 0.000149752,
        0.000199277, 0.000265180, 0.000352878, 0.000469579,
        0.000624875, 0.000831529, 0.00110653, 0.00147247,
        0.00195943, 0.00260744, 0.00346975, 0.00461724,
        0.00614421, 0.00817618, 0.0108801, 0.0144783,
        0.0192665, 0.0256382, 0.0341170, 0.0453999,
        0.0604142, 0.0803939, 0.106981, 0.142361,
        0.189442, 0.252093, 0.335463, 0.446404,
        0.594036, 0.790490, 1.05192, 1.39980,
        1.86273, 2.47875, 3.29851, 4.38936,
        5.84098, 7.77266, 10.3432, 13.7638,
        18.3156, 24.3728, 32.4332, 43.1593,
        57.4326, 76.4263, 101.701, 135.335,
        180.092, 239.651, 318.907, 424.373,
        564.718, 751.477, 1000.0
    ], dtype=np.float64)

    # Assuming alpha0 and palph are already defined as in the previous code

    # Constants
    seconds_in_day = 86400.0
    min_alpha0_value = 1.e-6
    conversion_factor = 1.e2
    
    # Loop through each element and apply the operations
    for k in range(nalph):
        alpha0[k] = alpha0[k] / seconds_in_day
        alpha0[k] = max(alpha0[k], min_alpha0_value)
        #palph[k] = palph[k] * conversion_factor
    

    # Assuming alpha0 and palph are already defined 
    
    # Create an interpolation function based on palph and alpha0
    interp_func = interp1d(palph, alpha0, kind='linear', bounds_error=False, fill_value='extrapolate')
    
    # Interpolate alpha0 onto the new pressure grid pref
    alpha0_interp = interp_func(pref)
    
    # Print the interpolated values (optional)
    # print("New pressure grid (pref):", pref)
    # print("Interpolated alpha0 values:", alpha0_interp)

    return alpha0_interp
    

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
        
        # Perform the difference in top-down CAM vertical style
        diffs = (arr[:,1:] - arr[:,:-1])
    elif ( len(arr.shape)==2 ):
            # Initialize the result array
        diffs = np.zeros((arr.shape[0] - 1), dtype=arr.dtype)
        
        # Perform the difference in top-down CAM vertical style
        diffs = (arr[1:] - arr[:-1])

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

    """
    
  ! Interior points use centered differences.
  ti(:,2:pver) = midpoint_interp(t)
  do k = 2, pver
     do i = 1, ncol
        rhoi(i,k) = p%ifc(i,k) / (rair*ti(i,k))
        dtdp = (t(i,k)-t(i,k-1)) * p%rdst(i,k-1)
        n2 = gravit*gravit/ti(i,k) * (1._r8/cpair - rhoi(i,k)*dtdp)
        ni(i,k) = sqrt(max(n2min, n2))
     end do
  end do
    """
    # Interior points use centered differences
    ti[:, 1:pver] = midpoint_interp(t)
    for k in range(1,pver):
        for i in range(ncol):
            rhoi[i, k] = p.ifc[i, k] / (rair * ti[i, k])
            dtdp = (t[i, k] - t[i, k - 1]) * p.rdst[i, k - 1]
            n2 = gravit * gravit / ti[i, k] * (1 / cpair - rhoi[i, k] * dtdp)
            ni[i, k] =  np.sqrt(max(n2min, n2))

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
    ubmc2mn = 0.01
    taumin = 1.e-10

    # For now
    #alpha=np.zeros( pver )
    pref_edge=1_000.*p.ifc[0,:]/p.ifc[0,93]
    alpha = rad_alphas(pref_edge)
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
                for i in np.arange( ncol ):
                    if ( (src_level[i]>=k) and ( ubmc[i] *(ubi[i, k + 1] - c[i, l])>0.) ):
                        tausat[i] = np.abs(
                            kwvrdg[i] * rhoi[i, k] * ubmc[i] ** 3 /
                            (satfac * ni[i, k])
                        )
            else:
                for i in np.arange( ncol ):
                    if ( (src_level[i]>=k) and ( ubmc[i] *(ubi[i, k + 1] - c[i, l])>0.) ):
                        tausat[i] = np.abs(
                            band.effkwv * rhoi[i, k] * ubmc[i] ** 3 /
                            (satfac * ni[i, k])
                        )
    
            if ro_adjust is not None:
                mask = src_level >= k
                tausat[mask] *= np.sqrt(ro_adjust[:, l, k])
    
            if kwvrdg is not None:
                for i in np.arange( ncol ):
                    if( src_level[i]>=k):
                        ubmc2[i] = np.maximum(ubmc[i] ** 2, ubmc2mn)
                        mi[i] = ni[i, k] / (2 * kwvrdg[i] * ubmc2[i]) * (alpha[k] + ni[i, k] ** 2 / ubmc2[i] * d[i])
                        wrk = -2 * mi * rog * t[i, k] * (piln[i, k + 1] - piln[i, k])
                        taudmp = tau[i, l, k + 1]
            
                        if (tausat[i]<=taumin):
                            tausat[i]=0.
                        if (taudmp[i]<=taumin):
                            taudmp[i]=0.
            
                        tau[i, l, k] = np.minimum(taudmp[i], tausat[i])
    
            else:
                for i in np.arange( ncol ):
                    if( src_level[i]>=k):
                        ubmc2[i] = np.maximum(ubmc[i] ** 2, ubmc2mn)
                        mi[i] = ni[i, k] / (2 * band.kwv * ubmc2[i]) * (alpha[k] + ni[i, k] ** 2 / ubmc2[i] * d[i])
                        wrk = -2 * mi * rog * t[i, k] * (piln[i, k + 1] - piln[i, k])
                        taudmp = tau[i, l, k + 1] * np.exp( wrk )
            
                        if (tausat[i]<=taumin):
                            tausat[i]=0.
                        if (taudmp[i]<=taumin):
                            taudmp[i]=0.
            
                        tau[i, l, k] = np.minimum(taudmp[i], tausat[i])

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


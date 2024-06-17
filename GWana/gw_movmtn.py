import numpy as np

class MovMtnSourceDesc:
    def __init__(self, storm_shift, k, min_hdepth, maxh, maxuh, hd=None, uh=None, mfcc=None):
        # Whether wind speeds are shifted to be relative to storm cells.
        self.storm_shift = storm_shift
        # Index for level where wind speed is used as the source speed. ->700hPa
        self.k = k
        # Heating depths below this value [m] will be ignored.
        self.min_hdepth = min_hdepth
        # Table bounds, for convenience. (Could be inferred from shape(mfcc).)
        self.maxh = maxh  # -bounds of the lookup table heating depths 
        self.maxuh = maxuh  # bounds of the lookup table wind
        # Heating depths [m].
        self.hd = hd if hd is not None else []
        self.uh = uh if uh is not None else []
        # Table of source spectra.
        self.mfcc = mfcc if mfcc is not None else []

# Example of how to initialize the class
# mov_mtn = MovMtnSourceDesc(storm_shift=True, k=5, min_hdepth=100.0, maxh=10, maxuh=10, hd=[0.0, 1.0, 2.0], uh=[0.1, 0.2, 0.3], mfcc=[[[0.0]*5]*5]*5)

#def gw_movmtn_src(ncol, lchnk , band, desc, u, v, netdt, netdt_shcu, xpwp_shcu, zm, src_level, tend_level, tau, ubm, ubi, xv, yv, c, hdepth):
#                          *             *
def gw_movmtn_src(ncol , band , u, v, netdt, netdt_shcu, xpwp_shcu, zm, src_level, tend_level, tau, ubm, ubi, xv, yv, c, hdepth):
    #-----------------------------------------------------------------------
    # Driver for multiple gravity wave drag parameterization.
    #-----------------------------------------------------------------------

    # Import necessary utility functions
    from gw_common import get_unit_vector, dot_2d, midpoint_interp

    #---------------------------Local Storage-------------------------------
    # Initialize arrays and local variables
    pver = u.shape[1]
    usteer = np.zeros(ncol)
    vsteer = np.zeros(ncol)
    uwavef = np.zeros((ncol, pver))
    vwavef = np.zeros((ncol, pver))
    steer_level = np.zeros(ncol)
    Cell_Retro_Speed = np.zeros(ncol)
    q0 = np.zeros(ncol)
    qj = np.zeros(ncol)
    xv_steer = np.zeros(ncol)
    yv_steer = np.zeros(ncol)
    umag_steer = np.zeros(ncol)
    boti = np.zeros(ncol, dtype=int)
    topi = np.zeros(ncol, dtype=int)
    hd_idx = np.zeros(ncol, dtype=int)
    uh = np.zeros(ncol)
    Umini = np.zeros(ncol, dtype=int)
    Umaxi = np.zeros(ncol, dtype=int)
    tau0 = np.zeros(2*band.ngwv+1)
    CS = np.zeros(ncol)
    CS1 = np.zeros(ncol)
    udiff = np.zeros(ncol)
    vdiff = np.zeros(ncol)
    ubmsrc = np.zeros(ncol)
    shift = 0
    ut = np.zeros(ncol)
    uc = np.zeros(ncol)
    umm = np.zeros(ncol)
    taumm = np.zeros(ncol)
    CF = 20.0  # Conversion factor
    AL = 1.0e5  # Averaging length
    hdmm_idx = np.zeros(ncol)
    uhmm_idx = np.zeros(ncol)
    c_idx = np.zeros((ncol, 2*band.ngwv+1))
    xpwp_src = np.zeros(ncol)
    Steer_k = pver - 2  #- 1

    # Initialize tau array and other arrays
    tau.fill(0.0)
    hdepth.fill(0.0)
    q0.fill(0.0)
    tau0.fill(0.0)

    # Calculate flux source from ShCu/PBL
    xpwp_src = xpwp_shcu # shcu_flux_src(xpwp_shcu, ncol, pver + 1)

    # Determine wind and unit vectors at the source (steering level)
    usteer = u[:, Steer_k]
    vsteer = v[:, Steer_k]
    steer_level.fill(Steer_k)
    xv_steer, yv_steer, umag_steer = get_unit_vector(usteer, vsteer)

    # Account for retrograde cell motion
    for i in range(ncol):
        Cell_Retro_Speed[i] = min(np.sqrt(usteer[i]**2 + vsteer[i]**2), 0.0)

    for i in range(ncol):
        usteer[i] -= xv_steer[i] * Cell_Retro_Speed[i]
        vsteer[i] -= yv_steer[i] * Cell_Retro_Speed[i]

    # Calculate heating depth
    boti.fill(pver)
    topi.fill(Steer_k - 10) #- 10)
    
    """
    hdepth = zm[:, topi] - zm[:, boti]
    hd_idx = np.searchsorted(desc.hd, hdepth, side='left')

    # Ensure heating depth is large enough
    hd_idx[hdepth < np.maximum(desc.min_hdepth, desc.hd[0])] = 0

    # Find the maximum heating rate
    for k in range(topi.min(), boti.max() + 1):
        q0 = np.maximum(q0, netdt[:, k])

    q0 *= CF
    qj = 9.81 / 285 * q0
    """
    
    CS1 = np.sqrt(usteer**2 + vsteer**2)
    CS = CS1 * xv_steer + CS1 * yv_steer

    # Calculate winds in the reference frame of the wave
    for i in range(ncol):
        udiff[i] = u[i, topi[i]] - usteer[i]
        vdiff[i] = v[i, topi[i]] - vsteer[i]
        uwavef[i, :] = u[i, :] - usteer[i]
        vwavef[i, :] = v[i, :] - vsteer[i]

    # Wave relative wind at source level
    for i in range(ncol):
        udiff[i] = uwavef[i, topi[i]]
        vdiff[i] = vwavef[i, topi[i]]

    # Unit vector components in the direction of the wavevector
    xv, yv, ubisrc = get_unit_vector(udiff, vdiff)

    """
    # --- CAM ouptut routines not available in Python ---
    # Output fields for diagnostics
    outfld('UCELL_MOVMTN', usteer, ncol, lchnk)
    outfld('VCELL_MOVMTN', vsteer, ncol, lchnk)
    outfld('CS_MOVMTN', CS, ncol, lchnk)
    outfld('CS1_MOVMTN', CS1, ncol, lchnk)
    outfld('STEER_LEVEL_MOVMTN', steer_level, ncol, lchnk)
    outfld('XPWP_SRC_MOVMTN', xpwp_src, ncol, lchnk)
    """

    # Project the local wave relative wind at midpoints onto the direction of the wavevector
    for k in range(pver):
        ubm[:, k] = dot_2d(uwavef[:, k], vwavef[:, k], xv, yv)

    # Source level on-crest wind
    for i in range(ncol):
        ubmsrc[i] = ubm[i, topi[i]]

    # Adjust the wave relative wind and unit vector components to be positive
    for k in range(pver):
        ubm[:, k] *= np.sign(ubmsrc)

    xv *= np.sign(ubmsrc)
    yv *= np.sign(ubmsrc)

    # Compute the interface wind projection by averaging the midpoint winds
    ubi[:, 1] = ubm[:, 0]
    ubi[:, 1:pver] = midpoint_interp(ubm)

    # Determine wind for the lookup table
    for i in range(ncol):
        ut[i] = ubm[i, topi[i]]
        uh[i] = ut[i] - CS[i]

    # Set phase speeds
    c[:, 0] = 0.0

    # Gravity wave sources
    for i in range(ncol):
        tau[i, 0, topi[i]:pver + 1] = xpwp_src[i]

    # Output the source level
    src_level[:] = topi
    tend_level[:] = topi


def index_of_nearest(x, grid):
    """
    Find the index of the nearest grid point for each value in x.

    Parameters:
    x (numpy.ndarray): 1D array of values.
    grid (numpy.ndarray): 1D array of grid points.

    Returns:
    numpy.ndarray: Array of indices of the nearest grid points.
    """
    n = len(grid)
    interfaces = (grid[:-1] + grid[1:]) / 2.0

    idx = np.ones(len(x), dtype=int)
    for i in range(n - 1):
        idx[np.where(x > interfaces[i])] = i + 1

    return idx

def shcu_flux_src(xpwp_shcu, ncol, pverx):
    """
    Calculate the flux source from shallow convection (ShCu).

    Parameters:
    xpwp_shcu (numpy.ndarray): 2D array of ShCu flux (ncol, pverx).
    ncol (int): Number of columns.
    pverx (int): Vertical extent.

    Returns:
    numpy.ndarray: 1D array of flux source from ShCu.
    """
    alpha_shcu_flux = 0.01  # Tunable parameter

    # Simple average over layers
    nlayers = 5
    xpwp_src = np.zeros(ncol)
    for k in range(nlayers):
        xpwp_src += xpwp_shcu[:, pverx - k - 1]

    xpwp_src = alpha_shcu_flux * xpwp_src / nlayers

    return xpwp_src



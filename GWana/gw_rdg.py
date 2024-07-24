import numpy as np

def gw_rdg_resid_src(ncol, band, p, u, v, t, mxdis, kwvrdg, zi, nm):
    # Initialize output arrays
    src_level = np.zeros(ncol, dtype=int)
    tend_level = np.zeros(ncol, dtype=int)
    tau = np.zeros((ncol, -band.ngwv, band.ngwv, p.shape[1] + 1))
    ubm = np.zeros((ncol, p.shape[1]))
    ubi = np.zeros((ncol, p.shape[1] + 1))
    xv = np.zeros(ncol)
    yv = np.zeros(ncol)
    ubmsrc = np.zeros(ncol)
    usrc = np.zeros(ncol)
    vsrc = np.zeros(ncol)
    nsrc = np.zeros(ncol)
    rsrc = np.zeros(ncol)
    m2src = np.zeros(ncol)
    c = np.zeros((ncol, -band.ngwv, band.ngwv))

    # Constants
    rair = 287.05  # specific gas constant for dry air
    Fcrit_res = 1.0
    orom2min = 0.01

    # Surface streamline displacement height (2*sgh).
    hdsp = mxdis  # no longer multiplied by 2
    src_level[:] = p.shape[1] + 1  # assuming pver is p.shape[1]

    tau[:, 0, :] = 0.0

    # Find depth of "source layer" for mountain waves
    # i.e., between ground and mountain top
    for k in range(p.shape[1] - 1, -1, -1):
        for i in range(ncol):
            if hdsp[i] >= zi[i, k + 1] and hdsp[i] < zi[i, k]:
                src_level[i] = k

    rsrc[:] = 0.0
    usrc[:] = 0.0
    vsrc[:] = 0.0
    nsrc[:] = 0.0

    for i in range(ncol):
        for k in range(p.shape[1] - 1, src_level[i] - 1, -1):
            rsrc[i] += p['mid'][i, k] / (rair * t[i, k]) * p['del'][i, k]
            usrc[i] += u[i, k] * p['del'][i, k]
            vsrc[i] += v[i, k] * p['del'][i, k]
            nsrc[i] += nm[i, k] * p['del'][i, k]

    dpsrc = np.zeros(ncol)
    for i in range(ncol):
        dpsrc[i] = p['ifc'][i, p.shape[1] + 1] - p['ifc'][i, src_level[i]]

    rsrc /= dpsrc
    usrc /= dpsrc
    vsrc /= dpsrc
    nsrc /= dpsrc

    # Get the unit vector components and magnitude at the surface.
    xv, yv, wmsrc = get_unit_vector(usrc, vsrc)
    ubmsrc = wmsrc

    # Project the local wind at midpoints onto the source wind.
    for k in range(p.shape[1]):
        ubm[:, k] = dot_2d(u[:, k], v[:, k], xv, yv)

    # Compute the interface wind projection by averaging the midpoint winds.
    ubi[:, 1] = ubm[:, 1]
    ubi[:, 2:p.shape[1]] = midpoint_interp(ubm)
    ubi[:, p.shape[1]] = ubm[:, p.shape[1] - 1]

    m2src = ((nsrc / (ubmsrc + 0.01)) ** 2 - kwvrdg ** 2) / ((nsrc / (ubmsrc + 0.01)) ** 2)

    ubi[:, 1] = ubm[:, 1]
    ubi[:, 2:p.shape[1]] = midpoint_interp(ubm)
    ubi[:, p.shape[1] + 1] = ubm[:, p.shape[1] - 1]

    tauoro = np.zeros(ncol)

    for i in range(ncol):
        if src_level[i] > 0 and m2src[i] > orom2min:
            sghmax = Fcrit_res * (ubmsrc[i] / nsrc[i]) ** 2
            tauoro[i] = 0.5 * kwvrdg[i] * min(hdsp[i] ** 2, sghmax) * rsrc[i] * nsrc[i] * ubmsrc[i]
        else:
            tauoro[i] = 0.0

    tend_level[:] = p.shape[1]
    c[:] = 0.0

    return src_level, tend_level, tau, ubm, ubi, xv, yv, ubmsrc, usrc, vsrc, nsrc, rsrc, m2src, c

# Utility functions
def dot_2d(u, v, xv, yv):
    return u * xv + v * yv

def midpoint_interp(ubm):
    return (ubm[:, :-1] + ubm[:, 1:]) / 2

def get_unit_vector(usrc, vsrc):
    magnitude = np.sqrt(usrc**2 + vsrc**2)
    xv = usrc / magnitude
    yv = vsrc / magnitude
    return xv, yv, magnitude

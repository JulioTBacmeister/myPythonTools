import numpy as np
from gw_utils import midpoint_interp

def gw_prof(ncol, p, cpair, t):
    # Initialize arrays for interface temperature, density, and Brunt-Vaisala frequencies
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
    for k in range(1, p.ver):
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
     kwvrdg, satfac_in, lapply_effgw_in, lapply_vdiff, tau_diag ):

    # Assuming all variables are defined and have the correct dimensions
    
    # Loop from bottom to top to get stress profiles
    for k in range(kbot_src, ktop + 1)[::-1]:
    
        # Determine the diffusivity for each column
        d = dback + kvtt[:, k]
    
        for l in range(-band.ngwv, band.ngwv + 1):
    
            # Determine the absolute value of the saturation stress
            # Define critical levels where the sign of (u-c) changes between interfaces
            ubmc = ubi[:, k] - c[:, l]
            tausat = np.zeros_like(ubmc)
    
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
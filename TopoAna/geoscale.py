import numpy as np
R_EARTH = 6371.22

def ridge_scale_from_fields(fields, sigmas_km, dx_km=3.336, gamma=0.75):
    """Index-space version (as now): derivatives in cube cells, 3.336 km assumed."""
    return _scale(fields, sigmas_km, gamma, J=None, dx_km=dx_km)

def local_jacobian(lon, lat):
    """km per cell of local east/north coordinates, per pixel: (dE/dx, dE/dy, dN/dx, dN/dy);
    x = array axis 1, y = array axis 0. lon/lat in degrees, any longitude convention."""
    lonc = np.nanmean(lon)
    lw = (lon - lonc + 180.) % 360. - 180.
    dlon_dy, dlon_dx = np.gradient(lw)
    dlat_dy, dlat_dx = np.gradient(lat)
    kE = R_EARTH * np.radians(1.) * np.cos(np.radians(lat))
    kN = R_EARTH * np.radians(1.)
    return kE*dlon_dx, kE*dlon_dy, kN*dlat_dx, kN*dlat_dy

def ridge_scale_from_fields_geo(fields, sigmas_km, lon, lat, gamma=0.75):
    """Same analysis with the curvatures converted to true local east/north distances.
    th_cross is then geographic: degrees from east toward north."""
    return _scale(fields, sigmas_km, gamma, J=local_jacobian(lon, lat))

def _scale(fields, sigmas_km, gamma, J=None, dx_km=3.336, aniso=0.5):
    shape = fields[0].shape
    R_best = np.full(shape, -np.inf); R1_best = np.full(shape, np.nan)
    sig_best = np.full(shape, np.nan); th_cross = np.full(shape, np.nan)
    if J is not None:
        a, b, c, d = J; det = a*d - b*c
        ia, ib, ic, idd = d/det, -b/det, -c/det, a/det        # inverse Jacobian: d(x,y)/d(E,N)
    for F, s_km in zip(fields, sigmas_km):
        Fy, Fx = np.gradient(F)
        Lyy = np.gradient(Fy, axis=0); Lxx = np.gradient(Fx, axis=1); Lxy = np.gradient(Fx, axis=0)
        if J is None:
            hxx, hyy, hxy, s = Lxx, Lyy, Lxy, s_km/dx_km             # per cell^2, sigma in cells
        else:                                                          # per km^2, sigma in km
            hxx = ia*ia*Lxx + 2*ia*ic*Lxy + ic*ic*Lyy                  # d2/dE2
            hyy = ib*ib*Lxx + 2*ib*idd*Lxy + idd*idd*Lyy               # d2/dN2
            hxy = ia*ib*Lxx + (ia*idd + ic*ib)*Lxy + ic*idd*Lyy        # d2/dEdN
            s = s_km
        half_tr = 0.5*(hxx + hyy); rad = np.sqrt(0.25*(hxx - hyy)**2 + hxy**2)
        lam_q, lam_p = half_tr - rad, half_tr + rad
        ridge_like = (lam_q < 0) & (np.abs(lam_p) < aniso*np.abs(lam_q))
        Rg = np.where(ridge_like, (s**2)**gamma * (-lam_q), -np.inf)
        R1 = (s**2) * (-lam_q)                                         # metres (~0.35 x height)
        better = Rg > R_best
        R_best = np.where(better, Rg, R_best); R1_best = np.where(better, R1, R1_best)
        sig_best = np.where(better, s_km, sig_best)
        th_p = 0.5*np.degrees(np.arctan2(2*hxy, hxx - hyy))
        th_cross = np.where(better, (th_p + 90.) % 180., th_cross)
    return sig_best, R1_best, th_cross
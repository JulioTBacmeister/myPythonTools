#!/usr/bin/env python3
"""
paint_tiles_nc.py
=================

Paint a per-tile quantity A[i, ir] back onto the high-resolution cubed sphere,
taking the ridge geometry from the **CAM-style tile topography netCDF**, e.g.

    ne30pg3_gmted2010_..._tiles.nc

rather than from Ridge_tile_map.dat. i is the SE column index (ncol), ir the
ridge index (nrdg), matching the arrays the gravity-wave driver writes.

Each tile is rendered as a convex polygon centred on its crest centroid,
oriented along the crest, with along-crest extent from CLNGT and cross-crest
extent from HWDTH. Where tiles overlap, the cell keeps the LARGER value of A.

Getting cube coordinates
------------------------
The tiles netCDF stores the crest centroid as latc/lonc in RADIANS, not as
cube indices. remapridge2tiles produced those by

    alph = ( xloid - 0.5 - ncube/2 ) * (pi/2) / ncube
    beta = ( yloid - 0.5 - ncube/2 ) * (pi/2) / ncube
    call CubedSphereRLLFromABP( alph, beta, panel, lonc, latc )

so this module inverts it: lonc/latc -> (alpha, beta, panel) by the same
gnomonic assignment CubedSphereABPFromRLL uses, then

    xloid = alpha * ncube/(pi/2) + 0.5 + ncube//2

Note ncube//2 is integer division, matching Fortran.

If the file already carries xloid/yloid/panel (cheaper and exact), those are
used instead and the inversion is skipped. Add them on the Fortran side with:

    call defvar('XLOID',2,xloidid)      ! and YLOID, PANEL
    call put2d(xloidid,'XLOID',clean(xloid_tiles))

after promoting xloid_tiles/yloid_tiles/panel_tiles to module scope the way
latc_tiles/lonc_tiles were.

Validity
--------
The writer's clean() zeroes unused tile slots, so there is no NumObjects on the
file. A slot counts as a real ridge when MXDIS > 0; paintridge2cube only paints
ridges with mxdis >= 1, so this is safe.

Usage
-----
::

    from paint_tiles_nc import load_tiles_nc, paint_tiles_nc

    tg = load_tiles_nc('..._tiles.nc', ncube=3000, grid_length_scale=3.3)
    A  = taudsw - tauoro                      # (ncol, nrdg)
    field = paint_tiles_nc(A, tg, shape='diamond', verbose=True)
    # field is (6, ncube, ncube)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["TileGeom", "load_tiles_nc", "paint_tiles_nc",
           "abp_from_rll", "rll_from_abp", "cube_index_from_lonlat",
           "lonlat_from_cube_index", "downwind_sign", "tile_polygon",
           "MODULATIONS"]

SHAPES = ("diamond", "rectangle", "triangle", "wedge_end")


# ---------------------------------------------------------------------------
# Cubed-sphere inverse transform, mirroring CubedSphereABPFromRLL
# ---------------------------------------------------------------------------
def abp_from_rll(lon, lat, rotate_cube=0.0):
    """(lon, lat) in radians -> (alpha, beta, panel), panel 1-6.

    Vectorised transcription of CubedSphereABPFromRLL in cube_to_target.F90.
    Panel numbering: 1 (+x), 2 (+y), 3 (-x), 4 (-y), 5 (-z), 6 (+z).
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)

    xx = np.cos(lon - rotate_cube) * np.cos(lat)
    yy = np.sin(lon - rotate_cube) * np.cos(lat)
    zz = np.sin(lat)

    pm = np.maximum.reduce([np.abs(xx), np.abs(yy), np.abs(zz)])
    ix = np.where(pm == np.abs(xx), np.sign(xx), 0.0).astype(int)
    iy = np.where(pm == np.abs(yy), np.sign(yy), 0.0).astype(int)
    iz = np.where(pm == np.abs(zz), np.sign(zz), 0.0).astype(int)

    panel = np.zeros(lon.shape, dtype=int)
    sx = np.zeros(lon.shape); sy = np.zeros(lon.shape); sz = np.zeros(lon.shape)

    def setp(mask, p, a, b, c):
        panel[mask] = p
        sx[mask] = a[mask]; sy[mask] = b[mask]; sz[mask] = c[mask]

    # Order matters: exactly as in the Fortran if/elseif chain.
    m = (iz == 1)
    setp(m, 6, yy, -xx, zz)
    m = (panel == 0) & (iz == -1)
    setp(m, 5, yy, xx, -zz)
    m = (panel == 0) & (ix == 1) & (iy != 1)
    setp(m, 1, yy, zz, xx)
    m = (panel == 0) & (ix == -1) & (iy != -1)
    setp(m, 3, -yy, zz, -xx)
    m = (panel == 0) & (iy == 1) & (ix != -1)
    setp(m, 2, -xx, zz, yy)
    m = (panel == 0) & (iy == -1) & (ix != 1)
    setp(m, 4, xx, zz, -yy)

    if np.any(panel == 0):
        raise ValueError(f"{int((panel==0).sum())} point(s) failed panel assignment")

    return np.arctan(sx / sz), np.arctan(sy / sz), panel


def rll_from_abp(alpha, beta, panel, rotate_cube=0.0):
    """(alpha, beta, panel) -> (lon, lat) in radians.

    Vectorised transcription of CubedSphereRLLFromABP. panel is 1-6 and may be
    an array. Inverse of abp_from_rll.
    """
    alpha = np.asarray(alpha, float); beta = np.asarray(beta, float)
    panel = np.asarray(panel, int)
    a1, b1 = np.tan(alpha), np.tan(beta)
    sz = (1.0 + a1*a1 + b1*b1)**-0.5
    sx, sy = sz*a1, sz*b1
    xx = np.zeros_like(sx); yy = np.zeros_like(sx); zz = np.zeros_like(sx)
    for p, (fx, fy, fz) in {6: (-sy,  sx,  sz), 5: ( sy,  sx, -sz),
                            1: ( sz,  sx,  sy), 3: (-sz, -sx,  sy),
                            2: (-sx,  sz,  sy), 4: ( sx, -sz,  sy)}.items():
        k = panel == p
        if k.any():
            xx[k], yy[k], zz[k] = fx[k], fy[k], fz[k]
    lat = np.arcsin(np.clip(zz, -1.0, 1.0))
    lon = np.mod(np.arctan2(yy, xx) + rotate_cube, 2.0*np.pi)
    return lon, lat


def lonlat_from_cube_index(x, y, panel, ncube, rotate_cube=0.0):
    """1-based cube index -> (lon, lat) radians. Inverse of
    cube_index_from_lonlat."""
    half = np.pi/2.0
    alpha = (np.asarray(x, float) - 0.5 - (ncube//2)) * half / ncube
    beta  = (np.asarray(y, float) - 0.5 - (ncube//2)) * half / ncube
    return rll_from_abp(alpha, beta, panel, rotate_cube)


def _wrap(d):
    """Longitude difference wrapped into (-pi, pi]."""
    return np.mod(d + np.pi, 2.0*np.pi) - np.pi


def downwind_sign(tg, usrc, vsrc, step=0.05, calm_cos=1.0e-4):
    """Which side of each crest is downwind, in cube index space.

    Returns (ncol, nrdg) of +1, -1 or 0.

      +1  the +n side (t > 0) is downwind
      -1  the -n side (t < 0) is downwind
       0  wind is along-crest or absent; side undefined (|cos| < calm_cos)

    where n = (cos(anglx), -sin(anglx)) is the cross-crest direction in cube
    index space.

    The sign cannot be reasoned out from anglx and angll alone: the map from
    (alpha,beta) to (east,north) is not conformal and its handedness varies by
    panel. So it is measured, in two parts.

    Stepping along the crest in cube space gives its geographic direction
    (ce,cn), whose perpendicular N = (cn,-ce) is the true across-crest
    direction. The wind's across-crest component is wN = wind.N -- this alone
    decides which geographic side is downwind, and is zero for along-crest
    flow.

    Separately, stepping along cube +n gives (de,dn); dN = (de,dn).N says
    whether cube +n points with or against N. The half to paint is
    sign(wN * dN).

    Note the wind must NOT simply be projected onto (de,dn): shear in the map
    means (de,dn) is not perpendicular to the crest, so an along-crest wind
    would give a spurious non-zero projection.
    """
    a = np.radians(tg.anglx)
    ux, uy =  np.sin(a), np.cos(a)      # along-crest, cube space
    nx, ny =  np.cos(a), -np.sin(a)     # cross-crest, cube space

    lon0, lat0 = lonlat_from_cube_index(tg.x, tg.y, tg.panel, tg.ncube)
    cl = np.cos(lat0)

    # geographic image of the ALONG-crest direction
    lonu, latu = lonlat_from_cube_index(tg.x + step*ux, tg.y + step*uy,
                                        tg.panel, tg.ncube)
    ce, cn = cl*_wrap(lonu - lon0), latu - lat0

    # geographic image of the cube-space CROSS-crest direction. Because the
    # map shears, this is NOT perpendicular to (ce,cn) in general, which is
    # why the wind cannot simply be projected onto it.
    lonn, latn = lonlat_from_cube_index(tg.x + step*nx, tg.y + step*ny,
                                        tg.panel, tg.ncube)
    de, dn = cl*_wrap(lonn - lon0), latn - lat0

    # true geographic normal to the crest: rotate (ce,cn) by -90 deg
    Ne, Nn = cn, -ce

    u = np.asarray(usrc, float); v = np.asarray(vsrc, float)
    wN = u*Ne + v*Nn          # wind component ACROSS the crest
    dN = de*Ne + dn*Nn        # does cube +n point along +N, or against it?

    # Undefined when the wind has no across-crest component. Scale-free, so it
    # does not depend on `step` or on ncube.
    speed = np.hypot(u, v)
    crest = np.hypot(ce, cn)
    with np.errstate(invalid="ignore", divide="ignore"):
        cosang = np.where((speed > 0) & (crest > 0), wN/np.maximum(speed*crest, 1e-300), 0.0)

    sgn = np.sign(wN * dN)
    sgn[np.abs(cosang) < calm_cos] = 0
    sgn[(speed <= 0) | (crest <= 0)] = 0
    return sgn.astype(np.int8)


def cube_index_from_lonlat(lon, lat, ncube, rotate_cube=0.0):
    """(lon, lat) radians -> (x, y, panel) in 1-based cube index space.

    Inverse of the forward transform in remapridge2tiles. Returns float
    indices; ncube//2 reproduces Fortran integer division.
    """
    alpha, beta, panel = abp_from_rll(lon, lat, rotate_cube)
    half = np.pi / 2.0
    x = alpha * ncube / half + 0.5 + (ncube // 2)
    y = beta * ncube / half + 0.5 + (ncube // 2)
    return x, y, panel


# ---------------------------------------------------------------------------
@dataclass
class TileGeom:
    """Ridge geometry in cube index space, ready for painting."""
    x: np.ndarray          # (ncol, nrdg) cube alpha index, 1-based
    y: np.ndarray          # (ncol, nrdg) cube beta index, 1-based
    panel: np.ndarray      # (ncol, nrdg) 1-6
    anglx: np.ndarray      # (ncol, nrdg) degrees, clockwise from b-axis
    halflen: np.ndarray    # (ncol, nrdg) along-crest half-length, cube cells
    halfwid: np.ndarray    # (ncol, nrdg) cross-crest half-width, cube cells
    valid: np.ndarray      # (ncol, nrdg) bool
    ncube: int
    grid_length_scale: float
    file: str = ""
    source: str = ""       # how x/y/panel were obtained

    def __repr__(self):
        return (f"TileGeom(ncol={self.x.shape[0]}, nrdg={self.x.shape[1]}, "
                f"ncube={self.ncube}, valid={int(self.valid.sum()):,}, "
                f"source={self.source!r})")


def load_tiles_nc(path, ncube, grid_length_scale, rotate_cube=0.0,
                  length_var="CLNGT", width_var="HWDTH", verbose=False):
    """Read ridge geometry from a CAM-style tiles topography netCDF.

    Parameters
    ----------
    path : str
    ncube : int
        Cube face dimension of the intermediate topo grid (3000 typically).
        NOT on the file, so it must be supplied.
    grid_length_scale : float
        km per cube cell. Written as a global attribute by wrtncdf_rll; for
        unstructured files pass it explicitly (earth_radius*(pi/2)/ncube,
        about 3.336 km at ncube=3000).
    length_var : str
        'CLNGT' is the in-column crest length in km, which is what you almost
        always want. 'CLNGT_RDG' is the whole-ridge length.
    width_var : str
        'HWDTH' is in km on the netCDF (unlike Ridge_tile_map.dat, where it is
        in cube cells).
    """
    import netCDF4 as nc

    with nc.Dataset(path) as d:
        # Fortran nf90_def_var lists dimids fastest-varying first, so a
        # variable declared (ncol,nrdg) in ridge_ana/cube_to_target comes back
        # from netCDF4-python as (nrdg,ncol). Read the dimension names and
        # transpose to a canonical (ncol,nrdg) rather than assuming either.
        def _orient(v):
            a = np.asarray(v[:], dtype=float)
            dims = v.dimensions
            if len(dims) != 2:
                raise ValueError(f"{path}: expected a 2-D variable, got {dims}")
            if dims[0] == "ncol":
                return a
            if dims[1] == "ncol":
                return a.T
            raise ValueError(f"{path}: neither dimension of {dims} is 'ncol'")

        def get(name):
            if name not in d.variables:
                raise KeyError(f"{path}: no variable {name!r}. "
                               f"have: {sorted(d.variables)}")
            return _orient(d.variables[name])

        _dims_of_mxdis = d.variables["MXDIS"].dimensions

        mxdis = get("MXDIS")
        anglx = get("ANGLX")
        length_km = get(length_var)
        width_km = get(width_var)

        have_xy = all(v in d.variables for v in ("XLOID", "YLOID", "PANEL"))
        if have_xy:
            x = get("XLOID"); y = get("YLOID")
            panel = get("PANEL").astype(int)
            source = "XLOID/YLOID/PANEL on file"
        else:
            latc = get("latc"); lonc = get("lonc")
            x, y, panel = cube_index_from_lonlat(lonc, latc, ncube, rotate_cube)
            source = "derived from lonc/latc"

    # clean() zeroed unused slots, so MXDIS>0 marks a real ridge
    valid = mxdis > 0.0
    if not have_xy:
        # lonc/latc of exactly 0 are fill from clean(), not the Gulf of Guinea
        valid &= ~((np.asarray(lonc) == 0.0) & (np.asarray(latc) == 0.0))
    valid &= np.isfinite(x) & np.isfinite(y)
    valid &= (panel >= 1) & (panel <= 6)
    valid &= length_km > 0.0

    halflen = 0.5 * length_km / grid_length_scale
    halfwid = 0.5 * width_km / grid_length_scale

    tg = TileGeom(x=x, y=y, panel=panel, anglx=anglx,
                  halflen=halflen, halfwid=halfwid, valid=valid,
                  ncube=ncube, grid_length_scale=grid_length_scale,
                  file=path, source=source)

    if verbose:
        print(f"loaded {path}")
        print(f"  file dims for MXDIS: {_dims_of_mxdis}  -> canonical (ncol,nrdg)")
        print(f"  ncol={x.shape[0]:,}  nrdg={x.shape[1]}  ncube={ncube}")
        print(f"  cube coords: {source}")
        print(f"  valid tiles: {int(valid.sum()):,} of {valid.size:,} "
              f"({100*valid.mean():.1f}%)")
        if valid.any():
            print(f"  half-length (cells): median {np.median(halflen[valid]):.1f}"
                  f"  max {halflen[valid].max():.1f}")
            print(f"  half-width  (cells): median {np.median(halfwid[valid]):.1f}"
                  f"  max {halfwid[valid].max():.1f}")
    return tg


# ---------------------------------------------------------------------------
MODULATIONS = ("uniform", "linear", "cosine", "exp", "ramp_downwind")


def _modulation(kind, s, t, L, W, decay=2.0):
    """Multiplier in [0,1] as a function of position within the polygon.

    s is along-crest, t across-crest, both in cube cells; L and W are the
    half-extents. Sign of t is as painted, so for side='downwind' the downwind
    half has t of one sign only -- use abs(t) if you want a shape that does not
    care which side it is on.

      'uniform'        1 everywhere (default; original behaviour)
      'linear'         1 at the crest, falling to 0 at |t| = W
      'cosine'         cos(pi/2 * t/W); smooth, zero slope at the crest
      'exp'            exp(-decay * |t|/W); e-folding W/decay
      'ramp_downwind'  0 at the crest rising to 1 at |t| = W, for things that
                       develop away from the ridge rather than at it
    """
    if kind == "uniform":
        return 1.0
    tn = np.abs(t) / W
    if kind == "linear":
        return np.clip(1.0 - tn, 0.0, 1.0)
    if kind == "cosine":
        return np.cos(0.5 * np.pi * np.clip(tn, 0.0, 1.0))
    if kind == "exp":
        return np.exp(-decay * tn)
    if kind == "ramp_downwind":
        return np.clip(tn, 0.0, 1.0)
    raise ValueError(f"modulation must be one of {MODULATIONS} or a callable, "
                     f"got {kind!r}")


def _inside(s, t, L, W, shape, half=0):
    """Inside-test in local crest coordinates.

    half =  0  whole polygon
         = +1  keep only t >= 0  (the +n side)
         = -1  keep only t <= 0
    """
    if shape == "diamond":
        m = np.abs(s) / L + np.abs(t) / W <= 1.0
    elif shape == "rectangle":
        m = (np.abs(s) <= L) & (np.abs(t) <= W)
    elif shape == "triangle":
        # Isoceles triangle whose BASE is the crest itself: vertices
        # (-L,0), (+L,0), (0,+/-W). One half of the diamond, hence always
        # one-sided; with side='downwind' the apex points downwind.
        m = np.abs(s) / L + np.abs(t) / W <= 1.0
        if half == 0:
            half = 1
    elif shape == "wedge_end":
        # The old 'triangle': base ACROSS the crest at one end, tapering to a
        # point at the other. Halving it with side= gives a RIGHT triangle.
        m = (np.abs(s) <= L) & (np.abs(t) <= W * (L - s) / (2.0 * L))
    else:
        raise ValueError(f"shape must be one of {SHAPES}, got {shape!r}")
    if half > 0:
        m &= t >= 0.0
    elif half < 0:
        m &= t <= 0.0
    return m


def tile_polygon(xc, yc, anglx_deg, L, W, shape="diamond"):
    """Polygon vertices in cube index space, for plotting or checking."""
    a = np.radians(anglx_deg)
    u = np.array([np.sin(a),  np.cos(a)])
    n = np.array([np.cos(a), -np.sin(a)])
    c = np.array([xc, yc])
    loc = {"diamond":   [( L, 0.), (0.,  W), (-L, 0.), (0., -W)],
           "rectangle": [( L,  W), (-L,  W), (-L, -W), ( L, -W)],
           "triangle":  [(-L, 0.), ( L, 0.), (0.,  W)],
           "wedge_end": [( L, 0.), (-L,  W), (-L, -W)]}[shape]
    return np.array([c + s * u + t * n for s, t in loc])


def paint_tiles_nc(A, tg, shape="diamond", fill=np.nan, ir_range=None,
                   side="both", usrc=None, vsrc=None, calm_policy="both",
                   modulation="uniform", decay=2.0,
                   min_halfwidth=0.5, length_scale=1.0, width_scale=1.0,
                   dtype=np.float32, return_owner=False, verbose=False):
    """Paint A[i, ir] onto a (6, ncube, ncube) cubed-sphere field.

    Parameters
    ----------
    A : (ncol, nrdg) array
        Quantity to paint. Non-finite entries are skipped.
    tg : TileGeom
        From load_tiles_nc.
    shape : {'diamond', 'rectangle', 'triangle', 'wedge_end'}
        'triangle' is half the diamond: isoceles with its BASE along the crest,
        vertices (-L,0), (+L,0), (0,W). Always one-sided, so with
        side='downwind' the apex points downwind. 'wedge_end' is the old
        'triangle' (base across the crest at one end); halving that with side=
        gives a right triangle.
    fill : float
        Value for cells no tile reached.
    ir_range : None, (lo, hi), slice, or sequence of int
        Restrict painting to a subset of the ridge dimension. A 2-tuple is
        HALF-OPEN and 0-based, matching Python slicing: (0,1) is the dominant
        ridge only, (10,20) is ridges 10..19. Tiles are sorted by importance
        (mxdis*clnin descending) in the topo file, so ir=0 is the largest
        contributor in each column. None uses every ridge.
    side : {'both', 'downwind', 'upwind'}
        Paint the whole polygon, or only the half on one side of the crest.
        'downwind' and 'upwind' require usrc and vsrc.
    usrc, vsrc : (ncol,nrdg) arrays, or (nrdg,ncol)
        Source-layer zonal and meridional wind, i.e. USRC_RDG and VSRC_RDG from
        the GW driver. Used only to decide which side of each crest is
        downwind; magnitude does not matter.
    modulation : str or callable
        Internal structure within each polygon. One of 'uniform' (default),
        'linear', 'cosine', 'exp', 'ramp_downwind' -- see _modulation -- or a
        callable f(s, t, L, W) returning a multiplier array, with s along-crest
        and t across-crest in cube cells. The painted value is A*f, so with a
        non-uniform modulation the "keep the larger" conflict rule compares
        MODULATED values: a strong ridge far away no longer automatically beats
        a weaker one whose crest passes through the cell.
    decay : float
        e-folding parameter for modulation='exp'; the scale is W/decay.
    calm_policy : {'both', 'skip'}
        What to do when the wind is along-crest, so the side is undefined.
        'both' paints the whole polygon, 'skip' paints nothing.
    min_halfwidth : float
        Floor on the cross-crest half-width in cells. Cell centres are tested,
        so without this a sub-cell-wide ridge paints nothing.
    length_scale, width_scale : float
        Multipliers on the extents, for tuning the rendering.
    dtype : numpy dtype
        float32 by default; a (6,3000,3000) float64 field is 0.43 GB.
    return_owner : bool
        Also return the flat tile index (i*nrdg + ir) that won each cell, or -1.
    """
    if shape not in SHAPES:
        raise ValueError(f"shape must be one of {SHAPES}, got {shape!r}")
    if not callable(modulation) and modulation not in MODULATIONS:
        raise ValueError(f"modulation must be a callable or one of {MODULATIONS}, "
                         f"got {modulation!r}")

    A = np.asarray(A, dtype=float)
    if A.shape == tg.x.shape:
        pass
    elif A.shape == tg.x.shape[::-1]:
        # A came straight from a Fortran-written netCDF as (nrdg,ncol)
        A = A.T
        if verbose:
            print(f"  A supplied as {A.shape[::-1]} (nrdg,ncol); transposed to (ncol,nrdg)")
    else:
        raise ValueError(
            f"A has shape {A.shape}; tile geometry is {tg.x.shape} (ncol,nrdg). "
            "If A came from the GW driver output, it is (time,nrdg,ncol) - pass "
            "A[itime] and it will be transposed automatically.")

    ncube = tg.ncube
    nrdg_all = A.shape[1]

    # ridge-index selection
    keep = np.zeros(nrdg_all, dtype=bool)
    if ir_range is None:
        keep[:] = True
        rng_txt = f"all {nrdg_all}"
    elif isinstance(ir_range, slice):
        keep[ir_range] = True
        rng_txt = f"slice {ir_range}"
    elif isinstance(ir_range, tuple) and len(ir_range) == 2:
        lo, hi = ir_range
        if not (0 <= lo <= hi <= nrdg_all):
            raise ValueError(f"ir_range {ir_range} outside [0,{nrdg_all}]")
        keep[lo:hi] = True
        rng_txt = f"[{lo},{hi})"
    else:
        idx = np.asarray(ir_range, dtype=int)
        if idx.min() < 0 or idx.max() >= nrdg_all:
            raise ValueError(f"ir_range indices outside [0,{nrdg_all-1}]")
        keep[idx] = True
        rng_txt = f"{list(idx)}"
    if not keep.any():
        raise ValueError(f"ir_range {ir_range} selects no ridges")

    use = tg.valid & np.isfinite(A) & keep[None, :]

    # which side of each crest to paint
    if side not in ("both", "downwind", "upwind"):
        raise ValueError("side must be 'both', 'downwind' or 'upwind'")
    if side == "both":
        sgn = None
    else:
        if usrc is None or vsrc is None:
            raise ValueError(f"side={side!r} needs usrc and vsrc "
                             "(USRC_RDG / VSRC_RDG from the GW driver)")
        u = np.asarray(usrc, float); v = np.asarray(vsrc, float)
        if u.shape == tg.x.shape[::-1]:
            u = u.T
        if v.shape == tg.x.shape[::-1]:
            v = v.T
        if u.shape != tg.x.shape or v.shape != tg.x.shape:
            raise ValueError(f"usrc/vsrc shapes {u.shape}/{v.shape} do not match "
                             f"tile geometry {tg.x.shape}")
        sgn = downwind_sign(tg, u, v)
        if side == "upwind":
            sgn = -sgn
        ncalm = int((sgn[use] == 0).sum())
        if calm_policy == "skip":
            use &= sgn != 0
        elif calm_policy != "both":
            raise ValueError("calm_policy must be 'both' or 'skip'")
        if verbose:
            print(f"  side={side}: {ncalm:,} tiles with along-crest wind "
                  f"({calm_policy})")

    field = np.full((6, ncube, ncube), -np.inf, dtype=dtype)
    owner = np.full((6, ncube, ncube), -1, dtype=np.int64) if return_owner else None

    ii, rr = np.nonzero(use)
    nrdg = A.shape[1]
    if verbose:
        print(f"painting {len(ii):,} tiles as {shape}s   ridges: {rng_txt}   side: {side}")

    npainted = nclipped = 0
    for i_, r_ in zip(ii, rr):
        L = max(float(tg.halflen[i_, r_]) * length_scale, 0.5)
        W = max(float(tg.halfwid[i_, r_]) * width_scale, min_halfwidth)
        xc = float(tg.x[i_, r_]); yc = float(tg.y[i_, r_])
        a = np.radians(float(tg.anglx[i_, r_]))
        ip = int(tg.panel[i_, r_]) - 1
        val = float(A[i_, r_])

        ca, sa = np.cos(a), np.sin(a)
        ex = abs(L * sa) + abs(W * ca)
        ey = abs(L * ca) + abs(W * sa)
        x0 = int(np.floor(xc - ex)); x1 = int(np.ceil(xc + ex))
        y0 = int(np.floor(yc - ey)); y1 = int(np.ceil(yc + ey))
        if x0 < 1 or x1 > ncube or y0 < 1 or y1 > ncube:
            nclipped += 1
        x0 = max(x0, 1); x1 = min(x1, ncube)
        y0 = max(y0, 1); y1 = min(y1, ncube)
        if x1 < x0 or y1 < y0:
            continue

        dx = np.arange(x0, x1 + 1, dtype=float)[None, :] - xc
        dy = np.arange(y0, y1 + 1, dtype=float)[:, None] - yc
        s = dx * sa + dy * ca
        t = dx * ca - dy * sa
        m = _inside(s, t, L, W, shape, 0 if sgn is None else int(sgn[i_, r_]))
        if not m.any():
            continue

        if callable(modulation):
            vals = val * np.asarray(modulation(s, t, L, W), dtype=float)
        elif modulation == "uniform":
            vals = val
        else:
            vals = val * _modulation(modulation, s, t, L, W, decay)

        sub = field[ip, y0 - 1:y1, x0 - 1:x1]
        if return_owner:
            cand = np.where(m, vals, -np.inf)
            better = m & (cand > sub)
            sub[better] = np.broadcast_to(cand, sub.shape)[better]
            owner[ip, y0 - 1:y1, x0 - 1:x1][better] = i_ * nrdg + r_
        else:
            np.maximum(sub, np.where(m, vals, -np.inf), out=sub)
        npainted += 1

    blank = ~np.isfinite(field)
    field[blank] = fill

    if verbose:
        print(f"  {npainted:,} tiles produced at least one cell")
        print(f"  {nclipped:,} clipped at a panel edge (no cross-panel wrapping)")
        print(f"  {100*(1-blank.mean()):.2f}% of the cubed sphere painted")

    return (field, owner) if return_owner else field


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__.splitlines()[3])
    p.add_argument("file")
    p.add_argument("--ncube", type=int, default=3000)
    p.add_argument("--gls", type=float, default=None,
                   help="km per cube cell (default: earth_radius*(pi/2)/ncube)")
    a = p.parse_args()
    gls = a.gls if a.gls is not None else 6371.0 * (np.pi / 2) / a.ncube
    load_tiles_nc(a.file, a.ncube, gls, verbose=True)
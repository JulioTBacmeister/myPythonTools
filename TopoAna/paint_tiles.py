#!/usr/bin/env python3
"""
paint_tiles.py
==============

Paint a per-tile quantity A(ntarget, nrdg) back onto the high-resolution
cubed sphere, so that ridge-scheme output can be viewed in the space where the
ridges were found.

Each tile is rendered as a convex polygon centred on its crest centroid,
oriented along the crest, with

    along-crest  extent  from clnin_tiles  (in-column crest length, km)
    cross-crest  extent  from hwdth_tiles  (or whatever you pass)

Where two tiles overlap, the cell keeps the LARGER value of A.

Geometry
--------
Everything happens in cube index space on a single panel: x is the alpha index,
y the beta index, both 1-based to match Fortran, and the centroids
xloid_tiles / yloid_tiles are already in those units.

anglx_tiles is "ridge orientation clockwise from the b-axis in the cubed sphere
panel", so the crest direction is

    u = ( sin(anglx), cos(anglx) )        anglx = 0  ->  crest along +y
    n = ( cos(anglx), -sin(anglx) )       cross-crest

A point P is tested in local coordinates s = (P-C).u along the crest and
t = (P-C).n across it, against a shape of half-length L and half-width W:

    'diamond'    |s|/L + |t|/W <= 1          rhombus, tapering to points
    'rectangle'  |s| <= L and |t| <= W
    'triangle'   isoceles: full width W at s = -L, tapering to a point at s = +L

Limitations
-----------
* No cross-panel wrapping. A polygon that would spill onto a neighbouring panel
  is clipped at the edge. Ridges near panel boundaries are therefore
  under-painted. Fixing this needs the cubed-sphere halo machinery.
* Cell centres are tested, so a polygon much narrower than one cube cell paints
  nothing. A minimum half-width is applied (see `min_halfwidth`).

Usage
-----
::

    from rdtilemap import rdtilemap
    from paint_tiles import paint_tiles

    t = rdtilemap('Ridge_tile_map.dat')
    A = taudsw - tauoro              # (ntarget, nrdg) from the GW driver

    field = paint_tiles(A, t, ncube=3000, grid_length_scale=3.3,
                        shape='diamond')
    # field is (6, ncube, ncube), fill where nothing was painted
"""

from __future__ import annotations

import numpy as np

__all__ = ["paint_tiles", "tile_polygon"]


SHAPES = ("diamond", "rectangle", "triangle")


def _inside(s, t, L, W, shape):
    """Boolean mask: is local point (s,t) inside the shape of half-size L,W?"""
    if shape == "diamond":
        return np.abs(s) / L + np.abs(t) / W <= 1.0
    if shape == "rectangle":
        return (np.abs(s) <= L) & (np.abs(t) <= W)
    if shape == "triangle":
        # full width at s=-L, point at s=+L
        frac = (L - s) / (2.0 * L)          # 1 at s=-L, 0 at s=+L
        return (np.abs(s) <= L) & (np.abs(t) <= W * frac)
    raise ValueError(f"shape must be one of {SHAPES}, got {shape!r}")


def tile_polygon(xc, yc, anglx_deg, L, W, shape="diamond"):
    """Vertices of one tile's polygon in cube index space.

    Returned as an (nv, 2) array of (x, y), for plotting or checking. Not used
    by the painter itself, which works with the analytic inside-test.
    """
    a = np.radians(anglx_deg)
    u = np.array([np.sin(a),  np.cos(a)])
    n = np.array([np.cos(a), -np.sin(a)])
    c = np.array([xc, yc])
    if shape == "diamond":
        loc = [( L, 0.0), (0.0,  W), (-L, 0.0), (0.0, -W)]
    elif shape == "rectangle":
        loc = [( L,  W), (-L,  W), (-L, -W), ( L, -W)]
    elif shape == "triangle":
        loc = [( L, 0.0), (-L,  W), (-L, -W)]
    else:
        raise ValueError(f"shape must be one of {SHAPES}, got {shape!r}")
    return np.array([c + s * u + t * n for s, t in loc])


def paint_tiles(A, tiles, ncube, grid_length_scale,
                width=None, shape="diamond",
                fill=np.nan, min_halfwidth=0.5, length_scale=1.0,
                width_scale=1.0, return_owner=False, verbose=False):
    """Paint A(ntarget, nrdg) onto a (6, ncube, ncube) cubed-sphere field.

    Parameters
    ----------
    A : (ntarget, nrdg) array
        Quantity to paint. NaN entries are skipped.
    tiles : TileMap
        From rdtilemap. Needs xloid_tiles/yloid_tiles (as `xloid`/`yloid`),
        `panel`, `anglx`, `clnin` (or `clngt`, see `length_from`) and
        `NumObjects`.
    ncube : int
        Cube face dimension (3000 for the standard intermediate grid).
    grid_length_scale : float
        km per cube cell, from the topo run. Used to convert clnin (km) to
        cells.
    width : None, float, or (ntarget, nrdg) array
        Cross-crest HALF-width in cube cells. None uses tiles.hwdth (which is
        in cube cells already). A scalar applies everywhere.
    shape : {'diamond', 'rectangle', 'triangle'}
    fill : float
        Value for cells no tile reached.
    min_halfwidth : float
        Floor on the cross-crest half-width, in cells. Without this, thin
        ridges paint nothing because no cell centre falls inside.
    length_scale, width_scale : float
        Multipliers on the along- and cross-crest extents, for tuning the
        rendering without changing the underlying data.
    return_owner : bool
        Also return an (6, ncube, ncube) int array holding the flat tile index
        (j*nrdg + ir) that won each cell, or -1.

    Returns
    -------
    field : (6, ncube, ncube) float array
    owner : optional (6, ncube, ncube) int array
    """
    if shape not in SHAPES:
        raise ValueError(f"shape must be one of {SHAPES}, got {shape!r}")

    A = np.asarray(A, dtype=float)
    ntarget, nrdg = A.shape

    xc_all = np.asarray(tiles.xloid, dtype=float)
    yc_all = np.asarray(tiles.yloid, dtype=float)
    pn_all = np.asarray(tiles.panel, dtype=int)
    an_all = np.asarray(tiles.anglx, dtype=float)

    # along-crest length: prefer the in-column length if present
    if hasattr(tiles, "clnin"):
        len_km = np.asarray(tiles.clnin, dtype=float)
    else:
        # clngt_tiles is the whole-ridge length in CUBE CELLS
        len_km = np.asarray(tiles.clngt, dtype=float) * grid_length_scale
    L_all = 0.5 * len_km / grid_length_scale * length_scale     # half-length, cells

    if width is None:
        W_all = 0.5 * np.asarray(tiles.hwdth, dtype=float)
    elif np.isscalar(width):
        W_all = np.full_like(L_all, float(width))
    else:
        W_all = np.asarray(width, dtype=float)
    W_all = W_all * width_scale

    # trim shapes to A in case the tile file has more slots than A
    def _fit(a):
        return a[:, :nrdg] if a.ndim == 2 and a.shape[1] >= nrdg else a

    xc_all, yc_all, pn_all, an_all = map(_fit, (xc_all, yc_all, pn_all, an_all))
    L_all, W_all = map(_fit, (L_all, W_all))

    valid = np.zeros(A.shape, dtype=bool)
    j = np.arange(nrdg)[None, :]
    valid |= j < np.asarray(tiles.NumObjects)[:, None]
    valid &= np.isfinite(A)
    valid &= (pn_all >= 1) & (pn_all <= 6)
    valid &= (xc_all > -999.) & (yc_all > -999.)
    valid &= (an_all > -900.)
    valid &= np.isfinite(L_all) & (L_all > 0)

    field = np.full((6, ncube, ncube), -np.inf, dtype=float)
    owner = np.full((6, ncube, ncube), -1, dtype=np.int64) if return_owner else None

    jj, rr = np.nonzero(valid)
    if verbose:
        print(f"painting {len(jj):,} tiles of {A.size:,} slots "
              f"({100*len(jj)/A.size:.1f}%) as {shape}s")

    npainted = 0
    nclipped = 0
    for j_, r_ in zip(jj, rr):
        L = max(float(L_all[j_, r_]), 0.5)
        W = max(float(W_all[j_, r_]), min_halfwidth)
        xc, yc = float(xc_all[j_, r_]), float(yc_all[j_, r_])
        a = np.radians(float(an_all[j_, r_]))
        ip = int(pn_all[j_, r_]) - 1
        val = float(A[j_, r_])

        # bounding box of the rotated shape, in 1-based cell indices
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

        xs = np.arange(x0, x1 + 1, dtype=float)
        ys = np.arange(y0, y1 + 1, dtype=float)
        dx = xs[None, :] - xc              # (ny, nx)
        dy = ys[:, None] - yc
        s = dx * sa + dy * ca              # along crest
        t = dx * ca - dy * sa              # across crest
        m = _inside(s, t, L, W, shape)
        if not m.any():
            continue

        sub = field[ip, y0 - 1:y1, x0 - 1:x1]
        if return_owner:
            better = m & (val > sub)
            sub[better] = val
            owner[ip, y0 - 1:y1, x0 - 1:x1][better] = j_ * nrdg + r_
        else:
            np.maximum(sub, np.where(m, val, -np.inf), out=sub)
        npainted += 1

    unpainted = ~np.isfinite(field)
    field[unpainted] = fill

    if verbose:
        frac = 100 * (1 - unpainted.mean())
        print(f"  {npainted:,} tiles produced at least one cell")
        print(f"  {nclipped:,} tiles clipped at a panel edge")
        print(f"  {frac:.2f}% of the cubed sphere painted")

    return (field, owner) if return_owner else field
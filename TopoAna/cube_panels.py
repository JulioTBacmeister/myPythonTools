"""
cube_panels.py  --  plot regions of the cubed-sphere topo grid that span one
                    or more cube faces (panels), with no seam at face edges.

    from myPythonTools.TopoAna import cube_panels as CP

    R = CP.CubeRegion(tlons, tlats, lonrange, latrange, MapProj)   # build once
    R.contourf(ax, dev, levels=topo_levels, cmap=cmap, norm=norm)
    R.contour (ax, raw, levels=[0.1, 1e4], colors='black')
    R.set_extent(ax)              # optional: frame the map on the box (needed around poles)
    R.faces                       # 1-based faces used, e.g. [1, 6] for Italy

Arrays are the notebook's (6, nc, nc) reshapes of the topo file: tlons, tlats
in degrees (lon 0..360 or -180..180, either works); fields likewise.  ax must
use the same projection as MapProj.  No matplotlib backend is set here.

How it works
------------
1. Faces are found automatically: any face with a cell inside the lon/lat box.
   Each face gets its own (y, x) index window (one xslice/yslice cannot serve
   two faces: the box covers a different index rectangle on each).
2. The points of all windows are transformed to MAP coordinates and joined in
   one Delaunay triangulation (built once, reused for every field).  Contouring
   on that triangulation runs straight across face edges, so there is no gap
   where two faces meet, and nothing breaks at the 0/360 longitude line or at a
   pole (faces 5 and 6).
3. Triangles much longer than the grid spacing (which only bridge concave
   corners outside the data) are masked, as are triangles touching NaNs.

stride > 1 thins the points (every stride-th cell) for big regions; the
triangulation of ~250k points (Italy at ncube=3000) takes a few seconds.
"""
import numpy as np
import matplotlib.tri as mtri
import cartopy.crs as ccrs

__all__ = ["wrap_lon", "in_box", "panel_slices", "CubeRegion"]

_PC = ccrs.PlateCarree()


def wrap_lon(lon, lon0):
    """Longitudes wrapped into [lon0-180, lon0+180)."""
    return (np.asarray(lon, dtype=float) - lon0 + 180.) % 360. - 180. + lon0


def in_box(lon, lat, lonrange, latrange):
    """Boolean mask of points inside the lon/lat box (either lon convention)."""
    lon0 = 0.5 * (lonrange[0] + lonrange[1])
    lo = wrap_lon(lon, lon0)
    return (lo > lonrange[0]) & (lo <= lonrange[1]) & (lat > latrange[0]) & (lat <= latrange[1])


def panel_slices(tlons, tlats, lonrange, latrange, faces=None, pad=3):
    """
    {face (0-based): (yslice, xslice)} for every face with cells in the box.
    faces : optional iterable of 0-based faces to consider (default all).
    pad   : extra cells around each window so contours reach the map edge.
    """
    nc = tlons.shape[-1]
    faces = range(tlons.shape[0]) if faces is None else np.atleast_1d(faces)
    out = {}
    for p in faces:
        p = int(p)
        m = in_box(tlons[p], tlats[p], lonrange, latrange)
        if not m.any():
            continue
        yy, xx = np.nonzero(m)
        out[p] = (slice(max(yy.min() - pad, 0), min(yy.max() + pad + 1, nc)),
                  slice(max(xx.min() - pad, 0), min(xx.max() + pad + 1, nc)))
    if not out:
        raise ValueError(f"no cube cells in lon {lonrange}, lat {latrange}")
    return out


class CubeRegion:
    """Points of every cube face inside a lon/lat box, triangulated in map coordinates."""

    def __init__(self, tlons, tlats, lonrange, latrange, map_proj,
                 faces=None, pad=3, stride=1, max_edge_factor=4.0, margin=0.15):
        # select a slightly larger box (margin = fraction of the box size) so the
        # corners of the map frame, which lie outside the lon/lat box, are filled
        dlon = margin * (lonrange[1] - lonrange[0]) if (lonrange[1] - lonrange[0]) < 360. else 0.
        dlat = margin * (latrange[1] - latrange[0])
        sel_lon = [lonrange[0] - dlon, lonrange[1] + dlon]
        sel_lat = [max(latrange[0] - dlat, -90.), min(latrange[1] + dlat, 90.)]
        self.slices = panel_slices(tlons, tlats, sel_lon, sel_lat, faces=faces, pad=pad)
        # map-coordinate bounds of the points inside the box itself (for set_extent)
        bx, by = [], []
        for p in self.slices:
            m = in_box(tlons[p], tlats[p], lonrange, latrange)
            if m.any():
                q = map_proj.transform_points(_PC, np.asarray(tlons[p])[m].astype(float),
                                              np.asarray(tlats[p])[m].astype(float))
                f = np.isfinite(q[:, 0]) & np.isfinite(q[:, 1])
                bx.append(q[f, 0]); by.append(q[f, 1])
        bx, by = np.concatenate(bx), np.concatenate(by)
        self.bounds = (bx.min(), bx.max(), by.min(), by.max())
        self.faces = [p + 1 for p in self.slices]
        self.stride = int(stride)
        lon = self._gather(tlons)
        lat = self._gather(tlats)
        pts = map_proj.transform_points(_PC, lon, lat)
        ok = np.isfinite(pts[:, 0]) & np.isfinite(pts[:, 1])   # e.g. far side of an Orthographic globe
        self._keep = ok
        x, y = pts[ok, 0], pts[ok, 1]
        self.tri = mtri.Triangulation(x, y)
        # mask triangles that bridge concave corners outside the data
        t = self.tri.triangles
        e = np.max(np.stack([np.hypot(x[t[:, i]] - x[t[:, j]], y[t[:, i]] - y[t[:, j]])
                             for i, j in ((0, 1), (1, 2), (2, 0))]), axis=0)
        self._long = e > max_edge_factor * np.median(e)
        self.tri.set_mask(self._long)
        self.npoints = x.size

    def _gather(self, field):
        s = self.stride
        parts = [np.asarray(field[p])[ys, xs][::s, ::s].ravel() for p, (ys, xs) in self.slices.items()]
        return np.concatenate(parts).astype(float)

    def values(self, field):
        """field (6, nc, nc) -> 1-D values at the triangulation points."""
        return self._gather(field)[self._keep]

    def _call(self, method, ax, field, **kw):
        kw.pop('transform', None)                 # points are already in map coordinates
        z = self.values(field)
        bad = ~np.isfinite(z)
        mask = self._long | bad[self.tri.triangles].any(axis=1)
        self.tri.set_mask(mask)
        try:
            return getattr(ax, method)(self.tri, np.where(bad, 0., z), **kw)
        finally:
            self.tri.set_mask(self._long)

    def set_extent(self, ax):
        """Frame the map on the box's points (works around poles, where
        cartopy's set_extent with a 0-360 longitude range collapses)."""
        x0, x1, y0, y1 = self.bounds
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)

    def contourf(self, ax, field, **kw):
        """ax.tricontourf over all faces in the region. Returns the ContourSet."""
        return self._call('tricontourf', ax, field, **kw)

    def contour(self, ax, field, **kw):
        """ax.tricontour over all faces in the region. Returns the ContourSet."""
        return self._call('tricontour', ax, field, **kw)
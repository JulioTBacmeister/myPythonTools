import numpy as np
from scipy.spatial import Delaunay


__all__ = ["SE_WRF_Mapper"]


def _wrap_lon(lon):
    return (lon + 180.0) % 360.0 - 180.0


def _project(lon, lat, lon0=-97.9, lat0=39.1):
    """
    Simple local CONUS projection.

    Sufficient for constructing a triangulation for illustrative
    interpolation. Coordinates are proportional to distance.
    """
    x = (_wrap_lon(lon - lon0)) * np.cos(np.deg2rad(lat0))
    y = lat - lat0

    return x, y


class SE_WRF_Mapper:
    """
    Fast piecewise-linear interpolation:

        CAM-SE unstructured grid --> WRF CONUS404 grid

    The expensive triangulation and interpolation weights are
    constructed once during initialization.
    """

    def __init__(
        self,
        se_lon,
        se_lat,
        wrf_lon,
        wrf_lat,
        degree_pad=3.0,
        chunk_size=200000,
        lon0=-97.9,            # <-- PATCH: centre of local projection
        lat0=39.1,             # <-- PATCH: (defaults = CONUS, as before)
    ):

        self.se_lon = _wrap_lon(np.asarray(se_lon))
        self.se_lat = np.asarray(se_lat)

        self.wrf_lon = _wrap_lon(np.asarray(wrf_lon))
        self.wrf_lat = np.asarray(wrf_lat)

        self.wrf_shape = self.wrf_lon.shape

        if self.se_lon.ndim != 1:
            raise ValueError("se_lon and se_lat should be 1-D arrays")

        if self.se_lon.shape != self.se_lat.shape:
            raise ValueError("se_lon and se_lat must have the same shape")

        if self.wrf_lon.shape != self.wrf_lat.shape:
            raise ValueError("wrf_lon and wrf_lat must have the same shape")

        #
        # ---------------------------------------------
        # Keep only CAM-SE points near the WRF domain
        # ---------------------------------------------
        #

        lat0 = self.wrf_lat.min() - degree_pad
        lat1 = self.wrf_lat.max() + degree_pad

        lon0 = self.wrf_lon.min() - degree_pad
        lon1 = self.wrf_lon.max() + degree_pad

        mask = (
            (self.se_lat >= lat0)
            & (self.se_lat <= lat1)
            & (self.se_lon >= lon0)
            & (self.se_lon <= lon1)
        )

        self.src_index = np.flatnonzero(mask)

        print(
            f"Using {len(self.src_index)} of "
            f"{len(self.se_lon)} CAM-SE columns"
        )

        #
        # ---------------------------------------------
        # Project source points into local x/y space
        # ---------------------------------------------
        #

        sx, sy = _project(
            self.se_lon[self.src_index],
            self.se_lat[self.src_index],
            lon0=lon0, lat0=lat0,      # <-- PATCH
        )

        src_xy = np.column_stack((sx, sy))

        #
        # ---------------------------------------------
        # Construct Delaunay triangulation
        # ---------------------------------------------
        #

        print("Constructing CAM-SE triangulation ...")

        tri = Delaunay(src_xy)

        #
        # Project WRF coordinates.
        #

        wx, wy = _project(
            self.wrf_lon,
            self.wrf_lat,
            lon0=lon0, lat0=lat0,      # <-- PATCH
        )

        wx = wx.ravel()
        wy = wy.ravel()

        ndst = wx.size

        #
        # For every WRF point we save:
        #
        #     3 source CAM-SE indices
        #     3 interpolation weights
        #
        # After this, the triangulation is no longer needed.
        #

        self.vertices = np.full(
            (ndst, 3),
            -1,
            dtype=np.int32,
        )

        self.weights = np.zeros(
            (ndst, 3),
            dtype=np.float32,
        )

        print(
            f"Generating interpolation weights for "
            f"{ndst} WRF points ..."
        )

        #
        # Work in chunks to keep temporary memory small.
        #

        for n0 in range(0, ndst, chunk_size):

            n1 = min(n0 + chunk_size, ndst)

            xy = np.column_stack((
                wx[n0:n1],
                wy[n0:n1],
            ))

            simplex = tri.find_simplex(xy)

            good = simplex >= 0

            if not np.any(good):
                continue

            sg = simplex[good]
            xyg = xy[good]

            #
            # Convert x/y into barycentric coordinates.
            #
            # First two weights:
            #

            delta = xyg - tri.transform[sg, 2]

            bary12 = np.einsum(
                "nij,nj->ni",
                tri.transform[sg, :2],
                delta,
            )

            #
            # Third weight:
            #

            bary3 = 1.0 - bary12.sum(axis=1)

            bary = np.column_stack((
                bary12,
                bary3,
            ))

            #
            # Three vertices of each triangle.
            #

            verts_local = tri.simplices[sg]

            #
            # Convert local subset indices back to
            # original global CAM-SE ncol indices.
            #

            verts_global = self.src_index[verts_local]

            ii = np.flatnonzero(good) + n0

            self.vertices[ii, :] = verts_global
            self.weights[ii, :] = bary.astype(np.float32)

        #
        # Record mapped WRF points.
        #

        self.valid = self.vertices[:, 0] >= 0

        print(
            f"Mapped {self.valid.sum()} of "
            f"{ndst} WRF points"
        )

        #
        # Triangulation itself can now disappear.
        #

        del tri


    def se_to_wrf(self, field):
        """
        Interpolate field(ncol) onto the CONUS404 grid.
        """

        field = np.asarray(field)

        if field.ndim != 1:
            raise ValueError(
                "field should currently have shape (ncol,)"
            )

        if field.shape[0] != self.se_lon.shape[0]:
            raise ValueError(
                "field length does not match SE grid"
            )

        out = np.full(
            self.vertices.shape[0],
            np.nan,
            dtype=np.result_type(field.dtype, np.float32),
        )

        v = self.vertices[self.valid]
        w = self.weights[self.valid]

        #
        # Three-point linear interpolation.
        #

        out[self.valid] = (
            w[:, 0] * field[v[:, 0]]
            + w[:, 1] * field[v[:, 1]]
            + w[:, 2] * field[v[:, 2]]
        )

        return out.reshape(self.wrf_shape)
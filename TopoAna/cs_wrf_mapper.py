import numpy as np
import esmpy


__all__ = ["CS_WRF_Mapper"]

def _wrap_lon(lon):
    return (lon + 180.0) % 360.0 - 180.0


def _make_grid(lon, lat):
    ny, nx = lon.shape

    grid = esmpy.Grid(
        np.array([nx, ny]),
        staggerloc=esmpy.StaggerLoc.CENTER,
        coord_sys=esmpy.CoordSys.SPH_DEG,
    )

    grid.get_coords(0)[...] = lon.T
    grid.get_coords(1)[...] = lat.T

    return grid


def _bounding_slice(mask, halo=2):
    rows = np.flatnonzero(np.any(mask, axis=1))
    cols = np.flatnonzero(np.any(mask, axis=0))

    if len(rows) == 0:
        return None

    j0 = max(rows[0] - halo, 0)
    j1 = min(rows[-1] + halo + 1, mask.shape[0])

    i0 = max(cols[0] - halo, 0)
    i1 = min(cols[-1] + halo + 1, mask.shape[1])

    return np.s_[j0:j1, i0:i1]


class CS_WRF_Mapper:

    def __init__(self, cs_lon, cs_lat, wrf_lon, wrf_lat,
                 degree_pad=1.0, halo=2):

        self.cs_lon = _wrap_lon(cs_lon)
        self.cs_lat = cs_lat

        self.wrf_lon = _wrap_lon(wrf_lon)
        self.wrf_lat = wrf_lat

        self.nface = cs_lon.shape[0]

        lat0 = self.wrf_lat.min() - degree_pad
        lat1 = self.wrf_lat.max() + degree_pad
        lon0 = self.wrf_lon.min() - degree_pad
        lon1 = self.wrf_lon.max() + degree_pad

        self.forward = []
        self.reverse = []

        for f in range(self.nface):

            lonf = self.cs_lon[f]
            latf = self.cs_lat[f]

            mask = (
                (latf >= lat0) & (latf <= lat1) &
                (lonf >= lon0) & (lonf <= lon1)
            )

            css = _bounding_slice(mask, halo=halo)

            if css is None:
                print(f"face {f}: no CONUS overlap")
                continue

            cslo = lonf[css]
            csla = latf[css]

            wmask = (
                (self.wrf_lat >= csla.min()) &
                (self.wrf_lat <= csla.max()) &
                (self.wrf_lon >= cslo.min()) &
                (self.wrf_lon <= cslo.max())
            )

            wrfs = _bounding_slice(wmask, halo=halo)

            if wrfs is None:
                continue

            wrflo = self.wrf_lon[wrfs]
            wrfla = self.wrf_lat[wrfs]

            print(
                f"face {f}: "
                f"CS {cslo.shape} --> WRF {wrflo.shape}"
            )

            # -----------------------------
            # CS --> WRF
            # -----------------------------

            csgrid = _make_grid(cslo, csla)
            wrfgrid = _make_grid(wrflo, wrfla)

            csfield = esmpy.Field(csgrid)
            wrffield = esmpy.Field(wrfgrid)

            regrid_fw = esmpy.Regrid(
                csfield,
                wrffield,
                regrid_method=esmpy.RegridMethod.BILINEAR,
                unmapped_action=esmpy.UnmappedAction.IGNORE,
            )

            csfield.data[...] = 1.0
            wrffield.data[...] = 0.0

            regrid_fw(
                csfield,
                wrffield,
                zero_region=esmpy.Region.SELECT,
            )

            valid_fw = wrffield.data.T > 0.5

            self.forward.append({
                "face": f,
                "cs_slice": css,
                "wrf_slice": wrfs,
                "src": csfield,
                "dst": wrffield,
                "regrid": regrid_fw,
                "valid": valid_fw,
            })

            # -----------------------------
            # WRF --> CS
            # -----------------------------

            wrffield_r = esmpy.Field(wrfgrid)
            csfield_r = esmpy.Field(csgrid)

            regrid_bw = esmpy.Regrid(
                wrffield_r,
                csfield_r,
                regrid_method=esmpy.RegridMethod.BILINEAR,
                unmapped_action=esmpy.UnmappedAction.IGNORE,
            )

            wrffield_r.data[...] = 1.0
            csfield_r.data[...] = 0.0

            regrid_bw(
                wrffield_r,
                csfield_r,
                zero_region=esmpy.Region.SELECT,
            )

            valid_bw = csfield_r.data.T > 0.5

            self.reverse.append({
                "face": f,
                "cs_slice": css,
                "wrf_slice": wrfs,
                "src": wrffield_r,
                "dst": csfield_r,
                "regrid": regrid_bw,
                "valid": valid_bw,
            })


    def cs_to_wrf(self, topo):

        out = np.zeros(self.wrf_lon.shape)
        count = np.zeros(self.wrf_lon.shape, dtype=np.int8)

        for M in self.forward:

            f = M["face"]
            css = M["cs_slice"]
            wrfs = M["wrf_slice"]

            src = M["src"]
            dst = M["dst"]

            src.data[...] = topo[f][css].T
            dst.data[...] = 0.0

            M["regrid"](
                src,
                dst,
                zero_region=esmpy.Region.SELECT,
            )

            tmp = dst.data.T
            ok = M["valid"]

            oo = out[wrfs]
            cc = count[wrfs]

            oo[ok] += tmp[ok]
            cc[ok] += 1

        good = count > 0

        out[good] /= count[good]
        out[~good] = np.nan

        return out


    def wrf_to_cs(self, field_wrf):

        out = np.full(self.cs_lon.shape, np.nan)

        for M in self.reverse:

            f = M["face"]
            css = M["cs_slice"]
            wrfs = M["wrf_slice"]

            src = M["src"]
            dst = M["dst"]

            src.data[...] = field_wrf[wrfs].T
            dst.data[...] = 0.0

            M["regrid"](
                src,
                dst,
                zero_region=esmpy.Region.SELECT,
            )

            tmp = dst.data.T
            ok = M["valid"]

            oo = out[f][css]
            oo[ok] = tmp[ok]

        return out
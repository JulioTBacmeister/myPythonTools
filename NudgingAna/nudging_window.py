"""
nudging_window.py

Python replica of the spatial "nudging window" logic in CESM/CAM's
nudging.F90 (Patrick Callaghan's module, as used in CAM6/7).

Faithfully reproduces:
  * nudging_readnl : the Invert -> (lo, hi) mapping
  * nudging_init   : the Nudge_Hwin_max / Nudge_Hwin_min normalization
                     constants evaluated at the window center and the
                     four "far corner" combinations
  * nudging_set_profile : the per-column horizontal Heaviside window
                     (product of four tanh transitions, rescaled by
                     max/min, mapped through lo/hi), the raw vertical
                     window in level-index space (with its constant-
                     window short circuits), and the final product
                     Wprof = Hcoef * Vwin

All arithmetic is float64, matching r8. Vectorized over lat/lon arrays
(e.g., MPAS cell centers) via numpy broadcasting.

Conventions (same as F90):
  * lon0 in [0, 360), lat0 in [-90, 90]
  * input lon may be anything; the code recenters lon-lon0 into
    (-180, 180] exactly as the F90 does
  * vertical window is in MODEL LEVEL INDEX space, k = 1..nlev
    (top-of-model k=1 in CAM)
"""

from dataclasses import dataclass
import numpy as np


@dataclass
class HWinParams:
    """Horizontal window namelist parameters (CAM defaults)."""
    lat0: float = 0.0
    lon0: float = 180.0
    latWidth: float = 9999.0
    lonWidth: float = 9999.0
    latDelta: float = 1.0
    lonDelta: float = 1.0
    invert: bool = False

    @property
    def lo(self) -> float:
        return 1.0 if self.invert else 0.0

    @property
    def hi(self) -> float:
        return 0.0 if self.invert else 1.0


@dataclass
class VWinParams:
    """Vertical window namelist parameters (CAM defaults for nlev levels)."""
    nlev: int = 93
    Lindex: float = 0.0       # LO transition model index
    Hindex: float = None      # HI transition model index (default nlev+1)
    Ldelta: float = 0.001
    Hdelta: float = 0.001
    invert: bool = False

    def __post_init__(self):
        if self.Hindex is None:
            self.Hindex = float(self.nlev + 1)

    @property
    def lo(self) -> float:
        return 1.0 if self.invert else 0.0

    @property
    def hi(self) -> float:
        return 0.0 if self.invert else 1.0


def _tanh_step(x):
    """(1 + tanh(x)) / 2 -- the parameterized Heaviside used throughout."""
    return (1.0 + np.tanh(x)) / 2.0


# ----------------------------------------------------------------------
# nudging_init : window normalization constants
# ----------------------------------------------------------------------
def hwin_norm(p: HWinParams):
    """
    Replicates the Nudge_Hwin_max / Nudge_Hwin_min computation in
    nudging_init. Returns (Hwin_max, Hwin_min, lonWidthH, latWidthH).

    max is the raw window evaluated at the window center (relative
    coords (0,0)); min is the smallest of the four far-corner values
    at relative lon = +/-180 and relative lat = (+90-lat0), (-90-lat0).
    """
    lonp, lon0r, lonn = 180.0, 0.0, -180.0
    latp = 90.0 - p.lat0
    lat0r = 0.0
    latn = -90.0 - p.lat0

    lonWH = p.lonWidth / 2.0
    latWH = p.latWidth / 2.0

    def quad(lonx, latx):
        return (_tanh_step((lonWH + lonx) / p.lonDelta)
                * _tanh_step((lonWH - lonx) / p.lonDelta)
                * _tanh_step((latWH + latx) / p.latDelta)
                * _tanh_step((latWH - latx) / p.latDelta))

    hmax = quad(lon0r, lat0r)
    hmin = min(quad(lonp, latn), quad(lonp, latp),
               quad(lonn, latn), quad(lonn, latp))
    return hmax, hmin, lonWH, latWH


# ----------------------------------------------------------------------
# nudging_set_profile : horizontal part
# ----------------------------------------------------------------------
def horizontal_window(lat, lon, p: HWinParams):
    """
    Horizontal window coefficient Hcoef at (lat, lon) in degrees.
    lat/lon may be scalars or numpy arrays (broadcast together).
    Replicates the Nudge_prof==2 horizontal branch of
    nudging_set_profile, including the constant-window short circuit.
    """
    lat = np.asarray(lat, dtype=np.float64)
    lon = np.asarray(lon, dtype=np.float64)

    hmax, hmin, lonWH, latWH = hwin_norm(p)

    if hmax <= hmin:
        # Constant horizontal window: Hcoef = max(lo, hi)
        return np.full(np.broadcast(lat, lon).shape,
                       max(p.lo, p.hi), dtype=np.float64)

    # Recenter relative to window center, wrap lon into (-180, 180]
    latx = lat - p.lat0
    lonx = lon - p.lon0
    lonx = np.where(lonx > 180.0, lonx - 360.0, lonx)
    lonx = np.where(lonx <= -180.0, lonx + 360.0, lonx)

    # RAW window value: product of four tanh transitions
    hcoef = (_tanh_step((lonWH + lonx) / p.lonDelta)
             * _tanh_step((lonWH - lonx) / p.lonDelta)
             * _tanh_step((latWH + latx) / p.latDelta)
             * _tanh_step((latWH - latx) / p.latDelta))

    # Rescale to [0,1] using the init-time max/min, then map lo<->hi
    hcoef = (hcoef - hmin) / (hmax - hmin)
    hcoef = (1.0 - hcoef) * p.lo + hcoef * p.hi
    return hcoef


# ----------------------------------------------------------------------
# nudging_set_profile : vertical part
# ----------------------------------------------------------------------
def vertical_window(p: VWinParams):
    """
    Vertical window Wprof(k), k = 1..nlev, replicating the raw tanh
    profile, and the two-way constant-window short circuit:
      if Vmax <= Vmin  OR  (Hindex >= nlev+1 and Lindex <= 0):
          Wprof(:) = max(lo, hi)
    """
    k = np.arange(1, p.nlev + 1, dtype=np.float64)
    lev_lo = (k - p.Lindex) / p.Ldelta
    lev_hi = (p.Hindex - k) / p.Hdelta
    wprof = _tanh_step(lev_lo) * _tanh_step(lev_hi)

    vmax, vmin = wprof.max(), wprof.min()
    if (vmax <= vmin) or (p.Hindex >= p.nlev + 1 and p.Lindex <= 0):
        return np.full(p.nlev, max(p.lo, p.hi), dtype=np.float64)

    wprof = (wprof - vmin) / (vmax - vmin)
    wprof = p.lo + wprof * (p.hi - p.lo)
    return wprof


# ----------------------------------------------------------------------
# Full 3D window (what multiplies the nudging coefficient)
# ----------------------------------------------------------------------
def window3d(lat, lon, hp: HWinParams, vp: VWinParams, prof: int = 2):
    """
    Full spatial window, shape lat/lon broadcast shape + (nlev,).
    prof: 0 = off (zeros), 1 = uniform (ones), 2 = Heaviside window.
    Matches nudging_set_profile's Wprof(:) = Hcoef * Vwin(:).
    """
    shape = np.broadcast(np.asarray(lat), np.asarray(lon)).shape
    if prof == 0:
        return np.zeros(shape + (vp.nlev,))
    if prof == 1:
        return np.ones(shape + (vp.nlev,))
    if prof != 2:
        raise ValueError(f"Unknown Nudge_prof value: {prof}")
    h = horizontal_window(lat, lon, hp)
    v = vertical_window(vp)
    return np.multiply.outer(np.asarray(h), v)


# ----------------------------------------------------------------------
if __name__ == "__main__":
    # ------------------------------------------------------------------
    # 1) Julio's actual namelist (ERA5 ne30pg3 L93 UV nudging)
    # ------------------------------------------------------------------
    hp = HWinParams(lat0=0.0, lon0=180.0,
                    latWidth=9999.0, lonWidth=9999.0,
                    latDelta=1.0, lonDelta=1.0, invert=False)
    vp = VWinParams(nlev=93, Lindex=0.0, Hindex=0.0,
                    Ldelta=0.001, Hdelta=0.001, invert=True)

    hmax, hmin, _, _ = hwin_norm(hp)
    print("=== Namelist as given ===")
    print(f"Hwin_max = {hmax:.15f}")
    print(f"Hwin_min = {hmin:.15f}")
    print(f"constant-horizontal branch taken: {hmax <= hmin}")

    lat = np.array([-90.0, -45.0, 0.0, 37.0, 90.0])
    lon = np.array([0.0, 90.0, 180.0, 251.3, 359.9])
    H = horizontal_window(lat, lon, hp)
    print("Hcoef at test points:", H)

    V = vertical_window(vp)
    print("Vwin: min,max =", V.min(), V.max(), " (nlev=93)")

    W = window3d(lat, lon, hp, vp)
    print("window3d shape:", W.shape, " unique values:", np.unique(W))

    # effective inverse timescale, matching Nudge_Xtau scaling
    nudge_step = 86400 // 24
    ucoef = 0.16667
    print(f"tau^-1 = Ucoef/Nudge_Step = {ucoef/nudge_step:.3e} s^-1 "
          f"(e-folding {nudge_step/ucoef/3600.0:.2f} h)")

    # ------------------------------------------------------------------
    # 2) A non-degenerate check: the doc-comment channel example
    #    (lat0=0, latWidth=30, latDelta=5, lon constant)
    # ------------------------------------------------------------------
    print("\n=== Channel example (latWidth=30, latDelta=5) ===")
    hp2 = HWinParams(lat0=0.0, lon0=180.0,
                     latWidth=30.0, lonWidth=999.0,
                     latDelta=5.0, lonDelta=1.0)
    hmax2, hmin2, _, _ = hwin_norm(hp2)
    print(f"Hwin_max = {hmax2:.15f}, Hwin_min = {hmin2:.6e}")
    lats = np.array([0.0, 10.0, 15.0, 20.0, 30.0, 60.0, 90.0])
    print("lat :", lats)
    print("Hcoef:", np.round(horizontal_window(lats, 180.0, hp2), 6))
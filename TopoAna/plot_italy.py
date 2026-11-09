# %% ============================================================
# CELL 0 - imports (unchanged from PlotNov28-2016, minus WRF bits)
# ============================================================
import os, sys, glob, importlib
workdir_ = '/glade/work/juliob'
if workdir_ not in sys.path:
    sys.path.append(workdir_)

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pathlib import Path

from myPythonTools.TopoAna import Plots as Po
from myPythonTools.TopoAna import rdtilemap as Ti
import paint_tiles_nc as Pti
from regrid_to_cube import read_cube_target_map, regrid_to_cube
from cs_wrf_mapper import CS_WRF_Mapper
from se_wrf_mapper import SE_WRF_Mapper
importlib.reload(Pti)

# %% ============================================================
# CELL 1 - case, one timestep to start with
# ============================================================
caseN   = 'ERA5-era5-2022-07-Rw050-Friuli'
arch    = f'/glade/derecho/scratch/juliob/archive/GW_UnitTest/{caseN}'

def step_file(date, secs):
    """date='2022-07-20', secs in {0,21600,43200,64800}"""
    return f'{arch}/{caseN}.h.{date}-{secs:05d}.nc'

fileN = step_file('2022-07-20', 0)          # first guess at bora peak
X = xr.open_dataset(fileN)
print(X.bnd_topo)

# %% ============================================================
# CELL 2 - per-ridge quantities -> displacement, u', T'  (your code)
# ============================================================
taudsw = X.TAUDSW.values
tauoro = X.TAUORO.values
usrc   = X.USRC_RDG.values
vsrc   = X.VSRC_RDG.values
nsrc   = X.NSRC_RDG.values
ubsrc  = X.UBMSRC_RDG.values
hwdth  = 1000.*X.HWDTH_RDG.values
mxd    = X.MXDIS_RDG.values
lonc   = X.lon.values
latc   = X.lat.values

rho  = 1.0
Tsrc = 280.
delprm     = np.sqrt( hwdth*taudsw/(rho*nsrc*ubsrc) )
delprm_oro = np.sqrt( hwdth*tauoro/(rho*nsrc*ubsrc) )
uprm       = nsrc*delprm
uprm_dsw   = nsrc*(delprm - delprm_oro)
warmy      = delprm * nsrc**2 * Tsrc / 9.8

# %% ============================================================
# CELL 3 - paint onto the 3000x3000x6 cube (your settings)
# ============================================================
tg = Pti.load_tiles_nc(X.bnd_topo, ncube=3000, grid_length_scale=3.3)

kw = dict(shape='rectangle', width_scale=1.0, verbose=True, ir_range=(0,43),
          side='downwind', usrc=usrc[0,:,:], vsrc=vsrc[0,:,:])
field_uprm     = np.nan_to_num(Pti.paint_tiles_nc(uprm[0,:,:],     tg, modulation='cosine',        **kw))
field_uprm_dsw = np.nan_to_num(Pti.paint_tiles_nc(uprm_dsw[0,:,:], tg, modulation='ramp_downwind', **kw))
field_warm     = np.nan_to_num(Pti.paint_tiles_nc(warmy[0,:,:],    tg, modulation='ramp_downwind', **kw))

# %% ============================================================
# CELL 4 - 3000x3000x6 topography (your code)
# ============================================================
topo_dir = str(Path(X.bnd_topo).parent)
FileN = glob.glob(f"{topo_dir}/topo*.nc")[-1]
print(FileN)
topo  = xr.open_dataset(FileN)
nc    = int(np.sqrt(topo['terr_dev'].size/6))
raw   = np.reshape((topo['terr_dev']+topo['terr_sm']).values, (6,nc,nc))
dev   = np.reshape(topo['terr_dev'].values, (6,nc,nc))
tlats = np.reshape(topo['lat'].values, (6,nc,nc))
tlons = np.reshape(topo['lon'].values, (6,nc,nc))

# %% ============================================================
# CELL 5 - NEW: a regular lat-lon plotting grid over NE Italy
#          (this replaces the WRF CONUS404 grid as the common target)
# ============================================================
# Box: Alps + Friuli + Carso + N Adriatic + Slovenian/Croatian coast.
# The Carso fire area is ~13.5-13.9E, 45.7-45.95N.  Trieste 13.77E 45.65N.
lon0, lon1 = 10.5, 16.0
lat0, lat1 = 44.0, 47.5
dlat = 0.010                     # ~1.1 km
dlon = dlat/np.cos(np.deg2rad(0.5*(lat0+lat1)))   # ~1.1 km at 45.75N

glon_1d = np.arange(lon0, lon1+1e-9, dlon)
glat_1d = np.arange(lat0, lat1+1e-9, dlat)
glon, glat = np.meshgrid(glon_1d, glat_1d)      # shape (ny, nx), like XLONG/XLAT
print('plot grid', glon.shape)

# %% ============================================================
# CELL 6 - NEW: mappers onto the plotting grid (one-off, ~ a minute)
# ============================================================
# Cube -> grid.  Loops over all 6 faces; Friuli straddles panel 1 (top rows)
# and panel 6 (north-polar).  The "no CONUS overlap" messages for the other
# four faces are just the module's old wording - harmless.
mapper_cs = CS_WRF_Mapper(tlons, tlats, glon, glat)

# SE -> grid.  Needs the local-projection centre moved from CONUS to Italy;
# see the small patch to se_wrf_mapper.py at the end of this file.
mapper_se = SE_WRF_Mapper(lonc, latc, glon, glat,
                          lon0=0.5*(lon0+lon1), lat0=0.5*(lat0+lat1))

# %% ============================================================
# CELL 7 - NEW: seam check.  Map the cube topography and look for a NaN
#          stripe along the panel-1 / panel-6 boundary (~44.2-44.6N).
# ============================================================
topo_grid = mapper_cs.cs_to_wrf(raw)
dev_grid  = mapper_cs.cs_to_wrf(dev)
print('unmapped plot-grid points:', np.isnan(topo_grid).sum(), 'of', topo_grid.size)

fig, ax = plt.subplots(figsize=(9,6), subplot_kw={'projection': ccrs.PlateCarree()})
ax.pcolormesh(glon, glat, np.isnan(topo_grid), cmap='Reds', transform=ccrs.PlateCarree())
ax.contour(glon, glat, np.nan_to_num(topo_grid), levels=[500,1000,2000], colors='k', linewidths=0.4, transform=ccrs.PlateCarree())
ax.coastlines('10m'); ax.add_feature(cfeature.BORDERS, linewidth=0.5)
ax.set_title('red = unmapped (seam/edge) ; contours = cube topo mapped to plot grid')

# %% ============================================================
# CELL 8 - NEW: map the fields for this timestep
# ============================================================
klev = 92                                         # same near-surface level as before
cwsp = np.sqrt(X.U.values**2 + X.V.values**2)
cwsp_grid = mapper_se.se_to_wrf(cwsp[0,klev,:])
te_grid   = mapper_se.se_to_wrf(X.T.values[0,klev,:])

uprm_grid     = np.nan_to_num(mapper_cs.cs_to_wrf(field_uprm))
uprm_dsw_grid = np.nan_to_num(mapper_cs.cs_to_wrf(field_uprm_dsw))
warm_grid     = np.nan_to_num(mapper_cs.cs_to_wrf(field_warm))

# %% ============================================================
# CELL 9 - NEW: the two-panel figure (Gatlinburg layout minus WRF)
# ============================================================
enh = 1.5                                   # same enhancement factor used for Fig. 2

def plot_friuli(cwsp_grid, uprm_grid, title_date,
                extent=(lon0, lon1, lat0, lat1), wsplvs=np.linspace(0,30,41),
                cmapN='OrRd', enh=enh, topo_grid=topo_grid, carso=True):
    data_crs = ccrs.PlateCarree()
    proj = ccrs.LambertConformal(central_longitude=13.3, central_latitude=45.8,
                                 standard_parallels=(44., 47.))
    fields = [cwsp_grid, cwsp_grid + enh*uprm_grid]
    titles = ['Mapped CAM-SE wind speed (ERA5, ne30pg3)',
              'Mapped CAM-SE wind speed + mountain-wave enhancement']
    fig, axs = plt.subplots(1, 2, figsize=(16, 7),
                            subplot_kw={'projection': proj}, constrained_layout=True)
    for ax, fld, title in zip(axs, fields, titles):
        cf = ax.contourf(glon, glat, np.ma.masked_invalid(fld), levels=wsplvs,
                         cmap=cmapN, extend='max', transform=data_crs)
        ax.contour(glon, glat, np.nan_to_num(topo_grid), levels=[1000, 2000],
                   colors='0.3', linewidths=0.4, transform=data_crs)
        ax.coastlines(resolution='10m', linewidth=0.7)
        ax.add_feature(cfeature.BORDERS, linewidth=0.6)
        ax.set_extent(extent, crs=data_crs)
        if carso:   # mark the Carso fire area
            ax.plot([13.5,13.9,13.9,13.5,13.5],[45.70,45.70,45.95,45.95,45.70],
                    color='k', lw=1.5, transform=data_crs)
        gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.6', linestyle='--')
        gl.top_labels = False; gl.right_labels = False
        ax.set_title(title, fontsize=13)
    cbar = fig.colorbar(cf, ax=axs, orientation='horizontal', shrink=0.5, pad=0.04)
    cbar.set_label(r'Wind speed (m s$^{-1}$)')
    fig.suptitle(f'Friuli / Carso bora event  {title_date}  ({caseN})', fontsize=15)
    return fig, axs

fig, axs = plot_friuli(cwsp_grid, uprm_grid, '2022-07-20 00Z')

# zoom on the Carso / Gulf of Trieste
fig, axs = plot_friuli(cwsp_grid, uprm_grid, '2022-07-20 00Z',
                       extent=(12.3, 14.8, 45.2, 46.6))

# %% ============================================================
# CELL 10 - NEW: loop over all timesteps.  Collect the grid fields,
#           then plot (a) max over the period, (b) a time series at Carso.
#           Painting is the slow step (~tens of s per step); 52 steps is fine
#           on a login node but go get coffee.
# ============================================================
dates = [f'2022-07-{d:02d}' for d in range(19, 32)]
secs  = [0, 21600, 43200, 64800]

t_lab, cw_t, up_t, upd_t, wm_t = [], [], [], [], []
for d in dates:
    for s in secs:
        fn = step_file(d, s)
        if not os.path.exists(fn):
            continue
        Xs = xr.open_dataset(fn)
        taudsw = Xs.TAUDSW.values; tauoro = Xs.TAUORO.values
        usrc = Xs.USRC_RDG.values; vsrc = Xs.VSRC_RDG.values
        nsrc = Xs.NSRC_RDG.values; ubsrc = Xs.UBMSRC_RDG.values
        hwdth = 1000.*Xs.HWDTH_RDG.values
        delprm     = np.sqrt(hwdth*taudsw/(rho*nsrc*ubsrc))
        delprm_oro = np.sqrt(hwdth*tauoro/(rho*nsrc*ubsrc))
        uprm = nsrc*delprm; uprm_dsw = nsrc*(delprm-delprm_oro)
        warmy = delprm*nsrc**2*Tsrc/9.8
        kw = dict(shape='rectangle', width_scale=1.0, verbose=False, ir_range=(0,43),
                  side='downwind', usrc=usrc[0,:,:], vsrc=vsrc[0,:,:])
        f_up  = np.nan_to_num(Pti.paint_tiles_nc(uprm[0,:,:],     tg, modulation='cosine',        **kw))
        f_upd = np.nan_to_num(Pti.paint_tiles_nc(uprm_dsw[0,:,:], tg, modulation='ramp_downwind', **kw))
        f_wm  = np.nan_to_num(Pti.paint_tiles_nc(warmy[0,:,:],    tg, modulation='ramp_downwind', **kw))
        cws = np.sqrt(Xs.U.values**2 + Xs.V.values**2)
        t_lab.append(f'{d} {s//3600:02d}Z')
        cw_t.append(mapper_se.se_to_wrf(cws[0,klev,:]))
        up_t.append(np.nan_to_num(mapper_cs.cs_to_wrf(f_up)))
        upd_t.append(np.nan_to_num(mapper_cs.cs_to_wrf(f_upd)))
        wm_t.append(np.nan_to_num(mapper_cs.cs_to_wrf(f_wm)))
        Xs.close()
        print(t_lab[-1], ' max uprm on grid %.1f m/s' % up_t[-1].max())

cw_t, up_t, upd_t, wm_t = map(np.array, (cw_t, up_t, upd_t, wm_t))   # (nt, ny, nx)

# %% ============================================================
# CELL 11 - NEW: period maximum and a Carso time series
# ============================================================
fig, axs = plot_friuli(cw_t.max(axis=0), up_t.max(axis=0), f'max over {t_lab[0]} - {t_lab[-1]}')

# time series averaged over the Carso box
cb = (glon>=13.5)&(glon<=13.9)&(glat>=45.70)&(glat<=45.95)
fig, ax = plt.subplots(figsize=(14,4))
ax.plot(cw_t[:,cb].mean(axis=1),                      label='CAM-SE wind')
ax.plot((cw_t+enh*up_t)[:,cb].mean(axis=1),           label=f'+ {enh} u\'')
ax.plot(upd_t[:,cb].mean(axis=1),                     label="u'_DSW only")
ax.set_xticks(range(len(t_lab))); ax.set_xticklabels(t_lab, rotation=90, fontsize=7)
ax.set_ylabel('m/s'); ax.set_title('Carso box mean'); ax.legend(); ax.grid(alpha=.3)

# which step is the strongest at Carso?
ibest = int(np.argmax((cw_t+enh*up_t)[:,cb].mean(axis=1)))
print('strongest step at Carso:', t_lab[ibest])
fig, axs = plot_friuli(cw_t[ibest], up_t[ibest], t_lab[ibest], extent=(12.3, 14.8, 45.2, 46.6))

# ================================================================
# PATCH to se_wrf_mapper.py  (backwards compatible; CONUS stays default)
# ================================================================
#   def _project(lon, lat, lon0=-97.9, lat0=39.1):      # unchanged
#
#   class SE_WRF_Mapper:
#       def __init__(self, se_lon, se_lat, wrf_lon, wrf_lat,
#                    degree_pad=3.0, chunk_size=200000,
#                    lon0=-97.9, lat0=39.1):              # <-- add these two
#           ...
#           sx, sy = _project(self.se_lon[self.src_index],
#                             self.se_lat[self.src_index],
#                             lon0=lon0, lat0=lat0)      # <-- pass them
#           ...
#           wx, wy = _project(self.wrf_lon, self.wrf_lat,
#                             lon0=lon0, lat0=lat0)      # <-- and here
"""
cam_hist.py -- open CAM history output by case, stream and time window.

    from cam_hist import open_cam_hist

    # Monthly-mean stream, DJF for 1990-1999 (Dec 1989 ... Feb 1999)
    ds = open_cam_hist("/glade/derecho/scratch/juliob/archive", "mycase", "h0a",
                       years=(1990, 1999), season="DJF")

    # Instantaneous stream, a few variables, JJA of one year
    ds = open_cam_hist(basedir, casename, "h1i", years=2005, season="JJA",
                       variables=["T", "KVH_CLUBB", "CLUBB_DTHLM_SFC"])

    # Just see which files would be read
    files = find_cam_hist_files(basedir, casename, "h0a", years=(1990, 1999))

Conventions assumed (standard CESM short-term archive):

    {basedir}/{casename}/atm/hist/{casename}.cam.{stream}.{date}.nc

where {date} is YYYY-MM for monthly files or YYYY-MM-DD-SSSSS otherwise, and
is the START of the period the file covers.  Override `hist_dir` for a run
directory or any other layout.

Seasons and years
-----------------
`season` is "DJF", "MAM", "JJA", "SON", "ANN", or any run of consecutive
month initials ("NDJFM", "JJAS").  Alternatively pass `months` as a list of
integers.  Seasons that cross the year boundary are labelled by the year of
their LAST month, so years=(1990, 1999), season="DJF" gives the ten complete
winters Dec 1989-Feb 1990 through Dec 1998-Feb 1999, and never an orphan
December or a winter missing its December.

Time stamps
-----------
Older CAM stamps a time-averaged sample at the END of its averaging interval,
so the January mean carries the date 1 February.  Selecting months on that
axis silently shifts everything by one month.  With `center_time=True` (the
default) the time coordinate is replaced by the midpoint of `time_bnds`
whenever the file has bounds, which is correct for both old and new CAM
conventions and is a no-op for instantaneous streams, which have no bounds.
"""

from __future__ import annotations

import difflib
import glob
import os
import re
import warnings
from typing import Iterable, Sequence

import numpy as np
import xarray as xr

__all__ = ["open_cam_hist", "find_cam_hist_files", "season_to_months"]

_MONTH_LETTERS = "JFMAMJJASOND"
_NAMED_SEASONS = {
    "DJF": [12, 1, 2],
    "MAM": [3, 4, 5],
    "JJA": [6, 7, 8],
    "SON": [9, 10, 11],
    "ANN": list(range(1, 13)),
}

# Small variables that are always worth carrying along when `variables` is
# given: time bounds, the hybrid-coordinate coefficients needed to rebuild
# pressure, and grid weights.
_ALWAYS_KEEP = (
    "time_bnds", "time_bounds", "date", "datesec",
    "hyam", "hybm", "hyai", "hybi", "P0", "gw", "area",
)

# Date in the file name: YYYY-MM[-DD[-SSSSS]] immediately before ".nc".
_DATE_RE = re.compile(r"\.(\d{4,})-(\d{2})(?:-(\d{2}))?(?:-(\d{5}))?\.nc$")


# ---------------------------------------------------------------------------
# Month / season helpers
# ---------------------------------------------------------------------------

def season_to_months(season: str) -> list[int]:
    """'DJF' -> [12, 1, 2]; 'NDJFM' -> [11, 12, 1, 2, 3]; 'ANN' -> 1..12."""
    s = season.strip().upper()
    if s in _NAMED_SEASONS:
        return list(_NAMED_SEASONS[s])
    if not s or len(s) > 12 or not set(s) <= set(_MONTH_LETTERS):
        raise ValueError(f"Unrecognised season {season!r}.")
    idx = (_MONTH_LETTERS * 2).find(s)
    if idx < 0:
        raise ValueError(
            f"Season {season!r} is not a run of consecutive month initials "
            f"(e.g. 'DJF', 'JJAS', 'NDJFM')."
        )
    return [(idx + i) % 12 + 1 for i in range(len(s))]


def _season_year_offset(months: Sequence[int]) -> dict[int, int]:
    """
    Map month -> offset to add to the calendar year to get the season year.

    Months before the year wrap (e.g. Nov, Dec in NDJFM) belong to the season
    that ends in the following calendar year, so they get +1.
    """
    wrap = next((i for i in range(1, len(months)) if months[i] < months[i - 1]), None)
    if wrap is None:
        return {m: 0 for m in months}
    return {m: (1 if i < wrap else 0) for i, m in enumerate(months)}


def _normalise_years(years) -> tuple[int, int] | None:
    if years is None:
        return None
    if isinstance(years, (int, np.integer)):
        return int(years), int(years)
    y0, y1 = years
    if y1 < y0:
        raise ValueError(f"years=({y0}, {y1}) is reversed.")
    return int(y0), int(y1)


def _resolve_months(season, months) -> list[int] | None:
    if season is not None and months is not None:
        raise ValueError("Give `season` or `months`, not both.")
    if season is not None:
        return season_to_months(season)
    if months is not None:
        months = [int(m) for m in ([months] if np.isscalar(months) else months)]
        bad = [m for m in months if not 1 <= m <= 12]
        if bad:
            raise ValueError(f"Invalid month(s): {bad}")
        return months
    return None


def _wanted(year: int, month: int, yr, months, offset) -> bool:
    """Is calendar (year, month) inside the requested season-years/months?"""
    if months is not None and month not in months:
        return False
    if yr is None:
        return True
    season_year = year + (offset.get(month, 0) if offset else 0)
    return yr[0] <= season_year <= yr[1]


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

def _hist_dir(basedir, casename, hist_dir):
    if hist_dir is not None:
        return os.fspath(hist_dir)
    return os.path.join(os.fspath(basedir), casename, "atm", "hist")


def _file_start(path: str) -> tuple[int, int, int, int] | None:
    m = _DATE_RE.search(os.path.basename(path))
    if m is None:
        return None
    y, mo, d, s = m.groups()
    return int(y), int(mo), int(d or 1), int(s or 0)


def _months_covered(start, nxt):
    """
    Calendar (year, month) pairs touched by a file covering [start, nxt).

    `nxt` is the next file's start.
    """
    y, m = start[0], start[1]
    ey, em = nxt[0], nxt[1]
    # If the next file starts exactly at the first instant of its month, that
    # month is not covered by this file.
    if nxt[2] == 1 and nxt[3] == 0:
        em -= 1
        if em == 0:
            ey, em = ey - 1, 12
    out = []
    while (y, m) <= (ey, em):
        out.append((y, m))
        m += 1
        if m == 13:
            y, m = y + 1, 1
    return out or [(start[0], start[1])]


def _add_months(y, m, k):
    t = y * 12 + (m - 1) + k
    return t // 12, t % 12 + 1


def _available_streams(hdir, casename, component):
    pat = os.path.join(hdir, f"{casename}.{component}.*.nc")
    rx = re.compile(rf"^{re.escape(casename)}\.{re.escape(component)}\.([^.]+)\.")
    streams = set()
    for f in glob.glob(pat):
        m = rx.match(os.path.basename(f))
        if m and _file_start(f) is not None:
            streams.add(m.group(1))
    return sorted(streams)


def find_cam_hist_files(
    basedir,
    casename: str,
    stream: str,
    years=None,
    season: str | None = None,
    months: Iterable[int] | int | None = None,
    component: str = "cam",
    hist_dir=None,
    pad_months: int = 0,
) -> list[str]:
    """
    Return the sorted history files for `stream` that can contain data in the
    requested years/season.  Files are chosen from the dates in their names;
    a file is kept if any month it covers is wanted.  Exact per-sample
    selection happens in open_cam_hist.  `pad_months` widens each file's
    assumed coverage on both sides, for time stamps that may sit outside the
    period named in the file (uncentred monthly means).
    """
    hdir = _hist_dir(basedir, casename, hist_dir)
    if not os.path.isdir(hdir):
        raise FileNotFoundError(f"History directory does not exist: {hdir}")

    pat = os.path.join(hdir, f"{casename}.{component}.{stream}.*.nc")
    dated = []
    for f in glob.glob(pat):
        st = _file_start(f)
        if st is not None:
            dated.append((st, f))

    if not dated:
        avail = _available_streams(hdir, casename, component)
        hint = f" Streams present: {', '.join(avail)}." if avail else ""
        close = difflib.get_close_matches(stream, avail, n=3)
        if close:
            hint += f" Did you mean {' or '.join(repr(c) for c in close)}?"
        raise FileNotFoundError(
            f"No files match {pat}.{hint or ' No CAM history files found there.'}"
        )

    dated.sort()
    yr = _normalise_years(years)
    mon = _resolve_months(season, months)
    offset = _season_year_offset(mon) if mon else None

    if yr is None and mon is None:
        return [f for _, f in dated]

    keep = []
    for i, (st, f) in enumerate(dated):
        if i + 1 < len(dated):
            covered = _months_covered(st, dated[i + 1][0])
        elif i > 0:
            # Last file: assume it spans as many months as the one before it.
            n = len(_months_covered(dated[i - 1][0], st))
            covered = [_add_months(st[0], st[1], k) for k in range(n)]
        else:
            covered = [(st[0], st[1])]
        if pad_months:
            y0, m0 = covered[0]
            y1, m1 = covered[-1]
            covered = ([_add_months(y0, m0, -k) for k in range(pad_months, 0, -1)]
                       + covered
                       + [_add_months(y1, m1, k) for k in range(1, pad_months + 1)])
        if any(_wanted(y, m, yr, mon, offset) for y, m in covered):
            keep.append(f)
    return keep


# ---------------------------------------------------------------------------
# Opening
# ---------------------------------------------------------------------------

def _cftime_kwargs() -> dict:
    """Ask xarray for cftime decoding in a way that suits old and new versions."""
    coders = getattr(xr, "coders", None)
    if coders is not None and hasattr(coders, "CFDatetimeCoder"):
        return {"decode_times": coders.CFDatetimeCoder(use_cftime=True)}
    return {"use_cftime": True}


def _make_preprocess(variables):
    if variables is None:
        return None
    wanted = list(dict.fromkeys(variables))

    def _pre(ds):
        missing = [v for v in wanted if v not in ds.variables]
        if missing:
            msgs = []
            for v in missing:
                close = difflib.get_close_matches(v, list(ds.data_vars), n=3, cutoff=0.6)
                msgs.append(f"{v!r}" + (f" (close: {', '.join(close)})" if close else ""))
            src = ds.encoding.get("source", "file")
            raise KeyError(f"Not in {os.path.basename(src)}: " + "; ".join(msgs))
        extra = [v for v in _ALWAYS_KEEP if v in ds.variables and v not in wanted]
        return ds[wanted + extra]

    return _pre


def _center_time(ds: xr.Dataset) -> xr.Dataset:
    bname = ds["time"].attrs.get("bounds") or ds["time"].encoding.get("bounds")
    if not bname or bname not in ds:
        bname = next((b for b in ("time_bnds", "time_bounds") if b in ds), None)
    if bname is None:
        return ds
    tb = ds[bname]
    if tb.dtype.kind in "iuf":
        warnings.warn(f"{bname} was not decoded to dates; time left uncentred.")
        return ds
    vals = tb.transpose("time", ...).values
    mid = np.array([a + (b - a) / 2 for a, b in zip(vals[:, 0], vals[:, 1])], dtype=object)
    attrs = dict(ds["time"].attrs)
    enc = dict(ds["time"].encoding)
    ds = ds.assign_coords(time=("time", mid))
    ds["time"].attrs.update(attrs)
    ds["time"].attrs["comment"] = f"Replaced by the midpoint of {bname} (cam_hist.open_cam_hist)."
    ds["time"].encoding.update(enc)
    return ds


def open_cam_hist(
    basedir,
    casename: str,
    stream: str,
    years=None,
    season: str | None = None,
    months: Iterable[int] | int | None = None,
    variables: Sequence[str] | str | None = None,
    component: str = "cam",
    hist_dir=None,
    center_time: bool = True,
    parallel: bool = False,
    chunks=None,
    **mfdataset_kwargs,
) -> xr.Dataset:
    """
    Open CAM history files for one case and stream, restricted to a time window.

    Parameters
    ----------
    basedir, casename
        Files are looked for in {basedir}/{casename}/atm/hist unless
        `hist_dir` is given, in which case `basedir` is ignored.
    stream
        History stream tag as it appears in the file names: 'h0', 'h0a',
        'h1i', 'h2a', ...
    years
        A single year, an inclusive (first, last) pair, or None for all.
        With a season that crosses the year boundary these are season years
        (the year of the season's last month).
    season, months
        One or the other: a season string ('DJF', 'NDJFM', ...) or a list of
        month numbers.  Either can be omitted for all months.
    variables
        Variable name(s) to keep.  Bounds, hybrid coefficients (hyam, hybm,
        hyai, hybi, P0) and grid weights (gw, area) come along when present.
        A misspelt name raises KeyError listing close matches.  None keeps all.
    component
        Model-component tag in the file name, 'cam' by default.
    center_time
        Replace time by the midpoint of time_bnds where bounds exist (see the
        module docstring).  Selection is done on the centred axis.
    parallel, chunks
        Passed to xarray.open_mfdataset.  parallel=True opens files through
        dask, which helps on a compute node with many files.
    **mfdataset_kwargs
        Anything else for xarray.open_mfdataset; overrides the defaults below.

    Returns
    -------
    xarray.Dataset
        Lazily loaded (dask-backed), times decoded to cftime.  attrs record
        the case, stream, selection and number of files read.
    """
    if isinstance(variables, str):
        variables = [variables]

    files = find_cam_hist_files(
        basedir, casename, stream, years=years, season=season, months=months,
        component=component, hist_dir=hist_dir,
        pad_months=0 if center_time else 1,
    )
    if not files:
        raise FileNotFoundError(
            f"Stream {stream!r} of {casename!r} has files, but none in "
            f"years={years}, season={season}, months={months}."
        )

    kw = dict(
        combine="nested",
        concat_dim="time",
        data_vars="minimal",
        coords="minimal",
        compat="override",
        join="override",
        parallel=parallel,
        chunks=chunks,
        preprocess=_make_preprocess(variables),
    )
    kw.update(_cftime_kwargs())
    kw.update(mfdataset_kwargs)

    ds = xr.open_mfdataset(files, **kw)

    if center_time and "time" in ds.dims:
        ds = _center_time(ds)

    # Exact per-sample selection.
    yr = _normalise_years(years)
    mon = _resolve_months(season, months)
    if (yr is not None or mon is not None) and "time" in ds.dims:
        t_year = ds["time"].dt.year.values
        t_mon = ds["time"].dt.month.values
        offset = _season_year_offset(mon) if mon else None
        mask = np.array(
            [_wanted(int(y), int(m), yr, mon, offset) for y, m in zip(t_year, t_mon)],
            dtype=bool,
        )
        if not mask.any():
            raise ValueError(
                f"Files were found but no time samples fall in years={years}, "
                f"season={season}, months={months}."
            )
        ds = ds.isel(time=np.flatnonzero(mask))

    ds.attrs.update(
        cam_hist_case=casename,
        cam_hist_stream=stream,
        cam_hist_selection=f"years={years}, season={season}, months={months}",
        cam_hist_nfiles=len(files),
        cam_hist_time_centered=int(bool(center_time)),
    )
    return ds
"""
rdgem.py  --  surface METARs from IEM mtarchive GEMPAK daily files

    from rdgem import fetch_gem, read_gem_metars
    from rdsynop import window6h

    files = fetch_gem(['2022-07-16', '2022-07-31'], outdir='gem')   # start, end (inclusive)
    S = read_gem_metars(files, station='LIPQ')     # same columns as read_metar / read_synop
    W = window6h(S, bora_sector=(20, 130))

The files live at
    https://mtarchive.geol.iastate.edu/YYYY/MM/DD/gempak/surface/sao/YYYYMMDD_sao.gem
(one UTC day each, ~20-25 MB, ~3000 stations worldwide).  They are binary
GEMPAK, but each report's raw METAR text is embedded as plain ASCII, so we
simply scan for "<STATION> DDHHMMZ ..." and decode the METAR.

Returned columns (same as the other readers):
    time (UTC)  station  dd (deg, NaN=calm/variable)  ff (m/s)  gust (m/s, NaN if none)
    T, Td (deg C)  p_msl (hPa; from Q group or altimeter A group)
"""
import os
import re
import glob
import datetime as dt
import urllib.request
import numpy as np
import pandas as pd

__all__ = ["fetch_gem", "read_gem_metars", "decode_metar"]

KT = 0.514444
MPS = 1.0

_URL = "https://mtarchive.geol.iastate.edu/{y}/{m:02d}/{d:02d}/gempak/surface/sao/{y}{m:02d}{d:02d}_sao.gem"


def fetch_gem(dates, outdir='gem', overwrite=False, verbose=True):
    """
    Download daily files.  dates = ['YYYY-MM-DD', 'YYYY-MM-DD'] -> every day
    from the first to the second inclusive; or a list of individual days.
    Returns the list of local file paths (existing files are not re-downloaded).
    """
    os.makedirs(outdir, exist_ok=True)
    if len(dates) == 2:
        d0 = dt.date.fromisoformat(dates[0]); d1 = dt.date.fromisoformat(dates[1])
        days = [d0 + dt.timedelta(n) for n in range((d1 - d0).days + 1)]
    else:
        days = [dt.date.fromisoformat(s) for s in dates]
    out = []
    for d in days:
        url = _URL.format(y=d.year, m=d.month, d=d.day)
        fn = os.path.join(outdir, os.path.basename(url))
        if os.path.exists(fn) and not overwrite:
            out.append(fn); continue
        if verbose: print('fetching', url)
        try:
            urllib.request.urlretrieve(url, fn)
            out.append(fn)
        except Exception as e:
            print('  FAILED:', e)
    return out


# ---- METAR decoding ---------------------------------------------------------
_WIND = re.compile(r'\b(\d{3}|VRB)(\d{2,3})(?:G(\d{2,3}))?(KT|MPS)\b')
_TEMP = re.compile(r'\b(M?\d{2})/(M?\d{2})\b')
_QNH  = re.compile(r'\bQ(\d{4})\b')
_ALT  = re.compile(r'\bA(\d{4})\b')


def _t(s):
    return -float(s[1:]) if s.startswith('M') else float(s)


def decode_metar(body):
    """body = everything after 'STATION DDHHMMZ'.  Returns dict."""
    dd = ff = gust = np.nan
    m = _WIND.search(body)
    if m:
        unit = KT if m.group(4) == 'KT' else MPS
        ff = int(m.group(2)) * unit
        if m.group(3):
            gust = int(m.group(3)) * unit
        if m.group(1) != 'VRB' and not (m.group(1) == '000' and ff == 0):
            dd = float(m.group(1))
    T = Td = np.nan
    m = _TEMP.search(body)
    if m:
        T, Td = _t(m.group(1)), _t(m.group(2))
    p = np.nan
    m = _QNH.search(body)
    if m:
        p = float(m.group(1))
    else:
        m = _ALT.search(body)
        if m:
            p = int(m.group(1)) / 100.0 * 33.8639
    return dict(dd=dd, ff=ff, gust=gust, T=T, Td=Td, p_msl=p)


def read_gem_metars(files, station):
    """
    files: list of paths or a glob pattern.  station: ICAO id, e.g. 'LIPQ'.
    Scans the binary files for raw METAR text of that station.
    """
    if isinstance(files, str):
        files = sorted(glob.glob(files))
    st = station.encode()
    pat = re.compile(rb'(?<![A-Z0-9])' + st + rb' (\d{2})(\d{2})(\d{2})Z ([ -~]{3,300}?)(?=[^ -~])')
    rows = []
    seen = set()
    for fn in files:
        m = re.search(r'(\d{4})(\d{2})(\d{2})_sao', os.path.basename(fn))
        if not m:
            raise ValueError(f'cannot get the date from file name {fn}')
        fy, fm, fd = map(int, m.groups())
        fdate = dt.date(fy, fm, fd)
        data = open(fn, 'rb').read()
        for h in pat.finditer(data):
            day, hh, mm = int(h.group(1)), int(h.group(2)), int(h.group(3))
            # The METAR carries only the day-of-month.  Normally the file's
            # year/month apply; a day-1 file also holds the last reports of
            # the previous month (e.g. '312350Z' in 20220801_sao.gem).
            try:
                d = dt.date(fy, fm, day)
                if d > fdate + dt.timedelta(days=1):
                    raise ValueError
            except ValueError:
                prev = fdate.replace(day=1) - dt.timedelta(days=1)
                d = dt.date(prev.year, prev.month, day)
            t = dt.datetime(d.year, d.month, d.day, hh, mm)
            body = h.group(4).decode(errors='replace').strip()
            key = (t, body[:20])
            if key in seen:
                continue
            seen.add(key)
            r = decode_metar(body)
            r.update(time=pd.Timestamp(t), station=station, raw=body)
            rows.append(r)
    if not rows:
        raise ValueError(f'no METARs for {station} in {len(files)} file(s)')
    S = pd.DataFrame(rows)[['time', 'station', 'dd', 'ff', 'gust', 'T', 'Td', 'p_msl', 'raw']]
    return S.sort_values('time').drop_duplicates('time').reset_index(drop=True)


if __name__ == '__main__':
    import sys
    S = read_gem_metars(sys.argv[1:-1], sys.argv[-1])
    pd.set_option('display.width', 200)
    print(S.drop(columns='raw').head(10).to_string())
    print(len(S), 'reports', S.time.min(), '->', S.time.max())
"""
rdsynop.py  --  read OGIMET SYNOP text dumps (FM-12 AAXX land reports)

    from rdsynop import read_synop, window6h
    S = read_synop('14106_jul2022.txt')          # pandas DataFrame, one row per report
    W = window6h(S, bora_sector=(20, 130))        # 6-hourly summary at 00/06/12/18Z

Decodes only what the DSW work needs: wind direction/speed, gust, T, Td,
station pressure, MSL pressure.  Everything else in the message is ignored.

Columns of S:
    time   (UTC, pandas Timestamp)      dd   wind direction, degrees (NaN if variable/missing)
    ff     wind speed, m/s              gust max gust in last hour (911ff) or 10 min (910ff), m/s
    T, Td  deg C                         p_stn, p_msl  hPa
    station  the IIiii index or, if absent, the WIGOS id

Units: the iw digit in the YYGGiw group says whether ff is m/s (0/1) or
knots (3/4); knots are converted so ff is always m/s.
"""
import re
import numpy as np
import pandas as pd

__all__ = ["read_synop", "window6h"]

# one report = "YYYYMMDDHHMM  WIGOS  IIiii  AAXX ... ==" possibly over several lines
_REPORT = re.compile(r'^(\d{12})\s+(\S+)\s+(\S+)\s+(AAXX .*?)==', re.M | re.S)


def _temp(group):
    """1sTTT / 2sTTT -> deg C; sign digit 0 = +, 1 = -"""
    if not re.fullmatch(r'[12]\d{4}', group):
        return np.nan
    s = -1.0 if group[1] == '1' else 1.0
    return s * int(group[2:]) / 10.0


def _pres(group):
    """3PPPP / 4PPPP -> hPa (PPPP in tenths, leading 1 dropped if > 1000)"""
    if not re.fullmatch(r'[34]\d{4}', group):
        return np.nan
    p = int(group[1:]) / 10.0
    return p + 1000.0 if p < 500.0 else p


def read_synop(path):
    txt = open(path, encoding='utf-8', errors='replace').read()
    rows = []
    for tstr, wigos, iiiii, body in _REPORT.findall(txt):
        g = [x.strip('=') for x in body.split()]
        if len(g) < 5 or g[0] != 'AAXX':
            continue
        iw = g[1][4]                      # YYGGiw
        knots = iw in '34'
        # section 1: g[3]=iRixhVV, g[4]=Nddff, then 1sTTT 2sTTT 3PPPP 4PPPP ...
        nddff = g[4]
        dd = ff = np.nan
        if re.fullmatch(r'\d{5}', nddff):
            d = int(nddff[1:3]); f = int(nddff[3:5])
            if nddff[1:3] == '00' and f == 0:
                dd = np.nan                 # calm
            elif d == 99:
                dd = np.nan                 # variable
            else:
                dd = d * 10.0
            ff = f * 0.514444 if knots else float(f)
        elif re.fullmatch(r'\d\d{2}//', nddff):
            pass
        T = Td = p_stn = p_msl = np.nan
        gust = np.nan
        sec = 1
        for x in g[5:]:
            if x == '333':
                sec = 3; continue
            if x in ('222', '444', '555'):
                sec = int(x[0]); continue
            if sec == 1:
                if x.startswith('1') and np.isnan(T):   T = _temp(x)
                elif x.startswith('2') and np.isnan(Td): Td = _temp(x)
                elif x.startswith('3') and np.isnan(p_stn): p_stn = _pres(x)
                elif x.startswith('4') and np.isnan(p_msl): p_msl = _pres(x)
            elif sec == 3:
                if re.fullmatch(r'91[01]\d\d', x):
                    gv = int(x[3:]) * (0.514444 if knots else 1.0)
                    gust = gv if np.isnan(gust) else max(gust, gv)
        station = iiiii if re.fullmatch(r'\d{5}', iiiii) else wigos
        rows.append(dict(time=pd.Timestamp(tstr), station=station, dd=dd, ff=ff,
                         gust=gust, T=T, Td=Td, p_stn=p_stn, p_msl=p_msl))
    S = pd.DataFrame(rows).sort_values('time').reset_index(drop=True)
    return S


def window6h(S, bora_sector=(20, 130), half_width='3h'):
    """
    Summarise reports in windows centred on 00/06/12/18Z (same instants as the
    GW_UnitTest steps).  bora_sector = (lo, hi) degrees, inclusive.
    Returns a DataFrame indexed by window centre with:
        n, ff_mean, ff_max, gust_max, bora_frac, ff_bora_mean, ff_bora_max,
        T_mean, p_msl_mean
    """
    hw = pd.Timedelta(half_width)
    t0 = S.time.min().floor('6h'); t1 = S.time.max().ceil('6h')
    centres = pd.date_range(t0, t1, freq='6h')
    lo, hi = bora_sector
    out = []
    for c in centres:
        w = S[(S.time >= c - hw) & (S.time < c + hw)]
        if len(w) == 0:
            continue
        b = w[(w.dd >= lo) & (w.dd <= hi) & (w.ff > 0)]
        out.append(dict(time=c, n=len(w),
                        ff_mean=w.ff.mean(), ff_max=w.ff.max(), gust_max=w.gust.max(),
                        bora_frac=len(b) / len(w),
                        ff_bora_mean=b.ff.mean() if len(b) else 0.0,
                        ff_bora_max=b.ff.max() if len(b) else 0.0,
                        T_mean=w['T'].mean(), p_msl_mean=w.p_msl.mean()))
    return pd.DataFrame(out).set_index('time')


if __name__ == '__main__':
    import sys
    S = read_synop(sys.argv[1])
    print(S.head()); print(len(S), 'reports from', S.station.iloc[0])
    pd.set_option('display.width', 160)
    print(window6h(S).round(1).to_string())
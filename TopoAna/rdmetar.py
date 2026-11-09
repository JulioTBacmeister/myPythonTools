"""
rdmetar.py  --  read Iowa Environmental Mesonet (IEM) ASOS/METAR CSV downloads

    from rdmetar import read_metar
    from rdsynop import window6h
    S = read_metar('LIPQ_jul2022.csv')          # DataFrame, one row per METAR
    W = window6h(S, bora_sector=(20, 130))      # same 6-h summary as for SYNOPs

Returns the SAME columns as rdsynop.read_synop so the two can be used
interchangeably:
    time   (UTC)                dd    wind direction, deg (NaN = calm/variable/missing)
    ff     wind speed, m/s      gust  gust, m/s (NaN when METAR reported none)
    T, Td  deg C (NaN if not downloaded)        p_msl  hPa (NaN if not downloaded)
    station  ICAO id, e.g. LIPQ

Handles whichever columns were ticked on the IEM download form:
    wind speed : 'sped' (mph) or 'sknt' (knots)
    gust       : 'gust_mph' (mph) or 'gust' (knots)
    temperature: 'tmpf' (F) or 'tmpc' (C);  dewpoint 'dwpf'/'dwpc'
    pressure   : 'mslp' (hPa) or 'alti' (inches Hg, converted approximately)
'M' is IEM's missing-value marker.
"""
import numpy as np
import pandas as pd

__all__ = ["read_metar"]

MPH = 0.44704
KT  = 0.514444


def _col(df, name, scale=1.0, offset=0.0):
    if name not in df.columns:
        return pd.Series(np.nan, index=df.index)
    return pd.to_numeric(df[name], errors='coerce') * scale + offset


def read_metar(path):
    df = pd.read_csv(path, na_values=['M'], comment='#')
    S = pd.DataFrame()
    S['time']    = pd.to_datetime(df['valid'])
    S['station'] = df['station']
    S['dd']      = _col(df, 'drct')
    # speed
    if 'sped' in df.columns:
        S['ff'] = _col(df, 'sped', MPH)
    elif 'sknt' in df.columns:
        S['ff'] = _col(df, 'sknt', KT)
    else:
        raise ValueError("no wind-speed column ('sped' or 'sknt') in file")
    # gust
    if 'gust_mph' in df.columns:
        S['gust'] = _col(df, 'gust_mph', MPH)
    elif 'gust' in df.columns:
        S['gust'] = _col(df, 'gust', KT)
    else:
        S['gust'] = np.nan
    # temperature / dewpoint
    S['T']  = _col(df, 'tmpc') if 'tmpc' in df.columns else (_col(df, 'tmpf') - 32.0) * 5.0 / 9.0
    S['Td'] = _col(df, 'dwpc') if 'dwpc' in df.columns else (_col(df, 'dwpf') - 32.0) * 5.0 / 9.0
    # pressure
    if 'mslp' in df.columns:
        S['p_msl'] = _col(df, 'mslp')
    elif 'alti' in df.columns:
        S['p_msl'] = _col(df, 'alti', 33.8639)      # altimeter setting, inHg -> hPa
    else:
        S['p_msl'] = np.nan
    # calm: direction 0 with zero speed -> NaN direction
    S.loc[(S.dd == 0) & (S.ff == 0), 'dd'] = np.nan
    return S.sort_values('time').reset_index(drop=True)


if __name__ == '__main__':
    import sys
    from rdsynop import window6h
    S = read_metar(sys.argv[1])
    print(S.head()); print(len(S), 'reports from', S.station.iloc[0],
                           '; gusts reported in', int(S.gust.notna().sum()))
    pd.set_option('display.width', 160)
    print(window6h(S).round(1).to_string())
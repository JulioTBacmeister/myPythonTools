#!/usr/bin/env python3
"""
regrid_to_cube.py
=================

Put a model field defined on the target grid (ne30pg3 columns, ncol) onto the
high-resolution cubed sphere (ncube x ncube x 6), so that model output can be
viewed alongside the topography and the painted ridge footprints.

This needs no interpolation machinery. `remapridge2target` already computed
`itrgtC`: the index of the target column that OWNS each cube cell. Regridding
is then one fancy-index operation,

    cube_field = A[itrgt - 1]

and the result is piecewise constant - every cube cell takes its column's
value. That is what you want for comparing against ridge footprints (it shows
exactly which cube cells a column covers) but it is blocky, and it is NOT
conservative or smooth. See "Limitations" below.

Where the map comes from
------------------------
`ridge_ana.F90`, subroutine `remapridge2target`, under `ldevelopment_diags`::

    OPEN (unit = 911, file='./output/grid_remap_nc<ncube>_<grid>.dat', form="UNFORMATTED")
    write(911) ncube, npeaks
    write(911) dble(itrgtC)      <- the map, 1-based, 0 where unassigned
    write(911) dble(itrgxC)      <- same but -1 where mxdisC <= 0.1
    write(911) xs, ys, xspk, yspk, peaks%i, peaks%j

So the run must have been made with -z / --development_diags. create_topo.py
passes that by default.

The file is large (~900 MB at ncube=3000, since both maps are stored as
float64). `read_cube_target_map` caches the result as a compact int32 .npy,
216 MB, which loads in a second thereafter.

Usage
-----
::

    from regrid_to_cube import read_cube_target_map, regrid_to_cube

    itrgt = read_cube_target_map('output/grid_remap_nc3000_ne30pg3.dat',
                                 cache=True, verbose=True)

    cube = regrid_to_cube(phis, itrgt)          # phis is (ncol,)
    # cube is (6, ncube, ncube), same layout as rdremap and paint_tiles_nc

Limitations
-----------
* Piecewise constant. Every cube cell in a column gets the identical value, so
  the result has visible column-sized blocks. That is faithful to what the
  model actually resolves, but do not mistake it for downscaling.
* Not conservative in the strict sense, though it is close: cube cells are
  assigned to exactly one column and cell areas vary only slightly across a
  panel, so an area-weighted mean of the cube field returns nearly the original
  column value.
* Cells with itrgt == 0 were never assigned to any column (outside the exchange
  grid). They get `fill`.
"""

from __future__ import annotations

import os
import struct

import numpy as np

__all__ = ["read_cube_target_map", "regrid_to_cube", "cube_cell_counts"]


class _F77:
    """Minimal Fortran sequential unformatted reader with byte-order sniffing."""

    def __init__(self, path, endian="auto"):
        self.path = path
        self._fh = open(path, "rb")
        self.endian = self._sniff() if endian == "auto" else endian

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self._fh.close()

    def _sniff(self):
        head = self._fh.read(4)
        for e in ("<", ">"):
            (n,) = struct.unpack(e + "i", head)
            if 0 < n < 1 << 34:
                self._fh.seek(4 + n)
                t = self._fh.read(4)
                if len(t) == 4 and struct.unpack(e + "i", t)[0] == n:
                    self._fh.seek(0)
                    return e
        self._fh.seek(0)
        raise IOError(f"{self.path}: not a Fortran sequential unformatted file")

    def _marker(self):
        raw = self._fh.read(4)
        if len(raw) < 4:
            raise EOFError("end of file at a record marker")
        return struct.unpack(self.endian + "i", raw)[0]

    def peek_size(self):
        here = self._fh.tell()
        tot = 0
        try:
            while True:
                h = self._marker()
                tot += abs(h)
                self._fh.seek(abs(h), 1)
                self._marker()
                if h >= 0:
                    return tot
        finally:
            self._fh.seek(here)

    def skip(self):
        while True:
            h = self._marker()
            self._fh.seek(abs(h), 1)
            self._marker()
            if h >= 0:
                return

    def read(self):
        chunks = []
        while True:
            h = self._marker()
            b = self._fh.read(abs(h))
            if len(b) != abs(h):
                raise EOFError("end of file inside a record")
            chunks.append(b)
            if self._marker() != h:
                raise IOError("record markers disagree")
            if h >= 0:
                return b"".join(chunks)


def read_cube_target_map(path, cache=True, cache_path=None, verbose=False):
    """Read `itrgtC` from a grid_remap_*.dat file.

    Returns (6, ncube, ncube) int32, 1-based column index, 0 where unassigned.

    Parameters
    ----------
    cache : bool
        Write / reuse a compact .npy alongside the .dat. The .dat stores the
        map as float64 and is ~900 MB at ncube=3000; the cache is int32.
    """
    if cache_path is None:
        cache_path = os.path.splitext(path)[0] + "_itrgt.npy"

    if cache and os.path.exists(cache_path) and \
            os.path.getmtime(cache_path) >= os.path.getmtime(path):
        a = np.load(cache_path)
        if verbose:
            print(f"loaded cached map {cache_path}  {a.shape}")
        return a

    with _F77(path) as f:
        i4 = np.dtype("i4").newbyteorder(f.endian)
        f8 = np.dtype("f8").newbyteorder(f.endian)
        hdr = np.frombuffer(f.read(), dtype=i4, count=2)
        ncube, npeaks = int(hdr[0]), int(hdr[1])
        ncells = 6 * ncube * ncube

        got = f.peek_size()
        if got != ncells * 8:
            raise IOError(f"{path}: itrgtC record is {got} bytes, expected "
                          f"{ncells*8} (6*{ncube}^2*8). Wrong file, or the write "
                          "order in remapridge2target has changed.")
        itrgt = np.frombuffer(f.read(), dtype=f8, count=ncells)
        itrgt = np.rint(itrgt).astype(np.int32).reshape(6, ncube, ncube)

    if cache:
        np.save(cache_path, itrgt)
        if verbose:
            print(f"cached map to {cache_path}")

    if verbose:
        nz = itrgt > 0
        print(f"read {path}")
        print(f"  ncube={ncube}  npeaks={npeaks:,}  cells={itrgt.size:,}")
        print(f"  assigned to a column: {100*nz.mean():.2f}%")
        print(f"  column index range: {itrgt[nz].min()} to {itrgt[nz].max()}")
    return itrgt


def cube_cell_counts(itrgt, ncol=None):
    """Cube cells owned by each column. Useful as a sanity check: at ncube=3000
    and ne30pg3 this should be ~1100 per column."""
    n = int(itrgt.max()) if ncol is None else ncol
    return np.bincount(itrgt.ravel(), minlength=n + 1)[1:]


def regrid_to_cube(A, itrgt, fill=np.nan, dtype=np.float32, verbose=False):
    """Scatter a per-column field onto the cubed sphere.

    Parameters
    ----------
    A : (ncol,) array
        Field on the target grid, in the same column order the topo run used.
    itrgt : (6, ncube, ncube) int32
        From read_cube_target_map.
    fill : float
        Value for cube cells with itrgt == 0.

    Returns
    -------
    (6, ncube, ncube) array
    """
    A = np.asarray(A, dtype=float).ravel()
    hi = int(itrgt.max())
    if hi > A.size:
        raise ValueError(
            f"map references column {hi} but A has only {A.size}. "
            "The grid_remap file and this field are from different target grids.")
    if hi < A.size and verbose:
        print(f"  note: map only references columns 1..{hi} of {A.size}")

    # index 0 holds the fill, so A can be indexed directly by itrgt
    lut = np.empty(A.size + 1, dtype=dtype)
    lut[0] = fill
    lut[1:] = A
    out = lut[itrgt]

    if verbose:
        bad = itrgt == 0
        print(f"  regridded to {out.shape}; {100*bad.mean():.2f}% unassigned -> fill")
    return out


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Inspect a grid_remap_*.dat cube->column map")
    p.add_argument("path")
    p.add_argument("--no-cache", action="store_true")
    a = p.parse_args()
    it = read_cube_target_map(a.path, cache=not a.no_cache, verbose=True)
    c = cube_cell_counts(it)
    print(f"  cube cells per column: min={c.min()} median={np.median(c):.0f} "
          f"max={c.max()} mean={c.mean():.1f}")
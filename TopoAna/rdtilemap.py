#!/usr/bin/env python3
"""
rdtilemap.py
============

Python translation of the IDL procedure ``rdtilemap`` that reads
``output/Ridge_tile_map.dat``, the file written by the Fortran subroutine
``remapridge2tiles``.

The Fortran writer (for reference)
----------------------------------
::

    OPEN (unit = 911, file = trim(ofile), form="UNFORMATTED")
    write(911) ntarget,nalloc,maxtiles,npeaks     ! record  1
    write(911) idxmap,idcoun                      ! record  2
    write(911) NumObjects                         ! record  3
    write(911) MyObject(:,1:maxtiles)             ! record  4
    write(911) WtObject(:,1:maxtiles)             ! record  5
    write(911) LnObject(:,1:maxtiles)             ! record  6
    write(911) xwoid_tiles                        ! record  7
    write(911) ywoid_tiles                        ! record  8
    write(911) xloid_tiles                        ! record  9
    write(911) yloid_tiles                        ! record 10
    write(911) panel_tiles                        ! record 11
    write(911) lonc_tiles                         ! record 12
    write(911) latc_tiles                         ! record 13
    write(911) lonw_tiles                         ! record 14
    write(911) latw_tiles                         ! record 15
    write(911) aniso_tiles                        ! record 16
    write(911) anglx_tiles                        ! record 17
    write(911) mxdis_tiles                        ! record 18
    write(911) hwdth_tiles                        ! record 19
    write(911) clngt_tiles                        ! record 20
    write(911) clext_tiles                        ! record 21
    write(911) clngt2                             ! record 22

``form="UNFORMATTED"`` with default (sequential) access means every ``write``
produces one logical record bracketed by 4-byte record-length markers, which is
exactly what IDL's ``openr,/f77`` handles and what ``F77SequentialReader``
below reproduces.  gfortran splits records larger than 2 GiB into subrecords
flagged by a negative leading marker; that case is handled too.

Shapes and index order
----------------------
Fortran declarations::

    integer  :: idxmap(nalloc, ntarget), idcoun(ntarget)
    integer  :: NumObjects(ntarget), MyObject(ntarget, nalloc)
    real(r8) :: WtObject(ntarget, nalloc), LnObject(ntarget, nalloc)
    real(r8) :: {x,y}{w,l}oid_tiles(ntarget, maxtiles)
    integer  :: panel_tiles(ntarget, maxtiles)
    real(r8) :: {lon,lat}{c,w}_tiles(ntarget, maxtiles)
    real(r8) :: aniso/anglx/mxdis/hwdth/clngt/clext_tiles(ntarget, maxtiles)

Note that ``MyObject``/``WtObject``/``LnObject`` are allocated with second
dimension ``nalloc`` but only the first ``maxtiles`` columns are written, so on
disk all the tile arrays are ``(ntarget, maxtiles)``.

IDL's first array index varies fastest, same as Fortran, so
``lonarr(nalloc,ntarget)`` in the IDL code is the same memory layout as the
Fortran ``idxmap(nalloc,ntarget)``.  NumPy defaults to C order, so every array
here is read flat and reshaped with ``order='F'``; the returned arrays have the
*same* subscript order as the Fortran/IDL ones (``idxmap[i, j]`` is cube cell
``i`` of target cell ``j``, ``mxdis[j, ir]`` is tile ``ir`` of target ``j``).

Known quirks carried over from the IDL
--------------------------------------
* The IDL declares ``clngt2 = fltarr(npeaks)`` (32-bit) while the Fortran
  module variable is very likely ``real(r8)``.  Rather than guess, this reader
  infers the word size of the last record from its length on disk.
* The IDL has a ``stop`` before ``close,1`` and before the ``r`` structure is
  built, so as written it drops into the debugger and never returns anything
  unless you type ``.cont``.  Presumably a leftover debugging line.
* The IDL structure ``r`` omits ``xwoid``, ``ywoid``, ``xloid``, ``yloid``,
  ``panel``, ``lonc``, ``latc``, ``lonw`` and ``latw``.  They are read but
  discarded.  This translation keeps them all.
* Missing tiles are flagged ``-9999.0`` (reals) / ``-1`` (panel) by the
  Fortran; see :meth:`TileMap.masked` for a convenience view.
* ``idxmap`` holds 1-based Fortran indices into the ``ncube*ncube*6`` cube
  arrays; subtract 1 before using them to subscript a NumPy array.

Usage
-----
::

    from rdtilemap import rdtilemap
    r = rdtilemap("output/Ridge_tile_map.dat")
    print(r.ntarget, r.maxtiles, r.mxdis.shape)
"""

from __future__ import annotations

import struct
from dataclasses import dataclass
from typing import BinaryIO

import numpy as np

__all__ = ["TileMap", "rdtilemap", "F77SequentialReader"]


# ---------------------------------------------------------------------------
# Minimal reader for Fortran sequential unformatted files (IDL's openr,/f77)
# ---------------------------------------------------------------------------
class F77SequentialReader:
    """Read records from a Fortran ``form="UNFORMATTED"``, sequential file.

    Each record is ``[int32 len][payload][int32 len]``.  A negative leading
    marker means the record was split into subrecords (gfortran's >2 GiB
    convention) and the pieces are concatenated transparently.
    """

    def __init__(self, path: str, endian: str = "auto"):
        self.path = path
        self._fh: BinaryIO = open(path, "rb")
        self.endian = self._detect_endian() if endian == "auto" else endian

    # -- context manager -------------------------------------------------
    def __enter__(self) -> "F77SequentialReader":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def close(self) -> None:
        if not self._fh.closed:
            self._fh.close()

    # -- internals -------------------------------------------------------
    def _detect_endian(self) -> str:
        """Sniff byte order from the first record marker.

        The first record is four default integers, i.e. 16 bytes, but we stay
        generic: accept whichever byte order gives a plausible, self-consistent
        record length.
        """
        head = self._fh.read(4)
        if len(head) < 4:
            raise IOError(f"{self.path}: file too short to be a Fortran unformatted file")
        for endian in ("<", ">"):
            (n,) = struct.unpack(endian + "i", head)
            if n <= 0 or n > 1 << 34:
                continue
            self._fh.seek(4 + n)
            tail = self._fh.read(4)
            if len(tail) == 4 and struct.unpack(endian + "i", tail)[0] == n:
                self._fh.seek(0)
                return endian
        self._fh.seek(0)
        raise IOError(
            f"{self.path}: does not look like a Fortran sequential unformatted "
            "file (record markers not self-consistent in either byte order)"
        )

    def _read_marker(self) -> int:
        raw = self._fh.read(4)
        if len(raw) < 4:
            raise EOFError("unexpected end of file while reading a record marker")
        return struct.unpack(self.endian + "i", raw)[0]

    # -- public API ------------------------------------------------------
    def read_record(self) -> bytes:
        """Return the payload of the next logical record."""
        chunks = []
        while True:
            head = self._read_marker()
            nbytes = abs(head)
            payload = self._fh.read(nbytes)
            if len(payload) != nbytes:
                raise EOFError("unexpected end of file inside a record payload")
            chunks.append(payload)
            tail = self._read_marker()
            if tail != head:
                raise IOError(
                    f"record markers disagree ({head} vs {tail}); "
                    "wrong byte order or corrupt file"
                )
            if head >= 0:  # last (or only) subrecord
                break
        return b"".join(chunks)

    def read_arrays(self, *specs) -> list[np.ndarray]:
        """Split one record into several arrays.

        Each spec is ``(dtype, shape)``.  Shapes are interpreted in Fortran
        order, matching the Fortran declaration.
        """
        buf = self.read_record()
        out, offset = [], 0
        for dtype, shape in specs:
            dt = np.dtype(dtype).newbyteorder(self.endian)
            count = int(np.prod(shape)) if shape else 1
            nbytes = count * dt.itemsize
            arr = np.frombuffer(buf, dtype=dt, count=count, offset=offset)
            out.append(np.reshape(arr, shape, order="F").astype(dt.newbyteorder("="), copy=True))
            offset += nbytes
        if offset != len(buf):
            raise IOError(
                f"record length mismatch: expected {offset} bytes, record holds {len(buf)}. "
                "Check integer/real kinds or the record ordering."
            )
        return out

    def peek_record_size(self) -> int:
        """Length in bytes of the next record, without consuming it."""
        here = self._fh.tell()
        total = 0
        try:
            while True:
                head = self._read_marker()
                total += abs(head)
                self._fh.seek(abs(head), 1)
                self._read_marker()
                if head >= 0:
                    return total
        finally:
            self._fh.seek(here)


# ---------------------------------------------------------------------------
# The container returned in place of the IDL anonymous structure `r`
# ---------------------------------------------------------------------------
@dataclass
class TileMap:
    """Everything in Ridge_tile_map.dat.

    Superset of the IDL structure ``r``: the IDL threw away the centroid and
    lon/lat arrays after reading them, this keeps them.
    """

    file: str

    # record 1 -- header
    ntarget: int
    nalloc: int
    maxtiles: int
    npeaks: int

    # record 2 -- cube-cell -> target-cell map
    idxmap: np.ndarray  # (nalloc, ntarget) int32, 1-based cube indices
    idcoun: np.ndarray  # (ntarget,)        int32, cube cells per target cell

    # records 3-6 -- object bookkeeping
    NumObjects: np.ndarray  # (ntarget,)          int32
    MyObject: np.ndarray  # (ntarget, maxtiles) int32, unique ridge/crest ID
    WtObject: np.ndarray  # (ntarget, maxtiles) f8, # cells in the wedge
    LnObject: np.ndarray  # (ntarget, maxtiles) f8, # cells in the crest

    # records 7-11 -- centroids on the cubed sphere, and panel
    xwoid: np.ndarray  # wedge-height-weighted centroid, cube x
    ywoid: np.ndarray  # wedge-height-weighted centroid, cube y
    xloid: np.ndarray  # crest centroid, cube x
    yloid: np.ndarray  # crest centroid, cube y
    panel: np.ndarray  # (ntarget, maxtiles) int32, cube panel 1-6

    # records 12-15 -- geographic positions (radians, from CubedSphereRLLFromABP)
    lonc: np.ndarray
    latc: np.ndarray
    lonw: np.ndarray
    latw: np.ndarray

    # records 16-21 -- ridge properties carried over from the ridge database
    aniso: np.ndarray
    anglx: np.ndarray
    angll: np.ndarray
    mxdis: np.ndarray
    hwdth: np.ndarray
    clngt: np.ndarray
    clext: np.ndarray  # crest extent computed in remapridge2tiles

    # record 22 -- accumulated crest length per peak
    clngt2: np.ndarray  # (npeaks,)

    MISSING_REAL = -9999.0
    MISSING_PANEL = -1

    # -- conveniences ----------------------------------------------------
    def valid(self) -> np.ndarray:
        """Boolean (ntarget, maxtiles) mask of slots that hold a real tile."""
        j = np.arange(self.maxtiles)[None, :]
        return j < self.NumObjects[:, None]

    def masked(self, name: str) -> np.ma.MaskedArray:
        """A tile array with unused slots and -9999 fill values masked out."""
        arr = getattr(self, name)
        if arr.shape != (self.ntarget, self.maxtiles):
            raise ValueError(f"{name!r} is not a (ntarget, maxtiles) tile array")
        bad = ~self.valid()
        if np.issubdtype(arr.dtype, np.floating):
            bad |= arr <= self.MISSING_REAL
        else:
            bad |= arr <= self.MISSING_PANEL
        return np.ma.masked_array(arr, mask=bad)

    def cube_cells(self, itarget: int) -> np.ndarray:
        """0-based 1D cube indices belonging to target cell ``itarget`` (0-based).

        The stored ``idxmap`` values are 1-based Fortran indices
        ``ii = (ip-1)*ncube**2 + (iy-1)*ncube + ix``.
        """
        n = int(self.idcoun[itarget])
        return self.idxmap[:n, itarget] - 1

    def __repr__(self) -> str:  # keep it short; the arrays are large
        return (
            f"TileMap(file={self.file!r}, ntarget={self.ntarget}, "
            f"nalloc={self.nalloc}, maxtiles={self.maxtiles}, npeaks={self.npeaks})"
        )

    def __getitem__(self, name):
        return getattr(self, name)

    def __contains__(self, name):
        return hasattr(self, name)

    def keys(self):
        return [k for k in vars(self) if not k.startswith("_")]
# ---------------------------------------------------------------------------
# The translation of `pro rdtilemap`
# ---------------------------------------------------------------------------
def rdtilemap(tile: str, endian: str = "auto", int_kind: int = 4,
              real_kind: int = 8, verbose: bool = False, read_angll = False ) -> TileMap:
    """Read a Ridge_tile_map.dat file.

    Parameters
    ----------
    tile
        Path to the file (the IDL keyword was also called ``tile``).
    endian
        ``"auto"`` (default), ``"<"`` or ``">"``.
    int_kind, real_kind
        Bytes per default Fortran integer and per ``r8`` real.  The defaults
        match a normal build; change ``int_kind`` to 8 if the writer was
        compiled with ``-fdefault-integer-8`` or similar.
    verbose
        Print a short summary of what was read.

    Returns
    -------
    TileMap
    """
    itype = f"i{int_kind}"
    rtype = f"f{real_kind}"

    with F77SequentialReader(tile, endian=endian) as f:
        # record 1 -----------------------------------------------------
        (hdr,) = f.read_arrays((itype, (4,)))
        ntarget, nalloc, maxtiles, npeaks = (int(v) for v in hdr)
        tile_shape = (ntarget, maxtiles)

        # record 2 -----------------------------------------------------
        idxmap, idcoun = f.read_arrays((itype, (nalloc, ntarget)), (itype, (ntarget,)))

        # record 3 -----------------------------------------------------
        (NumObjects,) = f.read_arrays((itype, (ntarget,)))

        # records 4-6 --------------------------------------------------
        (MyObject,) = f.read_arrays((itype, tile_shape))
        (WtObject,) = f.read_arrays((rtype, tile_shape))
        (LnObject,) = f.read_arrays((rtype, tile_shape))

        # records 7-11 -------------------------------------------------
        (xwoid,) = f.read_arrays((rtype, tile_shape))
        (ywoid,) = f.read_arrays((rtype, tile_shape))
        (xloid,) = f.read_arrays((rtype, tile_shape))
        (yloid,) = f.read_arrays((rtype, tile_shape))
        (panel,) = f.read_arrays((itype, tile_shape))

        # records 12-15 ------------------------------------------------
        (lonc,) = f.read_arrays((rtype, tile_shape))
        (latc,) = f.read_arrays((rtype, tile_shape))
        (lonw,) = f.read_arrays((rtype, tile_shape))
        (latw,) = f.read_arrays((rtype, tile_shape))

        # records 16-21/22 ------------------------------------------------
        (aniso,) = f.read_arrays((rtype, tile_shape))
        (anglx,) = f.read_arrays((rtype, tile_shape))
        if (read_angll==True):
            (angll,) = f.read_arrays((rtype, tile_shape))
        (mxdis,) = f.read_arrays((rtype, tile_shape))
        (hwdth,) = f.read_arrays((rtype, tile_shape))
        (clngt,) = f.read_arrays((rtype, tile_shape))
        (clext,) = f.read_arrays((rtype, tile_shape))

        # record 22/23 ----------------------------------------------------
        # The IDL assumed fltarr (4-byte); the Fortran module variable is more
        # likely r8.  Decide from the record length instead of guessing.
        nbytes = f.peek_record_size()
        if npeaks and nbytes == npeaks * 8:
            ctype = "f8"
        elif npeaks and nbytes == npeaks * 4:
            ctype = "f4"
        else:
            raise IOError(
                f"clngt2 record is {nbytes} bytes, which is neither "
                f"{npeaks}*4 nor {npeaks}*8"
            )
        (clngt2,) = f.read_arrays((ctype, (npeaks,)))

    if (read_angll==False):
        angll=np.zeros( tile_shape )

    
    r = TileMap(
        file=tile,
        ntarget=ntarget, nalloc=nalloc, maxtiles=maxtiles, npeaks=npeaks,
        idxmap=idxmap, idcoun=idcoun,
        NumObjects=NumObjects, MyObject=MyObject,
        WtObject=WtObject, LnObject=LnObject,
        xwoid=xwoid, ywoid=ywoid, xloid=xloid, yloid=yloid, panel=panel,
        lonc=lonc, latc=latc, lonw=lonw, latw=latw,
        aniso=aniso, anglx=anglx, angll=angll, mxdis=mxdis, hwdth=hwdth,
        clngt=clngt, clext=clext, clngt2=clngt2,
    )

    if verbose:
        print(f"read {tile}")
        print(f"  ntarget={ntarget}  nalloc={nalloc}  maxtiles={maxtiles}  npeaks={npeaks}")
        print(f"  objects per target cell: min={NumObjects.min()} "
              f"max={NumObjects.max()} mean={NumObjects.mean():.2f}")
        print(f"  cube cells per target cell: min={idcoun.min()} max={idcoun.max()}")
        print(f"  clngt2 stored as {clngt2.dtype}")

    return r


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description=__doc__.splitlines()[3])
    p.add_argument("file", nargs="?", default="output/Ridge_tile_map.dat")
    p.add_argument("--endian", default="auto", choices=["auto", "<", ">"])
    p.add_argument("--int-kind", type=int, default=4)
    args = p.parse_args()

    r = rdtilemap(args.file, endian=args.endian, int_kind=args.int_kind, verbose=True)
    mx = r.masked("mxdis")
    print(f"  mxdis over valid tiles: min={mx.min():.2f} max={mx.max():.2f}")
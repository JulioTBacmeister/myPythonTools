#!/usr/bin/env python3
"""
rdridgelist.py
==============

Reader for the ``Ridge_list_*.dat`` files written by ``find_ridges`` in
``ridge_ana.F90``.

Replaces the hand-unrolled sequence of ``f.read_record(...)`` calls with a
single record table generated from the Fortran write block, so that adding or
reordering a record upstream means editing one line here rather than
re-counting forty-nine calls.

Fortran writer (ridge_ana.F90, find_ridges, unit 31)::

    OPEN (unit = 31, file=ofile, form="UNFORMATTED")
    write(31) npeaks, NSW, PSW
    write(31) xs, ys, MyPanel, NSWx_diag
    write(31) mxvrx
    ... 49 records total ...
    write(31) hnodes_o
    CLOSE(31)

Shapes and orientation
----------------------
Per-ridge vectors are ``(npeaks,)``. The 2-D arrays are declared in Fortran as
``(n, npeaks)`` -- profile index first -- and are returned here transposed to
``(npeaks, n)``, so ``r.rdg_profiles[ipk]`` is one ridge's profile. That
matches how the arrays are normally used from Python and is what the previous
hand-written reader produced.

``n`` is ``PSW+1`` for the node/profile arrays, ``2*PSW+1`` for the extended
ones, and 3 for the wedge arrays.

Usage
-----
::

    from rdridgelist import rdridgelist, latest_ridge_list

    r = rdridgelist(latest_ridge_list(case_dir), verbose=True)
    print(r.npeaks, r.NSW, r.PSW)
    plt.plot(r.mxdis)
    plt.plot(r.rdg_profiles[1234])

By default only the small per-ridge vectors are read; the large profile arrays
are skipped. Pass ``profiles=True`` to load them -- they are several hundred MB
at ncube=3000.
"""

from __future__ import annotations

import glob
import os
import struct
from dataclasses import dataclass, field
from typing import BinaryIO

import numpy as np

__all__ = ["RidgeList", "rdridgelist", "latest_ridge_list", "F77SequentialReader"]


# ---------------------------------------------------------------------------
# Record table -- mirrors the write(31) block in find_ridges, in order.
#
# Each entry is (record_contents, load_by_default) where record_contents is a
# list of (name, dtype_code, shape_key). shape_key is:
#     'v'  -> (npeaks,)
#     'p'  -> (PSW+1,   npeaks)  transposed on return
#     'x'  -> (2*PSW+1, npeaks)  transposed on return
#     'w'  -> (3,       npeaks)  transposed on return
# A record may hold several arrays, e.g. write(31) xs,ys,MyPanel,NSWx_diag.
# ---------------------------------------------------------------------------
RECORDS = [
    ([("xs", "f8", "v"), ("ys", "f8", "v"),
      ("MyPanel", "i4", "v"), ("NSWx_diag", "i4", "v")], True),
    ([("mxvrx", "f8", "v")], True),
    ([("bsvar", "f8", "v")], True),
    ([("mxdis", "f8", "v")], True),
    ([("anglx", "f8", "v")], True),
    ([("aniso", "f8", "v")], True),
    ([("mnslp", "f8", "v")], True),
    ([("angll", "f8", "v")], True),
    ([("xspk", "f8", "v"), ("xs01", "f8", "v")], True),
    ([("yspk", "f8", "v"), ("ys01", "f8", "v")], True),
    ([("mxds0", "f8", "v")], True),
    ([("mxds1", "f8", "v")], True),
    ([("sft0", "f8", "v")], True),
    ([("sft1", "f8", "v")], True),
    ([("hwdth", "f8", "v")], True),
    ([("npks", "f8", "v")], True),     # count of peaks in profile; REAL in Fortran
    ([("mxvry", "f8", "v")], True),
    ([("nvls", "f8", "v")], True),
    ([("pkhts", "f8", "v")], True),
    ([("vldps", "f8", "v")], True),
    ([("rwpks", "f8", "v")], True),
    ([("rwvls", "f8", "v")], True),
    ([("lon0", "f8", "v")], True),
    ([("lon1", "f8", "v")], True),
    ([("lat0", "f8", "v")], True),
    ([("lat1", "f8", "v")], True),
    ([("uniqid", "f8", "v")], True),
    ([("riseq", "f8", "v")], True),
    ([("fallq", "f8", "v")], True),
    ([("clngth", "f8", "v")], True),
    ([("mxds2", "f8", "v")], True),
    ([("rdg_profiles", "f8", "p")], True),
    ([("crst_profiles", "f8", "p")], True),
    ([("crst_silhous", "f8", "p")], False),
    ([("isoht", "f8", "v")], True),
    ([("isowd", "f8", "v")], True),
    ([("isobs", "f8", "v")], True),
    ([("RefFac", "f8", "v")], True),
    ([("rdg_profiles_x", "f8", "x")], True),
    ([("xnodes_list", "i4", "p")], True),
    ([("hnodes_list", "f8", "p")], True),
    ([("dcenter_list", "i4", "v")], True),
    ([("xwedge_list", "i4", "w"), ("hwedge_list", "f8", "w")], True),
    ([("nnodes_list", "i4", "v")], True),
    ([("hwedge_x", "f8", "x")], True),
    ([("hnodes_x", "f8", "x")], True),
    ([("hwedge_o", "f8", "x")], True),
    ([("hnodes_o", "f8", "x")], True),
    ([("lonc", "f8", "v")], True),
    ([("latc", "f8", "v")], True),
]
OPTIONAL_RECORDS = []

class F77SequentialReader:
    """Fortran ``form="UNFORMATTED"`` sequential reader with byte-order sniffing.

    Handles gfortran's negative-marker subrecords, which ``scipy.io.FortranFile``
    does not. Records here stay under 2 GiB at ncube=3000, but that depends on
    npeaks and is not guaranteed.
    """

    def __init__(self, path: str, endian: str = "auto"):
        self.path = path
        self._fh: BinaryIO = open(path, "rb")
        self.endian = self._detect_endian() if endian == "auto" else endian

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def close(self):
        if not self._fh.closed:
            self._fh.close()

    def _detect_endian(self) -> str:
        head = self._fh.read(4)
        if len(head) < 4:
            raise IOError(f"{self.path}: too short to be a Fortran unformatted file")
        for e in ("<", ">"):
            (n,) = struct.unpack(e + "i", head)
            if n <= 0 or n > 1 << 34:
                continue
            self._fh.seek(4 + n)
            tail = self._fh.read(4)
            if len(tail) == 4 and struct.unpack(e + "i", tail)[0] == n:
                self._fh.seek(0)
                return e
        self._fh.seek(0)
        raise IOError(f"{self.path}: record markers not self-consistent in either byte order")

    def _marker(self) -> int:
        raw = self._fh.read(4)
        if len(raw) < 4:
            raise EOFError("end of file while reading a record marker")
        return struct.unpack(self.endian + "i", raw)[0]

    def peek_size(self) -> int:
        """Payload length of the next record, without consuming it."""
        here = self._fh.tell()
        total = 0
        try:
            while True:
                head = self._marker()
                total += abs(head)
                self._fh.seek(abs(head), 1)
                self._marker()
                if head >= 0:
                    return total
        finally:
            self._fh.seek(here)

    def skip_record(self) -> int:
        n = self.peek_size()
        while True:
            head = self._marker()
            self._fh.seek(abs(head), 1)
            self._marker()
            if head >= 0:
                return n

    def read_record(self) -> bytes:
        chunks = []
        while True:
            head = self._marker()
            payload = self._fh.read(abs(head))
            if len(payload) != abs(head):
                raise EOFError("end of file inside a record payload")
            chunks.append(payload)
            tail = self._marker()
            if tail != head:
                raise IOError(f"record markers disagree ({head} vs {tail}) in {self.path}")
            if head >= 0:
                return b"".join(chunks)


@dataclass
class RidgeList:
    """Contents of a Ridge_list_*.dat file. Fields are attached dynamically."""
    file: str
    npeaks: int
    NSW: int
    PSW: int
    loaded: list = field(default_factory=list)
    skipped: list = field(default_factory=list)

    def __repr__(self):
        return (f"RidgeList(file={os.path.basename(self.file)!r}, npeaks={self.npeaks}, "
                f"NSW={self.NSW}, PSW={self.PSW}, loaded={len(self.loaded)} arrays)")

    def __getitem__(self, name):
        return getattr(self, name)

    def __contains__(self, name):
        return hasattr(self, name)

    def keys(self):
        return [k for k in vars(self) if not k.startswith("_")]

    def valid(self):
        """Mask of ridges that survived thinout (mxdis is set to a negative
        flag for rejected peaks: -999, -222, -333, -444)."""
        return self.mxdis > 0.0


def latest_ridge_list(case_dir: str) -> str:
    """Most recently modified Ridge_list_*.dat under ``case_dir/output``.

    The original code used ``glob.glob(...)[-1]``, but glob does not sort --
    it returns directory order, so the 'last' entry was arbitrary. This picks
    by modification time.
    """
    pat = os.path.join(case_dir, "output", "Ridge_list*dat")
    hits = glob.glob(pat)
    if not hits:
        raise FileNotFoundError(f"no files matching {pat}")
    return max(hits, key=os.path.getmtime)


def rdridgelist(path: str, profiles: bool = False, endian: str = "auto",
                verbose: bool = False) -> RidgeList:
    """Read a Ridge_list_*.dat file.

    Parameters
    ----------
    path
        File to read.
    profiles
        Load the large 2-D profile/node arrays as well. Off by default: they
        are several hundred MB at ncube=3000.
    endian
        "auto" (default), "<" or ">".
    verbose
        Print a short summary.
    """
    with F77SequentialReader(path, endian=endian) as f:
        hdr = np.frombuffer(f.read_record(), dtype=np.dtype("i4").newbyteorder(f.endian), count=3)
        npeaks, NSW, PSW = (int(v) for v in hdr)

        shapes = {"v": (npeaks,),
                  "p": (PSW + 1, npeaks),
                  "x": (2 * PSW + 1, npeaks),
                  "w": (3, npeaks)}

        r = RidgeList(file=path, npeaks=npeaks, NSW=NSW, PSW=PSW)

        for entries, default, optional in ([(e, d, False) for e, d in RECORDS] +
                                           [(e, d, True) for e, d in OPTIONAL_RECORDS]):
            want = default or profiles
            expected = sum(int(np.prod(shapes[sk])) * np.dtype(dt).itemsize
                           for _, dt, sk in entries)
            try:
                got = f.peek_size()
            except EOFError:
                if optional:
                    break                       # older file without the trailing records
                raise
            if got != expected:
                names = ", ".join(n for n, _, _ in entries)
                raise IOError(
                    f"{path}: record for [{names}] is {got} bytes, expected {expected}. "
                    "The write order in find_ridges has probably changed; "
                    "update RECORDS in this module to match."
                )
            if not want:
                f.skip_record()
                r.skipped.extend(n for n, _, _ in entries)
                continue

            buf, off = f.read_record(), 0
            for name, dt, sk in entries:
                d = np.dtype(dt).newbyteorder(f.endian)
                shp = shapes[sk]
                cnt = int(np.prod(shp))
                a = np.frombuffer(buf, dtype=d, count=cnt, offset=off)
                a = np.reshape(a, shp, order="F").astype(d.newbyteorder("="), copy=True)
                if a.ndim == 2:
                    a = np.ascontiguousarray(a.T)     # -> (npeaks, n)
                setattr(r, name, a)
                r.loaded.append(name)
                off += cnt * d.itemsize

        trailing = 0
        try:
            while True:
                f.skip_record()
                trailing += 1
        except EOFError:
            pass
        if trailing:
            raise IOError(f"{path}: {trailing} unexpected record(s) after the last "
                          "known one -- find_ridges has written something new.")

    if verbose:
        print(f"read {path}")
        print(f"  npeaks={npeaks}  NSW={NSW}  PSW={PSW}")
        print(f"  loaded {len(r.loaded)} arrays, skipped {len(r.skipped)}"
              f"{' (pass profiles=True to load)' if not profiles else ''}")
        v = r.valid()
        print(f"  ridges surviving thinout: {int(v.sum())} of {npeaks} "
              f"({100*v.mean():.1f}%)")
        print(f"  mxdis  (valid): min={r.mxdis[v].min():.1f} max={r.mxdis[v].max():.1f} "
              f"mean={r.mxdis[v].mean():.1f}")
        print(f"  clngth (valid): min={r.clngth[v].min():.1f} max={r.clngth[v].max():.1f} "
              f"mean={r.clngth[v].mean():.1f}  [cube cells]")

    return r


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Summarise a Ridge_list_*.dat file")
    p.add_argument("path", help="file, or a case directory to search")
    p.add_argument("--profiles", action="store_true", help="also load the big 2-D arrays")
    a = p.parse_args()
    target = a.path if os.path.isfile(a.path) else latest_ridge_list(a.path)
    rdridgelist(target, profiles=a.profiles, verbose=True)

#!/usr/bin/env python3
"""
rdremap.py
==========

Reader for the ``remap_nc*_vX_*.dat`` files written by ``remapridge2cube`` in
``ridge_ana.F90``. These hold the ridge quantities *painted back onto the
high-resolution cubed sphere*, before any remapping to a target grid.

Companion to ``rdridgelist.py`` (the per-ridge database) and ``rdtilemap.py``
(the per-target-cell tiles). Same calling convention as both.

Fortran writer (ridge_ana.F90, remapridge2cube, unit 911)::

    write(911) ncube, npeaks
    write(911) uniqidC          ! crest ID painted along the ridgeline
    write(911) anisoC
    ... 19 records total ...
    write(911) repntC

MEMORY WARNING
--------------
Every field is ``ncube*ncube*6`` doubles: **432 MB each at ncube=3000**, and
there are 17 of them, so the whole file is ~7.4 GB. Fields are therefore
loaded only on request::

    r = rdremap(path)                        # header + peak coords only
    r = rdremap(path, fields=["uniqid"])     # one field
    r = rdremap(path, fields="all")          # everything - needs ~7.4 GB

Shapes
------
The Fortran arrays are flat, indexed ``ii = (ip-1)*ncube**2 + (iy-1)*ncube + ix``.
They are returned as ``(6, ncube, ncube)`` in C order, i.e. ``[panel, iy, ix]``,
with panel 0-based and iy/ix 0-based.

The peak-coordinate record (``xs, ys, xspk, yspk, peaks%i, peaks%j``) is small
and always read.

Field names drop the trailing 'C' of the Fortran names: ``uniqidC`` ->
``uniqid``, and so on.

Usage
-----
::

    from rdremap import rdremap, latest_remap, count_ids

    r = rdremap(latest_remap(case_dir), fields=["uniqid","uniqwg"], verbose=True)
    n = count_ids(r.uniqid, npeaks)      # cube cells per ridge ID
"""

from __future__ import annotations

import glob
import os
import struct
from dataclasses import dataclass, field
from typing import BinaryIO, Sequence

import numpy as np

__all__ = ["Remap", "rdremap", "latest_remap", "count_ids", "F77SequentialReader"]


# Record table, in write order. 'C' = a full cubed-sphere field,
# 'P' = the per-peak coordinate record.
CUBE_FIELDS = [
    "uniqid",   # crest ID, painted 1 cell wide along the ridgeline
    "aniso",
    "anglx",    # ridge orientation in (alpha,beta), degrees
    "mxdis",    # obstacle height, m
    "hwdth",
    "clngt",
    "block",    # fleshout_block
    "profi",    # fleshout_profi with rdg_profiles_x
    "nodes",    # fleshout_profi with hnodes_x
    "wedge",    # fleshout_profi with hwedge_x
    "nodos",    # fleshout_profi with hnodes_o
    "wedgo",    # fleshout_profi with hwedge_o  <- weights for xwoid/ywoid
    None,       # placeholder: the peak-coordinate record sits here
    "riseq",
    "fallq",
    "uniqwg",   # wedge ID, painted over the full footprint by color_on_profi
    "wedgi",
    "repnt",
]

PEAK_RECORD = [("xs", "f8"), ("ys", "f8"), ("xspk", "f8"), ("yspk", "f8"),
               ("peaks_i", "i4"), ("peaks_j", "i4")]


class F77SequentialReader:
    """Fortran sequential unformatted reader with byte-order sniffing.

    Handles gfortran subrecords (negative leading marker), which
    ``scipy.io.FortranFile`` does not. Records here are ~432 MB at ncube=3000,
    comfortably under the 2 GiB threshold, but that scales with ncube.
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
class Remap:
    """Painted cubed-sphere fields from remapridge2cube."""
    file: str
    ncube: int
    npeaks: int
    loaded: list = field(default_factory=list)
    skipped: list = field(default_factory=list)

    def __repr__(self):
        return (f"Remap(file={os.path.basename(self.file)!r}, ncube={self.ncube}, "
                f"npeaks={self.npeaks}, loaded={self.loaded})")

    def __getitem__(self, name):
        return getattr(self, name)

    def __contains__(self, name):
        return hasattr(self, name)

    def keys(self):
        return [k for k in vars(self) if not k.startswith("_")]

    @property
    def ncells(self) -> int:
        return 6 * self.ncube * self.ncube


def latest_remap(case_dir: str) -> str:
    """Most recently modified remap*dat under ``case_dir/output``.

    ``glob`` returns directory order, not sorted order, so ``glob(...)[-1]``
    picks an arbitrary file when several runs are present. This uses mtime.
    """
    pat = os.path.join(case_dir, "output", "remap*dat")
    hits = glob.glob(pat)
    if not hits:
        raise FileNotFoundError(f"no files matching {pat}")
    return max(hits, key=os.path.getmtime)


def count_ids(idfield: np.ndarray, npeaks: int) -> np.ndarray:
    """Cube cells carrying each ridge ID.

    Returns an array of length npeaks; element k-1 is the number of cells whose
    ID is k (IDs are 1-based, uniqid(ipk) = 1.0*ipk).

    Truncation, not rounding, is used, matching the Fortran:
        IdxOc(i) = uniqidC(...)      ! real -> integer assignment
    and only strictly positive values are counted, matching
        IdxPack = pack( IdxOc, (IdxOc>0) )
    """
    ids = np.asarray(idfield).astype(np.int64).ravel()
    ids = np.where(ids > 0, ids, 0)
    return np.bincount(ids, minlength=npeaks + 1)[1:]


def rdremap(path: str, fields: Sequence[str] | str = (), endian: str = "auto",
            verbose: bool = False) -> Remap:
    """Read a remap_*.dat file.

    Parameters
    ----------
    path
        File to read.
    fields
        Names of cubed-sphere fields to load (see CUBE_FIELDS), or "all".
        Default loads none -- each is 432 MB at ncube=3000.
    endian
        "auto" (default), "<" or ">".
    verbose
        Print a summary.
    """
    names = [f for f in CUBE_FIELDS if f]
    if fields == "all":
        want = set(names)
    else:
        want = set(fields)
        bad = want - set(names)
        if bad:
            raise ValueError(f"unknown field(s) {sorted(bad)}; available: {names}")

    with F77SequentialReader(path, endian=endian) as f:
        i4 = np.dtype("i4").newbyteorder(f.endian)
        f8 = np.dtype("f8").newbyteorder(f.endian)

        hdr = np.frombuffer(f.read_record(), dtype=i4, count=2)
        ncube, npeaks = int(hdr[0]), int(hdr[1])
        ncells = 6 * ncube * ncube
        r = Remap(file=path, ncube=ncube, npeaks=npeaks)

        for name in CUBE_FIELDS:
            if name is None:
                # xs, ys, xspk, yspk, peaks%i, peaks%j -- all length npeaks
                expected = npeaks * (4 * 8 + 2 * 4)
                got = f.peek_size()
                if got != expected:
                    raise IOError(
                        f"{path}: peak-coordinate record is {got} bytes, expected "
                        f"{expected}. Write order in remapridge2cube has changed."
                    )
                buf, off = f.read_record(), 0
                for nm, dt in PEAK_RECORD:
                    d = np.dtype(dt).newbyteorder(f.endian)
                    a = np.frombuffer(buf, dtype=d, count=npeaks, offset=off)
                    setattr(r, nm, a.astype(d.newbyteorder("="), copy=True))
                    off += npeaks * d.itemsize
                continue

            expected = ncells * 8
            got = f.peek_size()
            if got != expected:
                raise IOError(
                    f"{path}: record for '{name}' is {got} bytes, expected {expected} "
                    f"(6*{ncube}^2*8). Write order in remapridge2cube has changed; "
                    "update CUBE_FIELDS."
                )
            if name not in want:
                f.skip_record()
                r.skipped.append(name)
                continue
            a = np.frombuffer(f.read_record(), dtype=f8, count=ncells)
            a = a.astype(f8.newbyteorder("="), copy=True).reshape(6, ncube, ncube)
            setattr(r, name, a)
            r.loaded.append(name)

        trailing = 0
        try:
            while True:
                f.skip_record()
                trailing += 1
        except EOFError:
            pass
        if trailing:
            raise IOError(f"{path}: {trailing} unexpected record(s) after 'repnt' -- "
                          "remapridge2cube has written something new.")

    if verbose:
        gb = len(r.loaded) * ncells * 8 / 2**30
        print(f"read {path}")
        print(f"  ncube={ncube}  npeaks={npeaks}  cells={ncells:,}")
        print(f"  loaded: {r.loaded or '(none)'}  ({gb:.2f} GB)")
        print(f"  skipped {len(r.skipped)} field(s)")
        if "uniqid" in r.loaded:
            n = count_ids(r.uniqid, npeaks)
            print(f"  uniqid: {int((n>0).sum()):,} of {npeaks:,} ridges have painted "
                  f"crest cells; {int(n.sum()):,} cells total")

    return r


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Summarise a remap_*.dat file")
    p.add_argument("path", help="file, or a case directory to search")
    p.add_argument("--fields", default="uniqid",
                   help="comma-separated field list, or 'all', or 'none'")
    a = p.parse_args()
    target = a.path if os.path.isfile(a.path) else latest_remap(a.path)
    fl = () if a.fields == "none" else (a.fields if a.fields == "all"
                                        else a.fields.split(","))
    rdremap(target, fields=fl, verbose=True)
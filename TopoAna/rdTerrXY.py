#!/usr/bin/env python3
"""
rdterrxy.py
===========

Reader for the ``TerrXY_list_*.dat`` files written by ``find_ridges`` in
``ridge_ana.F90`` when the ``DEBUGOUTPUT`` cpp flag is defined.

Each record holds the topography sub-window used for one peak's ridge analysis,
both before and after rotation into the diagnosed crest-normal frame. This is
what you need to draw the ridge-finding process: the raw terrain block, the
same block rotated to the ridge orientation, and the cross-ridge profile that
falls out of averaging the rotated block along the crest.

Fortran writer (subset version)::

    write(32) ndiag , PSW , npeaks        ! ndiag records follow, NOT npeaks
    do idg = 1, ndiag
       ipk = idiag_list(idg)
       write(32) ipk, xs(ipk), ys(ipk), MyPanel(ipk), xspk(ipk), yspk(ipk)
       write(32) rt_diag(:,:,idg)         ! rotated window, (2*PSW+1, 2*PSW+1)
       write(32) suba_diag(:,:,idg)       ! unrotated window, same shape
    end do

The peak id is stored in every record, so entries join to ``Ridge_list`` by ID
rather than by position. That matters: with the subset patch only peaks inside
the chosen box are written, so record n is NOT peak n.

The ORIGINAL (pre-subset) format wrote ``npeaks, PSW`` as the header and had no
``ipk`` in the per-peak record. This reader detects which it is being given and
handles both.

Usage
-----
::

    from rdterrxy import rdterrxy, latest_terrxy
    from rdridgelist import rdridgelist, latest_ridge_list

    t = rdterrxy(latest_terrxy(case_dir), verbose=True)
    r = rdridgelist(latest_ridge_list(case_dir))

    k = t.index_of(123456)          # storage slot for peak 123456, or None
    win = t.suba[k]                 # unrotated terrain window
    rot = t.rt[k]                   # rotated into the crest frame
    prof = t.profile(k)             # cross-ridge profile from the rotated block

    # everything about that peak from the ridge list, by ID
    ipk = t.ipk[k]
    print(r.mxdis[ipk-1], r.anglx[ipk-1], r.clngth[ipk-1])
"""

from __future__ import annotations

import glob
import os
import struct
from dataclasses import dataclass, field
from typing import BinaryIO

import numpy as np

__all__ = ["TerrXY", "rdterrxy", "latest_terrxy", "F77SequentialReader"]


class F77SequentialReader:
    """Fortran sequential unformatted reader with byte-order sniffing.

    Handles gfortran's negative-marker subrecords. Same helper used by
    rdridgelist / rdremap / rdtilemap.
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
        raise IOError(f"{self.path}: record markers not self-consistent")

    def _marker(self) -> int:
        raw = self._fh.read(4)
        if len(raw) < 4:
            raise EOFError("end of file at a record marker")
        return struct.unpack(self.endian + "i", raw)[0]

    def peek_size(self) -> int:
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

    def at_eof(self) -> bool:
        here = self._fh.tell()
        b = self._fh.read(4)
        self._fh.seek(here)
        return len(b) < 4

    def read_record(self) -> bytes:
        chunks = []
        while True:
            h = self._marker()
            b = self._fh.read(abs(h))
            if len(b) != abs(h):
                raise EOFError("end of file inside a record payload")
            chunks.append(b)
            if self._marker() != h:
                raise IOError(f"record markers disagree in {self.path}")
            if h >= 0:
                return b"".join(chunks)


@dataclass
class TerrXY:
    """Per-peak analysis windows from a TerrXY_list file."""
    file: str
    PSW: int
    npeaks: int                 # total peaks in the run (not the number stored)
    ipk: np.ndarray             # (ndiag,) peak id of each stored window, 1-based
    xs: np.ndarray              # (ndiag,) window-centre cube x
    ys: np.ndarray              # (ndiag,) window-centre cube y
    panel: np.ndarray           # (ndiag,) cube panel 1-6
    xspk: np.ndarray            # (ndiag,) diagnosed crest x
    yspk: np.ndarray            # (ndiag,) diagnosed crest y
    rt: np.ndarray              # (ndiag, 2*PSW+1, 2*PSW+1) rotated window
    suba: np.ndarray            # (ndiag, 2*PSW+1, 2*PSW+1) unrotated window
    has_ids: bool = True        # False for the old format, where ipk is inferred
    _lookup: dict = field(default_factory=dict, repr=False)

    @property
    def ndiag(self) -> int:
        return self.ipk.size

    @property
    def nwin(self) -> int:
        return 2 * self.PSW + 1

    def __repr__(self):
        return (f"TerrXY(file={os.path.basename(self.file)!r}, ndiag={self.ndiag}, "
                f"window={self.nwin}x{self.nwin}, PSW={self.PSW}, "
                f"npeaks={self.npeaks}, has_ids={self.has_ids})")

    def __getitem__(self, name):
        return getattr(self, name)
    ##### added
    def __contains__(self, name):
        return hasattr(self, name)

    def keys(self):
        return [k for k in vars(self) if not k.startswith("_")]
    
    def index_of(self, peak_id: int):
        """Storage slot holding peak `peak_id`, or None if it was not retained."""
        return self._lookup.get(int(peak_id))

    def profile(self, k: int, ns0: int = None, ns1: int = None) -> np.ndarray:
        """Cross-ridge profile for slot k: the rotated window averaged along the
        crest, which is what ANISO_ANA computes as RTX.

        The Fortran averages over the inner sub-block ns0:ns1-1 of a
        (2*nsw+1) array, with ns0 = nsw/2+1 and ns1 = ns0+nsw+1. Here the
        window is (2*PSW+1) and nsw is not stored in the file, so by default
        the central half of the window is used. Pass ns0/ns1 to match the
        Fortran exactly if you know nsw.
        """
        a = self.rt[k]
        n = a.shape[0]
        if ns0 is None:
            ns0 = n // 4
        if ns1 is None:
            ns1 = n - n // 4
        return a[:, ns0:ns1].mean(axis=1)

    def valid_mask(self, k: int, fill: float = -9999999.0) -> np.ndarray:
        """Cells of the rotated window that came from inside the source block.
        rotby4 marks points outside the input with a large negative value."""
        return self.rt[k] > fill


def latest_terrxy(case_dir: str) -> str:
    """Most recently modified TerrXY_list_*.dat under case_dir/output.

    glob returns directory order, not sorted order, so pick by mtime.
    """
    pat = os.path.join(case_dir, "output", "TerrXY_list*dat")
    hits = glob.glob(pat)
    if not hits:
        raise FileNotFoundError(f"no files matching {pat}\n"
                                "  (was ridge_ana.F90 built with -DDEBUGOUTPUT?)")
    return max(hits, key=os.path.getmtime)


def rdterrxy(path: str, endian: str = "auto", max_records: int = None,
             verbose: bool = False) -> TerrXY:
    """Read a TerrXY_list_*.dat file.

    Parameters
    ----------
    path : str
    max_records : int or None
        Stop after this many windows. Each is 2*(2*PSW+1)^2*8 bytes, about
        258 kB at PSW=63, so a few thousand is fine, but use this if the file
        is large and you only want a handful for a figure.
    """
    with F77SequentialReader(path, endian=endian) as f:
        i4 = np.dtype("i4").newbyteorder(f.endian)
        f8 = np.dtype("f8").newbyteorder(f.endian)

        hdr_bytes = f.read_record()
        hdr = np.frombuffer(hdr_bytes, dtype=i4)
        if hdr.size == 3:
            ndiag, PSW, npeaks = (int(v) for v in hdr)
            has_ids = True
        elif hdr.size == 2:
            # original format: header was (npeaks, PSW) and records carried no id
            npeaks, PSW = int(hdr[0]), int(hdr[1])
            ndiag, has_ids = npeaks, False
        else:
            raise IOError(f"{path}: header record has {hdr.size} integers, "
                          "expected 2 (original) or 3 (subset format)")

        nwin = 2 * PSW + 1
        win_bytes = nwin * nwin * 8
        # per-peak scalar record: [ipk] + xs,ys + MyPanel + xspk,yspk
        scal_bytes = (4 if has_ids else 0) + 8 + 8 + 4 + 8 + 8

        ipk, xs, ys, pn, xspk, yspk = [], [], [], [], [], []
        rt, suba = [], []

        n = 0
        while not f.at_eof():
            if max_records is not None and n >= max_records:
                break
            got = f.peek_size()
            if got != scal_bytes:
                raise IOError(
                    f"{path}: record {n} scalar block is {got} bytes, expected "
                    f"{scal_bytes}. Header said PSW={PSW}; the write order in "
                    "find_ridges may have changed.")
            b = f.read_record()
            off = 0
            if has_ids:
                ipk.append(int(np.frombuffer(b, i4, 1, off)[0])); off += 4
            xs.append(float(np.frombuffer(b, f8, 1, off)[0])); off += 8
            ys.append(float(np.frombuffer(b, f8, 1, off)[0])); off += 8
            pn.append(int(np.frombuffer(b, i4, 1, off)[0])); off += 4
            xspk.append(float(np.frombuffer(b, f8, 1, off)[0])); off += 8
            yspk.append(float(np.frombuffer(b, f8, 1, off)[0])); off += 8

            for store in (rt, suba):
                got = f.peek_size()
                if got != win_bytes:
                    raise IOError(f"{path}: window record is {got} bytes, expected "
                                  f"{win_bytes} ({nwin}x{nwin} f8)")
                a = np.frombuffer(f.read_record(), dtype=f8, count=nwin*nwin)
                # Fortran (2*PSW+1, 2*PSW+1) column-major
                store.append(np.reshape(a, (nwin, nwin), order="F")
                             .astype(f8.newbyteorder("="), copy=True))
            n += 1

    if not has_ids:
        ipk = list(range(1, n + 1))     # original format: record k IS peak k

    t = TerrXY(file=path, PSW=PSW, npeaks=npeaks,
               ipk=np.array(ipk, dtype=np.int64),
               xs=np.array(xs), ys=np.array(ys),
               panel=np.array(pn, dtype=np.int32),
               xspk=np.array(xspk), yspk=np.array(yspk),
               rt=np.array(rt), suba=np.array(suba),
               has_ids=has_ids)
    t._lookup = {int(p): k for k, p in enumerate(t.ipk)}

    if verbose:
        print(f"read {path}")
        print(f"  format: {'subset (ids stored)' if has_ids else 'original (positional)'}")
        print(f"  windows stored: {t.ndiag:,}"
              + (f" of {npeaks:,} peaks in the run" if has_ids else ""))
        print(f"  window size: {nwin} x {nwin}   PSW={PSW}")
        print(f"  memory: {(t.rt.nbytes + t.suba.nbytes)/2**20:.1f} MB")
        if t.ndiag:
            print(f"  peak ids: {t.ipk.min()} to {t.ipk.max()}")
            print(f"  panels:   {sorted(set(t.panel.tolist()))}")
            print(f"  cube x:   {t.xs.min():.0f} to {t.xs.max():.0f}")
            print(f"  cube y:   {t.ys.min():.0f} to {t.ys.max():.0f}")
    return t


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Summarise a TerrXY_list_*.dat file")
    p.add_argument("path", help="file, or a case directory to search")
    p.add_argument("--max", type=int, default=None, help="stop after N windows")
    a = p.parse_args()
    target = a.path if os.path.isfile(a.path) else latest_terrxy(a.path)
    rdterrxy(target, max_records=a.max, verbose=True)
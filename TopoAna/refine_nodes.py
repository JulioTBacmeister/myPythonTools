#!/usr/bin/env python3
"""
refine_nodes.py -- add one best-fit node between each pair of rnodes nodes.

rnodes gives the local minima and maxima of a ridge profile. The polyline
through those nodes is a crude description of the profile: between two nodes it
is a straight segment, so any bow in the flank is invisible. This inserts ONE
extra node in each interval, positioned in both x and h so that the resulting
two-segment polyline minimises the rms error against the profile over that
interval.

For a fixed insertion position xm the polyline is LINEAR in the inserted height
hm, so hm has a closed form and only xm has to be searched:

    poly(x) = a(x) + b(x)*hm
    hm      = sum( b * (profile - a) ) / sum( b*b )

Positions are 1-BASED throughout, matching rnodes and the Fortran.

    from refine_nodes import refine_segment, refine_profile
    out = refine_profile(profile, xnodes, nnode)
    out['x'], out['h']          # refined node list, original + inserted
    out['inserted']             # bool, one per refined node
    out['rms_before'], out['rms_after']
"""
import numpy as np

__all__ = ["refine_segment", "refine_profile"]


def refine_segment(prof, x0, x1):
    """
    Best single inserted node between 1-based positions x0 and x1.

    Returns (xm, hm, rms_before, rms_after), or None if the interval is too
    short to hold an interior point.
    """
    prof = np.asarray(prof, dtype=float)
    h0, h1 = prof[x0 - 1], prof[x1 - 1]
    if x1 - x0 < 2:
        return None

    xs = np.arange(x0, x1 + 1)              # 1-based sample positions
    ps = prof[x0 - 1:x1]                    # the profile over the interval

    # straight segment, for reference
    straight = h0 + (h1 - h0) * (xs - x0) / (x1 - x0)
    rms_before = float(np.sqrt(np.mean((ps - straight) ** 2)))

    best = None
    for xm in range(x0 + 1, x1):            # interior positions only
        left = xs <= xm
        # poly = a + b*hm on each side
        a = np.empty_like(xs, dtype=float)
        b = np.empty_like(xs, dtype=float)

        tl = (xs[left] - x0) / (xm - x0)
        a[left] = h0 * (1.0 - tl)
        b[left] = tl

        tr = (x1 - xs[~left]) / (x1 - xm)
        a[~left] = h1 * (1.0 - tr)
        b[~left] = tr

        bb = float(np.dot(b, b))
        if bb <= 0.0:
            continue
        hm = float(np.dot(b, ps - a) / bb)
        resid = ps - (a + b * hm)
        rms = float(np.sqrt(np.mean(resid ** 2)))
        if best is None or rms < best[2]:
            best = (xm, hm, rms)

    if best is None:
        return None
    return best[0], best[1], rms_before, best[2]


def refine_profile(prof, xnodes, nnode=None, hnodes=None):
    """
    Insert one best-fit node in every interval of the rnodes node list.

    prof   : 1-D profile (the rdg_profiles_x segment)
    xnodes : either the 1-based node positions from rnodes (padded with -1 is
             fine), or the whole dict rnodes() returned, in which case nnode
             is taken from it
    nnode  : how many positions are real; required only if xnodes is an array
    hnodes : ignored, accepted so the rnodes outputs can be passed straight
             through. Heights are read from prof at xnodes, which is what
             rnodes does, so passing them changes nothing.

    Returns a dict with
        x, h       : the refined node list, originals and inserted interleaved
        inserted   : bool per node, True for the ones this routine added
        rms_before : rms of the original node polyline against prof
        rms_after  : rms of the refined polyline against prof
    """
    if isinstance(xnodes, dict):                 # the rnodes() output
        out = xnodes
        nnode = int(out["nnode"])
        xnodes = out["xnodes"]
    if nnode is None:
        raise ValueError("nnode is required when xnodes is an array")

    prof = np.asarray(prof, dtype=float)
    xn = np.asarray(xnodes, dtype=int)[:nnode]

    xs_out, hs_out, ins_out = [], [], []
    rb, ra, nseg = [], [], 0

    for k in range(nnode):
        xs_out.append(int(xn[k]))
        hs_out.append(float(prof[xn[k] - 1]))
        ins_out.append(False)
        if k == nnode - 1:
            break
        got = refine_segment(prof, int(xn[k]), int(xn[k + 1]))
        if got is None:
            continue
        xm, hm, r0, r1 = got
        xs_out.append(int(xm))
        hs_out.append(float(hm))
        ins_out.append(True)
        rb.append(r0)
        ra.append(r1)
        nseg += 1

    x = np.array(xs_out, dtype=int)
    h = np.array(hs_out, dtype=float)
    ins = np.array(ins_out, dtype=bool)

    # whole-profile rms, original node polyline vs refined polyline
    span = np.arange(int(xn[0]), int(xn[-1]) + 1)
    p = prof[span - 1]
    before = np.interp(span, xn, prof[xn - 1])
    after = np.interp(span, x, h)
    return {
        "x": x, "h": h, "inserted": ins,
        "rms_before": float(np.sqrt(np.mean((p - before) ** 2))),
        "rms_after": float(np.sqrt(np.mean((p - after) ** 2))),
        "seg_rms_before": np.array(rb), "seg_rms_after": np.array(ra),
        "nseg": nseg,
    }

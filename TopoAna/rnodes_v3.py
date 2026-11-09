"""
rnodes_v3.py  --  staircase wedge: rms_wedge + Julio's slope-sign refinement

    from rnodes_v3 import rnodes_v3
    w = rnodes_v3(ridge)                       # ridge: the nsw+1 real profile points

    w['xnodes_L'], w['hnodes_L']   # outer->inner nodes on the left : [x_L, ..., x_0]
    w['xnodes_R'], w['hnodes_R']   # inner->outer nodes on the right: [x_1, ..., x_R]
    w['xwedge'],   w['hwedge']     # full polyline: left nodes + x_c + right nodes
    w['disp_L'], w['wdth_L']       # obstacle height / half-width on the left  = h_c - h(foot_L), x_c - foot_L
    w['disp_R'], w['wdth_R']
    w['mxdis'], w['hwdth'], w['riseq'], w['fallq']
    w['x0'], w['x1']               # the FIRST knees (same as rms_wedge's 3-point wedge)
    w['foot_L'], w['foot_R']       # where the height is measured to
    w['ok']

Algorithm
  Steps 1-3 as in rms_wedge (central peak, active window [x_L,x_R], best
  5-point wedge [x_L, x_0, x_c, x_1, x_R] by rms over the profile).
  Step 4 (per side, shown for the left):
      outer, inner = x_L, x_0
      loop:
         s_out = slope of fitted segment outer->inner
         s_in  = slope of fitted segment inner->(next node toward the crest)
         if sign(s_out) != sign(s_in)  (case A: a genuine foot)       -> STOP, foot = inner
         if no interior cell between outer and inner                  -> STOP, foot = outer
         else (case B: ramp continues) insert a node k in (outer,inner)
              choosing k to minimise the rms of the polyline over [x_L, x_c]
              with all other nodes fixed; then inner = k and repeat.
      disp_L = h_c - h(foot), wdth_L = x_c - foot
  A zero outer slope (flat apron) counts as a sign change (stop).

Positions are 1-based, like rnodes.py / the Fortran.
"""
import numpy as np

__all__ = ["rnodes_v3", "plot_wedge"]


def _central_peak(h, ns2):
    n = len(h)
    best = -1; dmin = None
    for i in range(2, n):                     # 1-based 2..n-1
        if h[i-1] > h[i-2] and h[i-1] > h[i]:
            d = abs(i - ns2)
            if dmin is None or d < dmin or (d == dmin and h[i-1] > h[best-1]):
                dmin = d; best = i
    return best


def _window(h, ic):
    n = len(h); hc = h[ic-1]
    i = ic
    while True:
        if i == 1: iL = 1; break
        if h[i-2] < hc: i -= 1; continue
        iL = i - 1; break
    i = ic
    while True:
        if i == n: iR = n; break
        if h[i] < hc: i += 1; continue
        iR = i + 1; break
    return iL, iR


def _rms(xs, hs, xp, hp):
    return float(np.sqrt(np.mean((np.interp(xs, xp, hp) - hs) ** 2)))


def rnodes_v3(ridge, ns2=None, max_nodes_per_side=6):
    h = np.asarray(ridge, dtype=float); n = len(h); nsw = n - 1
    if ns2 is None: ns2 = nsw // 2
    H = lambda k: h[k-1]
    empty = dict(ok=False)

    ic = _central_peak(h, ns2)
    if ic < 0: return empty
    hc = H(ic)
    iL, iR = _window(h, ic)
    xs = np.arange(iL, iR + 1); hs = h[iL-1:iR]

    # ---- 3) first knees by rms (same as rms_wedge) ------------------------
    lefts = list(range(iL + 1, ic)) or [iL]
    rights = list(range(ic + 1, iR)) or [iR]
    best = None
    for i0 in lefts:
        for i1 in rights:
            xp = np.array([iL, i0, ic, i1, iR], float); hp = np.array([H(iL), H(i0), hc, H(i1), H(iR)])
            r = _rms(xs, hs, xp, hp)
            if best is None or r < best[0]: best = (r, i0, i1)
    _, x0, x1 = best

    # node lists (positions), left stored outer->inner, right inner->outer
    L = [iL, x0] if x0 != iL else [iL]
    R = [x1, iR] if x1 != iR else [iR]

    def polyline():
        xp = np.array(L + [ic] + R, float)
        hp = np.array([H(k) for k in L] + [hc] + [H(k) for k in R])
        return xp, hp

    # ---- 4) staircase refinement --------------------------------------------
    def staircase(side):
        """side='L' works outward from x_0 toward x_L; 'R' from x_1 toward x_R.
        Returns ordered node list (outer->inner for L, inner->outer for R) and foot."""
        if side == 'L':
            outer_lim, knee = iL, x0
        else:
            outer_lim, knee = iR, x1
        nodes = [knee]                 # inner-most first; we prepend outward
        inner_ref = ic                 # node on the crest side of the current knee
        cur = knee
        if cur == outer_lim:
            return [outer_lim], outer_lim
        while True:
            s_in = (hc - H(cur)) / (ic - cur) if inner_ref == ic else \
                   (H(inner_ref) - H(cur)) / (inner_ref - cur)
            s_out = (H(cur) - H(outer_lim)) / (cur - outer_lim)
            # orient slopes "toward the crest" so a ramp has s_out and s_in same sign
            if side == 'R':
                s_in, s_out = -s_in, -s_out
            if s_out == 0 or np.sign(s_out) != np.sign(s_in) or len(nodes) >= max_nodes_per_side:
                nodes.insert(0, outer_lim) if side == 'L' else nodes.append(outer_lim)
                return nodes, cur                                   # case A
            cands = list(range(min(outer_lim, cur) + 1, max(outer_lim, cur)))
            if not cands:
                nodes.insert(0, outer_lim) if side == 'L' else nodes.append(outer_lim)
                return nodes, outer_lim                             # ramp reaches the limit
            bestr = None
            for k in cands:
                if side == 'L':
                    Lt = [outer_lim, k] + nodes; Rt = R
                else:
                    Lt = L; Rt = nodes + [k, outer_lim]
                xp = np.array(Lt + [ic] + Rt, float)
                hp = np.array([H(q) for q in Lt] + [hc] + [H(q) for q in Rt])
                r = _rms(xs, hs, xp, hp)
                if bestr is None or r < bestr[0]: bestr = (r, k)
            k = bestr[1]
            if side == 'L': nodes.insert(0, k)
            else:           nodes.append(k)
            inner_ref, cur = cur, k

    L_nodes, foot_L = staircase('L')
    R_nodes, foot_R = staircase('R')
    # final polyline
    xw = np.array(L_nodes + [ic] + R_nodes, int)
    hw = np.array([H(k) for k in xw])
    rms = _rms(xs, hs, xw.astype(float), hw)

    disp_L = hc - H(foot_L); wdth_L = ic - foot_L
    disp_R = hc - H(foot_R); wdth_R = foot_R - ic
    return dict(ok=True, xwedge=xw, hwedge=hw, rms=rms,
                xnodes_L=np.array(L_nodes), hnodes_L=np.array([H(k) for k in L_nodes]),
                xnodes_R=np.array(R_nodes), hnodes_R=np.array([H(k) for k in R_nodes]),
                xc=ic, x0=x0, x1=x1, xL=iL, xR=iR, foot_L=foot_L, foot_R=foot_R,
                disp_L=disp_L, wdth_L=wdth_L, disp_R=disp_R, wdth_R=wdth_R,
                riseq=disp_L, fallq=disp_R, mxdis=max(disp_L, disp_R), hwdth=wdth_L + wdth_R,
                mxdis3=hc - min(H(x0), H(x1)), hwdth3=x1 - x0)


def plot_wedge(ax, ridge, w, title=None):
    """Profile (grey), staircase polyline (blue), feet (red) and height bars."""
    h = np.asarray(ridge, float); x = np.arange(1, len(h) + 1)
    ax.plot(x, h, color='0.5', lw=2, label='profile')
    if not w.get('ok'):
        ax.set_title((title or '') + '  (no wedge)'); return
    ax.plot(w['xwedge'], w['hwedge'], 'o-', color='tab:blue', lw=2, ms=5, label='staircase wedge')
    ax.plot([w['x0'], w['xc'], w['x1']], [h[w['x0']-1], h[w['xc']-1], h[w['x1']-1]], 's--',
            color='tab:orange', ms=6, lw=1, label='3-pt wedge')
    for foot, disp in ((w['foot_L'], w['disp_L']), (w['foot_R'], w['disp_R'])):
        ax.plot(foot, h[foot-1], 'o', color='red', ms=9, mfc='none', mew=2)
        ax.vlines(foot, h[foot-1], h[w['xc']-1], color='red', lw=1, ls=':')
        ax.text(foot, h[w['xc']-1], f"{disp:.0f} m", color='red', ha='center', va='bottom', fontsize=9)
    ax.axvline(w['xL'], color='green', lw=0.8, ls='--'); ax.axvline(w['xR'], color='green', lw=0.8, ls='--')
    ax.set_title((title or '') + f"  mxdis3={w['mxdis3']:.0f}  stair={w['mxdis']:.0f}  hwdth {w['hwdth3']}->{w['hwdth']}")
    ax.legend(fontsize=8)


if __name__ == '__main__':
    import sys
    sys.path.insert(0, __file__.rsplit('/', 1)[0] or '.')
    x = np.arange(33, dtype=float)
    # ramp case modelled on Julio's 2026-11-07 figure (nsw=32 here for the demo)
    ramp = np.interp(x, [0, 8, 12, 18, 25, 26, 27, 32], [-240, 0, 680, 320, 450, 500, 0, 0])
    ramp[:8] = np.linspace(-240, -5, 8); ramp[1:3] -= 10
    knob = np.array([-56.88,-44.72,-31.79,-19.91,-9.17,-14.90,-1.90,12.49,34.64,61.95,75.96,91.04,110.82,110.38,
                     108.44,101.81,92.94,103.74,123.47,130.13,122.91,127.45,142.43,150.16,172.37,191.10])
    for name, p in [('ramp (figure)', ramp[:27]), ('knob (31 Oct)', knob)]:
        w = rnodes_v3(p)
        print(f"{name}: xwedge={w['xwedge']}  h={np.round(w['hwedge'])}")
        print(f"   3-pt: mxdis={w['mxdis3']:.0f} hwdth={w['hwdth3']}   staircase: disp_L={w['disp_L']:.0f} "
              f"disp_R={w['disp_R']:.0f} mxdis={w['mxdis']:.0f} hwdth={w['hwdth']}  feet=({w['foot_L']},{w['foot_R']}) rms={w['rms']:.1f}")
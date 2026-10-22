"""Three-turning-point WKB (Bohr-Sommerfeld + Gamow tunnelling) for the trapped modes of
the horizonless Hayward metric, V_l = f [ l(l+1)/r^2 + f'/r ], units ell = 1.
Writes output/wkb_output.txt (everything printed, including the rows of Table 1) and output/wkb.json.

Gamow formula (corrected, see paper Sec. 5.4): with E = w^2 and
  T_s = 2 * int dr / (f sqrt(E-V))          (dimension length^2),
  T   = w_R * T_s                            (true round-trip time; group velocity sqrt(w^2-V)/w),
  Im(w) = -exp(-2*S_barrier) / (2 T).

Numerics: scipy.integrate.quad / scipy.optimize.brentq.  All three integrals have square-root endpoint
behaviour at the turning points; they are evaluated after the exact substitution
r = (a+b)/2 - (b-a)/2 cos(theta), theta in [0, pi], which makes the integrands smooth.
Independent cross-check of these integrals with mpmath (tanh-sinh, 30 digits): crosscheck_wkb_mpmath.py.
Usage: python3 wkb.py        (a few seconds)

PROVENANCE / CORRECTION NOTICE (added after v1.0.1; no code below is changed by this notice).
The author's ORIGINAL wkb.py (kept verbatim as wkb_original.py) used
    Im(w) = -(1/Ts) exp(-2 S),   Ts = 2 int dr/(f sqrt(E-V))
which is not dimensionally a rate (Ts has dimension length^2 when ell=1 is restored).  This file uses the round-trip
TIME T = w_R * Ts and Im(w) = -exp(-2S)/(2T), i.e. Im_new = Im_old/(2 w_R)  (verified numerically to 1e-7 on every
Table-1 row; Re(w) and S are unchanged).  Table 1 and Sec. 5.4 of the manuscript therefore differ from the first version
(l=1, M=0.90 M*: Q 41.5 -> 23.0; l=2: 92.2 -> 79.3; see ../PROVENANCE.md and Sec. 5.4 of the paper).  The quadrature substitution
and the analytic f' are numerically equivalent to the original (relative effect <= 1e-7 on Im(w), 1e-10 on V).
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from _util import tee_stdout, save_json

l_scale = 1.0
Mstar = 3*np.sqrt(3)/4*l_scale


def f_of(M):
    return lambda r: 1 - 2*M*r**2/(r**3+2*M*l_scale**2)


def fp_of(M):
    """analytic f'(r)"""
    def fp(r):
        D = r**3 + 2*M*l_scale**2
        return -2*M*(2*r*D - 3*r**4)/D**2
    return fp


def V_of(M, l):
    f = f_of(M); fp = fp_of(M)
    return lambda r: f(r)*(l*(l+1)/r**2 + fp(r)/r)


def well_barrier(M, l):
    V = V_of(M, l)
    r = np.linspace(1e-4, 40, 300000); Vv = V(r); dV = np.gradient(Vv, r)
    idx = np.where(dV[:-1]*dV[1:] < 0)[0]
    pts = [brentq(lambda x: (V(x+1e-5)-V(x-1e-5))/2e-5, r[i], r[i+1]) for i in idx]
    assert len(pts) == 2
    rmin, rmax = pts
    return rmin, V(rmin), rmax, V(rmax)


def turning_points(V, E, rmin, rmax, rhi=60.0):
    r1 = brentq(lambda x: V(x)-E, 1e-6, rmin)
    r2 = brentq(lambda x: V(x)-E, rmin, rmax)
    r3 = brentq(lambda x: V(x)-E, rmax, rhi)
    return r1, r2, r3


def theta_quad(g, a, b):
    """int_a^b g(r) dr with the substitution r = (a+b)/2 - (b-a)/2 cos(theta)."""
    m, h = 0.5*(a+b), 0.5*(b-a)
    val, _ = quad(lambda th: g(m - h*np.cos(th))*h*np.sin(th), 0.0, np.pi, limit=200, epsabs=1e-10, epsrel=1e-8)
    return val


def bohr_sommerfeld_residual(E, V, f, rmin, rmax, n):
    r1, r2, _ = turning_points(V, E, rmin, rmax)
    val = theta_quad(lambda r: np.sqrt(max(E-V(r), 0.0))/f(r), r1, r2)
    return val-(n+0.5)*np.pi


def solve_mode(M, l, n=0):
    f = f_of(M); V = V_of(M, l)
    rmin, Vmin, rmax, Vmax = well_barrier(M, l)
    lo, hi = Vmin*1.0001, Vmax*0.9999
    Elo = bohr_sommerfeld_residual(lo, V, f, rmin, rmax, n)
    Ehi = bohr_sommerfeld_residual(hi, V, f, rmin, rmax, n)
    if Elo*Ehi > 0:
        return None  # this (n,l) does not fit in the well
    E = brentq(lambda E: bohr_sommerfeld_residual(E, V, f, rmin, rmax, n), lo, hi, xtol=1e-10)
    wR = np.sqrt(E)
    r1, r2, r3 = turning_points(V, E, rmin, rmax)
    barrier = theta_quad(lambda r: np.sqrt(max(V(r)-E, 0.0))/f(r), r2, r3)
    Ts = 2*theta_quad(lambda r: 1.0/(f(r)*np.sqrt(max(E-V(r), 1e-300))), r1, r2)   # 2 int dr/(f sqrt(E-V)), E=w^2
    T = wR*Ts                    # true round-trip time in the well
    wI = -np.exp(-2*barrier)/(2.0*T)
    return dict(wR=wR, wI=wI, T=T, barrier_action=barrier, r1=r1, r2=r2, r3=r3, Vmin=Vmin, Vmax=Vmax, E=E)


def fmt_sci(x):
    m, e = f"{x:.1e}".split("e")
    return f"${m}\\times10^{{{int(e)}}}$"


TABLE1 = [(0.90, 1), (0.90, 2), (0.90, 4), (0.95, 1), (0.95, 2), (0.99, 1)]

if __name__ == "__main__":
    tee_stdout("wkb_output.txt")
    J = {"modes": {}, "gap_l5": {}, "table1_latex": []}
    print(f"{'M/M*':>6} {'l':>3} {'n':>2}  {'Re(w) l':>10} {'Im(w) l':>12} {'Q=Re/2|Im|':>12} {'tau=1/|Im|':>12} {'barrier action':>15}")
    for Mfrac in [0.90, 0.95, 0.99]:
        M = Mfrac*Mstar
        for l in [1, 2, 3, 4, 5, 10, 20, 30]:
            for n in [0, 1]:
                res = solve_mode(M, l, n)
                key = f"{Mfrac:.2f}_{l}_{n}"
                if res is None:
                    print(f"{Mfrac:6.2f} {l:3d} {n:2d}   (no mode fits: well too shallow for this n)")
                    J["modes"][key] = None
                    continue
                Q = res['wR']/(2*abs(res['wI'])); tau = 1.0/abs(res['wI'])
                print(f"{Mfrac:6.2f} {l:3d} {n:2d}  {res['wR']:10.5f} {res['wI']:12.3e} {Q:12.4e} {tau:12.4e} {res['barrier_action']:15.4f}")
                J["modes"][key] = dict(wR=res['wR'], wI=res['wI'], Q=Q, tau=tau, barrier_action=res['barrier_action'])

    print("\n--- well-to-barrier gap Vmax-Vmin at l=5 (units ell^-2) ---")
    for Mfrac in [0.90, 0.99]:
        rmin, Vmin, rmax, Vmax = well_barrier(Mfrac*Mstar, 5)
        print(f"M/M*={Mfrac:.2f}: {Vmax-Vmin:.3f}")
        J["gap_l5"][f"{Mfrac:.2f}"] = Vmax-Vmin

    print("\n--- Table 1 rows (LaTeX, n=0) ---")
    for Mfrac, l in TABLE1:
        res = solve_mode(Mfrac*Mstar, l, 0)
        Q = res['wR']/(2*abs(res['wI'])); tau = 1.0/abs(res['wI'])
        row = f"${Mfrac:.2f}$ & ${l}$ & ${res['wR']:.3f}$ & {fmt_sci(Q)} & {fmt_sci(tau)} \\\\"
        print(row)
        J["table1_latex"].append(row)
    save_json("wkb.json", J)

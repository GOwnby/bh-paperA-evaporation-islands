"""Direct (non-WKB) check of the trapped-mode frequencies: integrate
   f d/dr( f dpsi/dr ) + (w^2 - V_l) psi = 0           (G=c=1, ell=1)
for complex w.  The regular solution (psi ~ r^{l+1}) is integrated outward from r=1e-3 to a matching point r_m
in the well; the outgoing solution (psi ~ e^{i w x}) is started at the end of a rotated contour
r = R0 + s e^{i theta} (theta=0.35, s in [0,60]), where it decays, and integrated inward to r_m (the unwanted
incoming solution is suppressed by e^{-2 Re(w) Im(r)}); w is found by a secant iteration on the mismatch of
the logarithmic derivatives.  Standard complex-contour/shooting method (no WKB, no ansatz).
Writes output/qnm_shooting.txt and .json.  Usage: python3 qnm_shooting.py [REFINE]   (REFINE = step-size refinement factor, default 1; pure-Python RK4; see README for run time)
"""
import cmath
import numpy as np
from _util import tee_stdout, save_json
from wkb import solve_mode

Ms = 3*np.sqrt(3)/4


def mk(Mf, l):
    M = Mf*Ms
    f = lambda r: 1 - 2*M*r**2/(r**3+2*M)
    def fp(r):
        D = r**3+2*M
        return -2*M*(2*r*D-3*r**4)/D**2
    V = lambda r: f(r)*(l*(l+1)/r**2 + fp(r)/r)
    return f, fp, V


def make_step(Mf, l, w):
    """RK4 step for (psi, dpsi/dr) with plain complex scalars (fast)."""
    M = Mf*Ms; ll = l*(l+1); w2 = w*w
    def deriv(r, y0, y1):
        D = r**3 + 2*M
        fr = 1 - 2*M*r**2/D
        fpr = -2*M*(2*r*D - 3*r**4)/D**2
        V = fr*(ll/r**2 + fpr/r)
        return y1, -(fpr/fr)*y1 - (w2 - V)/fr**2*y0
    def step(r, y0, y1, dr):
        a0, a1 = deriv(r, y0, y1)
        b0, b1 = deriv(r+dr/2, y0+dr/2*a0, y1+dr/2*a1)
        c0, c1 = deriv(r+dr/2, y0+dr/2*b0, y1+dr/2*b1)
        d0, d1 = deriv(r+dr, y0+dr*c0, y1+dr*c1)
        return y0 + dr/6*(a0+2*b0+2*c0+d0), y1 + dr/6*(a1+2*b1+2*c1+d1)
    return step


def mismatch(w, Mf, l, rm=2.0, R0=7.0, S=60.0, th=0.35, h=0.01, refine=1):
    h = h/refine
    step = make_step(Mf, l, w)
    M = Mf*Ms
    r = 1e-3; y0 = complex(r**(l+1)); y1 = complex((l+1)*r**l)
    n = int((rm-r)/(1e-3/refine)); dr = (rm-r)/n
    for i in range(n):
        y0, y1 = step(r, y0, y1, dr); r += dr
        a = abs(y0)
        if a > 1e50: y0, y1 = y0/a, y1/a
    reg = y1/y0
    e = cmath.exp(1j*th); rend = R0+S*e
    fend = 1 - 2*M*rend**2/(rend**3+2*M)
    o0, o1 = 1.0+0j, 1j*w/fend; r = rend
    ns = int(S/h); dr = -h*e
    for i in range(ns):
        o0, o1 = step(r, o0, o1, dr); r += dr
        a = abs(o0); o0, o1 = o0/a, o1/a
    ns = int((R0-rm)/(0.002/refine)); dr = -(R0-rm)/ns; r = complex(R0)
    for i in range(ns):
        o0, o1 = step(r, o0, o1, dr); r += dr
        a = abs(o0); o0, o1 = o0/a, o1/a
    return reg - o1/o0


def solve(w0, Mf, l, **kw):
    w1 = w0*(1+1e-3) + 1e-4j
    f0 = mismatch(w0, Mf, l, **kw); f1 = mismatch(w1, Mf, l, **kw)
    for it in range(12):
        w2 = w1 - f1*(w1-w0)/(f1-f0)
        w0, f0 = w1, f1; w1 = w2; f1 = mismatch(w1, Mf, l, **kw)
        if abs(w1-w0) < 1e-9: break
    return w1


if __name__ == "__main__":
    import sys
    REF = int(sys.argv[1]) if len(sys.argv) > 1 else 1      # step-size refinement factor (verify.sh --full uses 2)
    tee_stdout("qnm_shooting.txt" if REF == 1 else f"qnm_shooting_refine{REF}.txt")
    if REF != 1:
        print(f"# all RK4 step sizes divided by {REF}")
    J = {"direct": {}, "robustness_rm1.8": {}}
    print("M/M*  l | direct: omega            Q=Re/2|Im|   | WKB (wkb.py): Re, Im      Q       | |Im|_direct/|Im|_WKB")
    for Mf, l in [(0.90, 1), (0.90, 2), (0.90, 3), (0.90, 4), (0.95, 1), (0.95, 2), (0.99, 1)]:
        wk = solve_mode(Mf*Ms, l, 0)
        w = solve(complex(wk['wR'], 0.5*wk['wI']), Mf, l, refine=REF)
        print(f"{Mf:.2f} {l:2d} | {w.real:.5f}{w.imag:+.3e}i  {w.real/(2*abs(w.imag)):.3e} | {wk['wR']:.5f} {wk['wI']:+.3e}  {wk['wR']/(2*abs(wk['wI'])):.3e} | {abs(w.imag)/abs(wk['wI']):.2f}", flush=True)
        J["direct"][f"{Mf:.2f}_{l}"] = dict(re=w.real, im=w.imag, Q=w.real/(2*abs(w.imag)), ratio_to_wkb=abs(w.imag)/abs(wk['wI']))
    print("\nrobustness: matching point r_m = 1.8 instead of 2.0 (M/M*=0.90, l=2; 0.95, l=1)")
    for Mf, l in [(0.90, 2), (0.95, 1)]:
        wk = solve_mode(Mf*Ms, l, 0)
        w = solve(complex(wk['wR'], 0.5*wk['wI']), Mf, l, rm=1.8, refine=REF)
        print(f"{Mf:.2f} {l:2d} | {w.real:.5f}{w.imag:+.3e}i", flush=True)
        J["robustness_rm1.8"][f"{Mf:.2f}_{l}"] = dict(re=w.real, im=w.imag)
    save_json("qnm_shooting.json" if REF == 1 else f"qnm_shooting_refine{REF}.json", J)

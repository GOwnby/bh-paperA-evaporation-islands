"""Independent check of the scipy quadratures in wkb.py.  For the six rows of Table 1 (n=0) the turning points
and E=w_R^2 found by wkb.py (scipy) are polished with mpmath.findroot, and the Bohr-Sommerfeld action, the
barrier action and the round-trip time are recomputed with mpmath tanh-sinh quadrature (30 digits) directly
on the singular endpoint integrands (no variable substitution).  Reports the Bohr-Sommerfeld residual and the
relative differences of Re(w), Im(w), Q, tau between the two implementations.
Writes output/crosscheck_wkb_mpmath.txt and .json.   Usage: python3 crosscheck_wkb_mpmath.py  (about 1 min; needs mpmath)
"""
import numpy as np
import mpmath as mp
from _util import tee_stdout, save_json
from wkb import solve_mode, Mstar, TABLE1

mp.mp.dps = 30
tee_stdout("crosscheck_wkb_mpmath.txt")


def run(Mfrac, l):
    M = mp.mpf(Mfrac)*3*mp.sqrt(3)/4
    D = lambda r: r**3 + 2*M
    f = lambda r: 1 - 2*M*r**2/D(r)
    fp = lambda r: -2*M*(2*r*D(r) - 3*r**4)/D(r)**2
    V = lambda r: f(r)*(l*(l+1)/r**2 + fp(r)/r)
    sc = solve_mode(float(M), l, 0)
    E = mp.mpf(sc['E'])
    r1, r2, r3 = [mp.findroot(lambda x: V(x)-E, mp.mpf(sc[k])) for k in ('r1', 'r2', 'r3')]
    act = mp.quad(lambda r: mp.sqrt(max(E-V(r), 0))/f(r), [r1, r2])          # Bohr-Sommerfeld integral
    resid = act - mp.pi/2
    S = mp.quad(lambda r: mp.sqrt(max(V(r)-E, 0))/f(r), [r2, r3])
    Ts = 2*mp.quad(lambda r: 1/(f(r)*mp.sqrt(E-V(r))), [r1, r2])
    wR = mp.sqrt(E); T = wR*Ts
    wI = -mp.exp(-2*S)/(2*T)
    Q = wR/(2*abs(wI)); tau = 1/abs(wI)
    scQ = sc['wR']/(2*abs(sc['wI'])); sctau = 1/abs(sc['wI'])
    rel = lambda a, b: float(abs(a-b)/abs(b))
    return dict(M_over_Mstar=Mfrac, l=l, BS_residual=float(resid), wR_mp=float(wR), wI_mp=float(wI), Q_mp=float(Q), tau_mp=float(tau),
                rel_wI=rel(sc['wI'], wI), rel_Q=rel(scQ, Q), rel_tau=rel(sctau, tau), rel_barrier=rel(sc['barrier_action'], S))


rows = []
print("M/M*  l | BS residual (should be ~0; E solved to xtol 1e-10) | rel diff scipy vs mpmath: Im(w)   Q   tau   barrier action")
for Mfrac, l in TABLE1:
    r = run(Mfrac, l); rows.append(r)
    print(f"{Mfrac:.2f} {l:2d} | {r['BS_residual']:+.2e} | {r['rel_wI']:.2e} {r['rel_Q']:.2e} {r['rel_tau']:.2e} {r['rel_barrier']:.2e}", flush=True)
save_json("crosscheck_wkb_mpmath.json", {"rows": rows})

"""Added after the first archived version (v1.0.1).  Justification of the corrected Gamow prefactor in wkb.py.

(1) wkb_original.py is the author's original script (verbatim; its functions are loaded here by executing the part of the file
    before its print loops).  On the six rows of Table 1 it reproduces the Table-1 values of v1.0.1 (v101 rows below), and
    Im(w)_corrected * 2 w_R = Im(w)_original to ~1e-7 (the original lacks the factor 1/(2 w_R): its Ts = 2 int dr/(f sqrt(E-V)) has
    dimension length^2, the round-trip time is T = w_R * Ts).  Re(w) and the barrier action are unchanged.
(2) Evidence that the corrected prefactor is the right one: the direct (non-WKB) complex-frequency integration of qnm_shooting.py,
    for M = 0.90 M*, l = 1..7, divided by |Im w| of each WKB version.  The ratio to the corrected WKB stays within ~1.1-1.35, the
    ratio to the original WKB drifts by a factor ~5 with l (it scales like 1/(2 w_R)).
Writes output/compare_wkb_original.txt and .json.   Usage: python3 compare_wkb_original.py   (about 40 s)
"""
import os, warnings
import numpy as np
from _util import HERE, tee_stdout, save_json
import wkb as new
import qnm_shooting as q

warnings.filterwarnings("ignore")      # the original emits an IntegrationWarning (roundoff in quad), see PROVENANCE.md
tee_stdout("compare_wkb_original.txt")

src = open(os.path.join(HERE, "wkb_original.py")).read().split("print(f\"{'M/M*'")[0]
old = {}
exec(src, old)

V101_TABLE = [   # rows of Table 1 of v1.0.1 (Re(w), Q, tau as printed in that paper)
    (0.90, 1, "0.277", "4.2e+01", "3.0e+02"), (0.90, 2, "0.430", "9.2e+01", "4.3e+02"), (0.90, 4, "0.732", "1.2e+03", "3.1e+03"),
    (0.95, 1, "0.200", "8.9e+02", "8.9e+03"), (0.95, 2, "0.309", "1.7e+04", "1.1e+05"), (0.99, 1, "0.090", "1.3e+06", "3.0e+07")]
J = {"rows": {}}
print("(1) original wkb.py functions vs corrected wkb.py, Table-1 rows (n=0)")
print("M/M*  l | Re(w) | original: Im, Q, tau | corrected: Im, Q, tau | Im_new*2wR/Im_old - 1 | v1.0.1 table row reproduced")
ok_all = True; dmax = 0.0; dre = 0.0; ds = 0.0
for Mf, l, sre, sq, st in V101_TABLE:
    o = old['solve_mode'](Mf*old['Mstar'], l, 0); n = new.solve_mode(Mf*new.Mstar, l, 0)
    Qo = o['wR']/(2*abs(o['wI'])); to = 1/abs(o['wI']); Qn = n['wR']/(2*abs(n['wI'])); tn = 1/abs(n['wI'])
    d = abs(n['wI']*2*o['wR']/o['wI'] - 1)
    dmax = max(dmax, d); dre = max(dre, abs(n['wR']/o['wR']-1)); ds = max(ds, abs(n['barrier_action']/o['barrier_action']-1))
    ok = (f"{o['wR']:.3f}" == sre and f"{Qo:.1e}" == sq and f"{to:.1e}" == st)
    ok_all = ok_all and ok
    print(f"{Mf:.2f} {l} | {o['wR']:.5f} | {o['wI']:.4e} {Qo:.3e} {to:.3e} | {n['wI']:.4e} {Qn:.3e} {tn:.3e} | {d:.1e} | {ok}")
    J["rows"][f"{Mf:.2f}_{l}"] = dict(Im_old=o['wI'], Im_new=n['wI'], Q_old=Qo, Q_new=Qn, tau_old=to, tau_new=tn, rel_prefactor=d, v101_row_ok=ok)
J["v101_table_reproduced"] = {"ok": bool(ok_all)}
J["max_rel_Im_prefactor"] = dmax; J["max_rel_wR"] = dre; J["max_rel_barrier"] = ds
print(f"v1.0.1 Table 1 reproduced by the original functions: {ok_all}; max |Im_new*2wR/Im_old-1| = {dmax:.2e}; max rel. diff Re(w) {dre:.1e}, S {ds:.1e}")

print("\n(2) direct complex-frequency integration (qnm_shooting.py) vs both WKB versions, M = 0.90 M*")
print(" l | direct omega | |Im|_direct/|Im|_corrected | |Im|_direct/|Im|_original")
rn, ro = [], []
J["direct"] = {}
for l in range(1, 8):
    Mf = 0.90
    wk = new.solve_mode(Mf*new.Mstar, l, 0); wo = old['solve_mode'](Mf*old['Mstar'], l, 0)
    w = q.solve(complex(wk['wR'], 0.5*wk['wI']), Mf, l)
    a, b = abs(w.imag)/abs(wk['wI']), abs(w.imag)/abs(wo['wI'])
    rn.append(a); ro.append(b)
    print(f"{l:2d} | {w.real:.5f}{w.imag:+.4e}i | {a:.3f} | {b:.3f}", flush=True)
    J["direct"][str(l)] = dict(re=w.real, im=w.imag, ratio_corrected=a, ratio_original=b)
J["spread_ratio_corrected"] = max(rn)/min(rn); J["spread_ratio_original"] = max(ro)/min(ro)
J["min_ratio_corrected"] = min(rn); J["max_ratio_corrected"] = max(rn)
print(f"max/min over l: corrected {J['spread_ratio_corrected']:.2f}, original {J['spread_ratio_original']:.2f}")
save_json("compare_wkb_original.json", J)

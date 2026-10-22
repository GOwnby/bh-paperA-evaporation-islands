"""Numerical (mpmath/numpy) checks of the analytic statements and quoted numbers of Sections 2-3, 4.4, 5.1-5.3 of the paper.
The symbolic (exact, sympy) counterparts are in check_symbolic.py.  Does not need SciPy.
Writes output/analytic_checks.txt and output/analytic_checks.json (read by ../verify.py).
Usage: python3 check_analytic.py   (a few seconds; needs numpy, mpmath)
"""
import numpy as np, mpmath as mp
from _util import tee_stdout, save_json
tee_stdout("analytic_checks.txt", with_scipy=False)
print(f"# mpmath {mp.__version__}")
mp.mp.dps = 30
J = {}
# CODATA 2018 SI values
hbar, c, G, Msun, kB = 1.054571817e-34, 2.99792458e8, 6.67430e-11, 1.98847e30, 1.380649e-23
yr = 365.25*86400.0
mP = np.sqrt(hbar*c/G); lP = np.sqrt(hbar*G/c**3)
Ms = 3*mp.sqrt(3)/4                       # M*/ell, G=c=1
f = lambda r, M: 1 - 2*M*r**2/(r**3+2*M)  # ell = 1


def hor_roots(M, extra=60):
    return sorted([float(mp.re(x)) for x in mp.polyroots([1, -2*M, 0, 2*M], maxsteps=300, extraprec=extra)
                   if abs(mp.im(x)) < 1e-12 and mp.re(x) > 0])


print("== Sec. 2.2: critical mass and fold scalings (ell=1) ==")
Mstar_kg = float(Ms)*lP*c**2/G
print(f"M*(ell=l_P) = {Mstar_kg:.3e} kg = {Mstar_kg/mP:.2f} m_P   (m_P = {mP:.4e} kg)")
J["Mstar_kg"] = Mstar_kg; J["Mstar_over_mP"] = Mstar_kg/mP; J["mP_kg"] = mP
print("(r_pm - r_h)/[2 sqrt2 3^(1/4)/3 sqrt(eps)]  (expect -1, +1 as eps->0):")
J["fold"] = {}
for eps in [1e-6, 1e-8]:
    M = Ms+eps
    rs = hor_roots(M)
    pred = 2*mp.sqrt(2)*mp.mpf(3)**0.25/3*mp.sqrt(eps)
    ratios = [float((r-mp.sqrt(3))/pred) for r in (rs[0], rs[-1])]
    print(f"  eps={eps:g}: {ratios}")
    k = mp.diff(lambda r: f(r, M), rs[-1])/2
    kr = float(k/((2*mp.sqrt(2)*mp.mpf(3)**0.25/9)*mp.sqrt(eps)))
    print(f"  kappa/(c_kappa sqrt eps) = {kr:.5f}")
    J["fold"][f"{eps:g}"] = dict(r_minus_ratio=ratios[0], r_plus_ratio=ratios[1], kappa_ratio=kr)
ck = 2*mp.sqrt(2)*mp.mpf(3)**0.25/9
print(f"c_kappa = {float(ck):.4f}")
J["c_kappa"] = float(ck)
fpp = float(mp.diff(lambda r: f(r, Ms), mp.sqrt(3), 2)/2*3)
print(f"(1/2) f''(r_h) * 3 ell^2 at M* = {fpp:.6f}  (expect 1)")
J["half_fpp_times_3"] = fpp
print("eta expansions: [exact, series]")
J["eta"] = {}
for M in [1e3, 1e4]:
    eta = 1/M**2
    rp = hor_roots(M, 80)[-1]
    k = mp.diff(lambda r: f(r, M), rp)/2
    a = [rp/(2*M)-1, -eta/4-eta**2/8]; b = [float(k*4*M-1), -eta/2-3*eta**2/8]
    print(f"  M={M:g}: r+/2M-1 = [{a[0]:.6e}, {a[1]:.6e}]   kappa/kappaS-1 = [{b[0]:.6e}, {b[1]:.6e}]")
    J["eta"][f"{M:g}"] = dict(r_exact=a[0], r_series=a[1], kappa_exact=b[0], kappa_series=b[1])
eta10 = (lP/(G*10*Msun/c**2))**2
print(f"eta(10 Msun, ell=l_P) = {eta10:.2e}")
J["eta_10Msun"] = eta10

print("\n== Sec. 3: Hawking law and adiabaticity ==")
P = lambda M: hbar*c**6/(15360*np.pi*G**2*M**2)
A = lambda M: hbar*c/(7680*np.pi*G*M**2)
TH = lambda M: hbar*c**3/(8*np.pi*G*M*kB)
tev = lambda M: 5120*np.pi*G**2*M**3/(hbar*c**4)   # from dM/dt = -P
M10 = 10*Msun
print(f"10 Msun: P = {P(M10):.2e} W, T_H = {TH(M10):.2e} K, A = {A(M10):.2e}, r_s/c = {2*G*M10/c**3:.2e} s")
print(f"evaporation time (Hawking law, no greybody) = {tev(M10)/yr:.2e} yr")
chk = P(M10)*(2*G*M10/c**3)/(M10*c**2)
print(f"check A = P (r_s/c)/(M c^2): {chk:.3e}")
J.update(P_10Msun_W=P(M10), TH_10Msun_K=TH(M10), A_10Msun=A(M10), rs_over_c_10Msun_s=2*G*M10/c**3, tevap_10Msun_yr=tev(M10)/yr,
         A_from_definition_10Msun=chk)
Mb = np.sqrt(hbar*c/(7680*np.pi*G))
print(f"M_break = {Mb:.2e} kg = {Mb/mP:.4f} m_P; M*/M_break = {Mstar_kg/Mb:.0f}; A(M*) = {A(Mstar_kg):.2e}")
J.update(M_break_kg=Mb, M_break_over_mP=Mb/mP, Mstar_over_Mbreak=Mstar_kg/Mb, A_Mstar=A(Mstar_kg), A_at_Mbreak=A(Mb))
rp_dot = 2*G*P(M10)/c**4; kap = c**3/(4*G*M10); rp = 2*G*M10/c**2   # |dr+/du|, kappa (1/s), r+ (SI)
print(f"far-field identity (SI, 10 Msun): delta/r+ = |r+'|/(kappa r+) = {rp_dot/kap/rp:.3e},  2A = {2*A(M10):.3e}")
J.update(far_field_delta_over_rp=rp_dot/kap/rp, far_field_2A=2*A(M10))
eps_c = (float(ck)/(2*np.pi))**2
print(f"eps_c = (c_kappa/2pi)^2 = {eps_c:.2e} Planck masses (ell = l_P)")
J["eps_c_mP"] = eps_c
Sbits = 4*np.pi*G*M10**2/(hbar*c)/np.log(2)
print(f"S_BH(10 Msun) = {Sbits:.2e} bits")
J["S_BH_bits_10Msun"] = Sbits
E = Mstar_kg*c**2
print(f"E = M* c^2 = {E:.2e} J ; m_P c^2 = {mP*c**2:.2e} J")
bek = 2*np.pi*E*lP/(hbar*c)/np.log(2)    # Bekenstein bound S <= 2 pi k R E/(hbar c), R = l_P
print(f"Bekenstein bound for E = M* c^2 confined to R = l_P: {bek:.1f} bits")
J.update(E_Mstar_J=E, mP_c2_J=mP*c**2, bekenstein_bits_Mstar_in_lP=bek)

print("\n== Sec. 5.1: horizonless geometry (ell=1) ==")
J["rmin"] = {}
for Mf in [0.5, 0.9, 0.99]:
    M = Mf*float(Ms); r = np.linspace(0.05, 6, 600001); fv = f(r, M); i = fv.argmin()
    pr, pf = 2**(2/3)*M**(1/3), 1-Mf**(2/3)
    print(f"M/M*={Mf}: r_min numeric {r[i]:.4f} vs 2^(2/3) M^(1/3) = {pr:.4f};  f_min {fv[i]:.5f} vs 1-(M/M*)^(2/3) = {pf:.5f}")
    J["rmin"][str(Mf)] = dict(r_num=r[i], r_formula=pr, f_num=fv[i], f_formula=pf)

print("\n== Sec. 5.3: stable inner light ring = local minimum of f/r^2 ==")
def stat(M):
    r = np.linspace(0.05, 30, 600000); g = f(r, M)/r**2; d = np.diff(g)
    s = np.where(d[:-1]*d[1:] < 0)[0]
    return [(float(r[i+1]), float(g[i+1]), 'min' if d[i] < 0 else 'max') for i in s]
def has_min(mf):
    M = mf*float(Ms); r = np.linspace(1.5, 3.5, 400001); g = f(r, M)/r**2; d = np.diff(g); return bool(np.any(d[:-1]*d[1:] < 0))
lo, hi = 0.80, 0.82
for _ in range(40):
    m = (lo+hi)/2
    if has_min(m): hi = m
    else: lo = m
print(f"inner light ring exists for M/M* >= {hi:.4f}")
J["lightring_threshold"] = hi
J["lightring"] = {}
for mf in [0.90, 0.95, 0.99, 1.0]:
    st = stat(mf*float(Ms))
    print(f"M/M*={mf}: " + ", ".join(f"r={a:.4f}, f/r^2={b:.5f} ({k})" for a, b, k in st))
    J["lightring"][f"{mf}"] = [dict(r=a, g=b, kind=k) for a, b, k in st]
save_json("analytic_checks.json", J)

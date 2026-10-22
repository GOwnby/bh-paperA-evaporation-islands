"""Exact (sympy) verification of the analytic statements of the paper: critical mass and horizon structure,
fold scalings, c_kappa, eta expansions, time-dilation closed form, Hawking-law identities, near-extremal scaling
chain, horizonless geometry (r_min, f(r_min), no horizon for M<M*), light-ring condition and its threshold mass,
and the reduction of the scalar wave equation to the Regge-Wheeler form of Eq. (RWpotential).
Each check prints PASS/FAIL with the derived expression; results go to output/symbolic_checks.txt and .json.
Does not need SciPy.  Usage: python3 check_symbolic.py   (about 1 min; needs sympy)
"""
import traceback
import sympy as sp
from _util import tee_stdout, save_json

tee_stdout("symbolic_checks.txt", with_scipy=False)
print(f"# sympy {sp.__version__}")
RES = {}


def check(name, fn):
    """run fn() -> (bool, detail str); record and print"""
    try:
        ok, detail = fn()
    except Exception as e:                                           # a failing derivation must be visible, not fatal
        ok, detail = False, "EXCEPTION " + repr(e) + "\n" + traceback.format_exc(limit=3)
    RES[name] = dict(ok=bool(ok), detail=str(detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


r, x, s, eta = sp.symbols('r x s eta', positive=True)
M, L, eps = sp.symbols('M ell epsilon', positive=True)
f = lambda rr, MM, LL=L: 1 - 2*MM*rr**2/(rr**3 + 2*MM*LL**2)
P = lambda rr, MM, LL=L: rr**3 - 2*MM*rr**2 + 2*MM*LL**2          # horizon polynomial, f = P/D
D = lambda rr, MM, LL=L: rr**3 + 2*MM*LL**2
Mstar = 3*sp.sqrt(3)/4*L
rh = sp.sqrt(3)*L
ck = 2*sp.sqrt(2)*sp.Integer(3)**sp.Rational(1, 4)/9

# ---------------------------------------------------------------- Sec 2.2
def c_horizon():
    ok1 = sp.simplify(f(r, M) - P(r, M)/D(r, M)) == 0
    ok2 = sp.simplify(P(rh, Mstar)) == 0
    ok3 = sp.simplify(sp.diff(P(r, M), r).subs({r: rh, M: Mstar})) == 0
    P2 = sp.simplify(sp.diff(P(r, M), r, 2).subs({r: rh, M: Mstar}))
    PM = sp.simplify(sp.diff(P(r, M), M).subs({r: rh, M: Mstar}))
    ok4 = (P2 == 3*sp.sqrt(3)*L) and (PM == -4*L**2)
    return ok1 and ok2 and ok3 and ok4, f"f=P/D; P(r_h)=P'(r_h)=0 at M*; P''(r_h)={P2} != 0, dP/dM={PM} != 0"
check("horizon_polynomial_double_root_at_Mstar", c_horizon)


def c_disc():
    d = sp.factor(sp.discriminant(P(r, M), r))
    sols = [sol for sol in sp.solve(d, M) if sol.is_positive]
    ok = len(sols) == 1 and sp.simplify(sols[0] - Mstar) == 0
    sgn_hi = sp.sign(d.subs({M: sp.Rational(11, 10)*Mstar, L: 1}).evalf())
    sgn_lo = sp.sign(d.subs({M: sp.Rational(9, 10)*Mstar, L: 1}).evalf())
    ok = ok and sgn_hi > 0 and sgn_lo < 0
    return ok, f"discriminant_r(P) = {d}; its only positive root is M* = {sols}; sign>0 at 1.1M* (3 real roots), <0 at 0.9M* (1 real root)"
check("discriminant_unique_critical_mass", c_disc)


def c_nohorizon():
    Pm = sp.simplify(P(4*M/3, M))          # P'(r) = r(3r-4M): global minimum of P on r>0 is at r=4M/3
    ok_der = sp.simplify(sp.diff(P(r, M), r) - r*(3*r - 4*M)) == 0
    cond = sp.simplify(Pm - 2*M*(L**2 - 16*M**2/27)) == 0
    sol = sp.solve(sp.Eq(L**2, 16*M**2/27), M)
    ok = ok_der and cond and sp.simplify(sol[0] - Mstar) == 0
    return ok, f"P'(r)=r(3r-4M); min_(r>0) P = P(4M/3) = {sp.factor(Pm)} > 0 iff M < 3 sqrt3 ell/4 = M*: no horizon for M<M*, two for M>M*"
check("no_horizon_below_Mstar_exact", c_nohorizon)

check("horizon_area_at_Mstar", lambda: (sp.simplify(4*sp.pi*rh**2 - 12*sp.pi*L**2) == 0, "A(M*) = 4 pi r_h^2 = 12 pi ell^2"))

a1, a2, a3 = sp.symbols('a1 a2 a3')


def c_fold():
    one = sp.Integer(1)   # units ell = 1 (all lengths scale with ell)
    Ps = sp.expand(P(rh.subs(L, 1) + a1*s + a2*s**2 + a3*s**3, (Mstar + s**2).subs(L, 1), one))
    c2 = sp.simplify(Ps.coeff(s, 2)); c3 = sp.simplify(Ps.coeff(s, 3))
    sols1 = sp.solve(c2, a1)
    pred = 2*sp.sqrt(2)*sp.Integer(3)**sp.Rational(1, 4)/3
    ok = (len(sols1) == 2) and all(sp.simplify(q**2 - pred**2) == 0 for q in sols1)
    # next order: a2 from the s^3 coefficient
    a2s = [sp.simplify(sp.solve(c3.subs(a1, q), a2)[0]) for q in sols1]
    return ok, f"eps=s^2: r_pm - r_h = a1 s + a2 s^2: a1 = {sols1} (= +-2 sqrt2 3^(1/4)/3 sqrt(ell eps)), a2 = {a2s} (same for both branches)"
check("fold_scaling_r_pm", c_fold)


def c_kappa():
    one = sp.Integer(1)
    q = (2*sp.sqrt(2)*sp.Integer(3)**sp.Rational(1, 4)/3)
    rp = rh.subs(L, 1) + q*s                       # r_+ to leading order; kappa = f'(r_+)/2 = P'(r_+)/(2 D(r_+)) on the root
    ex = sp.diff(P(r, Mstar.subs(L, 1), one), r)/(2*D(r, Mstar.subs(L, 1), one))
    kap = sp.series(ex.subs(r, rp), s, 0, 2).removeO()
    lead = sp.simplify(kap.coeff(s, 1))
    ok = sp.simplify(lead - ck) == 0
    f2 = sp.simplify(sp.diff(f(r, Mstar, one), r, 2).subs(r, rh.subs(L, 1)))
    return ok, f"kappa = f'(r_+)/2 = {lead} sqrt(eps) = c_kappa sqrt(eps), c_kappa = 2 sqrt2 3^(1/4)/9 = {sp.N(ck, 6)}"
check("c_kappa_value", c_kappa)


def c_ads2():
    one = sp.Integer(1)
    ser = sp.series(f(rh.subs(L, 1) + x, Mstar.subs(L, 1), one), x, 0, 4).removeO()
    c0, c1, c2c = [sp.simplify(ser.coeff(x, k)) for k in range(3)]
    ok = (c0 == 0) and (c1 == 0) and (sp.simplify(c2c - sp.Rational(1, 3)) == 0)
    return ok, f"at M*: f = {sp.simplify(ser)} + ..., i.e. f ~ (r-r_h)^2/(3 ell^2); AdS2 radius^2 = 2/f'' = 3 ell^2 = r_h^2"
check("near_horizon_AdS2_at_Mstar", c_ads2)


def c_eta():
    b1, b2, b3 = sp.symbols('b1 b2 b3')
    rho = 2*(1 + b1*eta + b2*eta**2 + b3*eta**3)            # rho = r_+/M, P/M^3 = rho^3 - 2 rho^2 + 2 eta  (ell^2 = eta M^2)
    eq = sp.expand(rho**3 - 2*rho**2 + 2*eta)
    sol1 = sp.solve(eq.coeff(eta, 1), b1)[0]
    sol2 = sp.solve(eq.coeff(eta, 2).subs(b1, sol1), b2)[0]
    sol3 = sp.solve(eq.coeff(eta, 3).subs({b1: sol1, b2: sol2}), b3)[0]
    rho_s = rho.subs({b1: sol1, b2: sol2, b3: sol3})
    kk = sp.series(2*(3*rho_s**2 - 4*rho_s)/(rho_s**3 + 2*eta), eta, 0, 4).removeO()   # kappa/kappa_Schw = 4M f'(r_+)/2
    ok = (sol1 == -sp.Rational(1, 4)) and (sol2 == -sp.Rational(1, 8)) and sp.simplify(kk - (1 - eta/2 - 3*eta**2/8 + kk.coeff(eta, 3)*eta**3)) == 0
    return ok, f"r_+ = 2M(1 + ({sol1}) eta + ({sol2}) eta^2 + ({sol3}) eta^3 ...); kappa/kappa_Schw = {sp.expand(kk)}"
check("eta_expansion_r_plus_and_kappa", c_eta)

check("time_dilation_closed_form", lambda: (
    sp.simplify(f(2*M, M) - L**2/(4*M**2 + L**2)) == 0,
    "f(r=2M) = ell^2/(4M^2+ell^2), hence Gamma = sqrt(f(2M)) = ell/sqrt(4M^2+ell^2) (Fig. 1 curve)"))

# ---------------------------------------------------------------- Sec 3
hb, c_, G_, kB, Mm = sp.symbols('hbar c G k_B M_m', positive=True)
Ph = hb*c_**6/(15360*sp.pi*G_**2*Mm**2)
TH = hb*c_**3/(8*sp.pi*G_*Mm*kB)
sigma = sp.pi**2*kB**4/(60*hb**3*c_**2)
rs = 2*G_*Mm/c_**2


def c_hawking():
    ok1 = sp.simplify(sigma*4*sp.pi*rs**2*TH**4 - Ph) == 0
    Acal = Ph*(rs/c_)/(Mm*c_**2)
    ok2 = sp.simplify(Acal - hb*c_/(7680*sp.pi*G_*Mm**2)) == 0
    Mb = sp.solve(sp.Eq(hb*c_/(7680*sp.pi*G_*Mm**2), 1), Mm)
    ok3 = len(Mb) == 1 and sp.simplify(Mb[0] - sp.sqrt(hb*c_/(7680*sp.pi*G_))) == 0
    return ok1 and ok2 and ok3, f"P = sigma (4 pi r_s^2) T_H^4 [Stefan-Boltzmann]; A = P (r_s/c)/(M c^2) = hbar c/(7680 pi G M^2); A=1 at M_break = {Mb[0]}"
check("hawking_power_adiabaticity_Mbreak", c_hawking)


def c_cuberoot():
    t, M0 = sp.symbols('t M0', positive=True)
    k = hb*c_**4/(15360*sp.pi*G_**2)                          # dM/dt = -P/c^2 = -k/M^2
    Mt = (M0**3 - 3*k*t)**sp.Rational(1, 3)
    ok1 = sp.simplify(Mt**2*sp.diff(Mt, t) + k) == 0 and sp.simplify(Mt.subs(t, 0) - M0) == 0
    tev = sp.solve(sp.Eq(M0**3 - 3*k*t, 0), t)[0]
    ok2 = sp.simplify(tev - 5120*sp.pi*G_**2*M0**3/(hb*c_**4)) == 0
    ok3 = sp.simplify(Ph/c_**2 - k/Mm**2) == 0
    return ok1 and ok2 and ok3, f"M(t) = (M0^3 - 3 k t)^(1/3), k = hbar c^4/(15360 pi G^2) [i.e. M ~ (t_f - t)^(1/3)]; t_evap = {tev}"
check("hawking_cube_root_law", c_cuberoot)


def c_scaling_chain():
    u_, K_, C_ = sp.symbols('u K C', positive=True)
    e = sp.Function('epsilon')
    sol = sp.dsolve(sp.Eq(e(u_).diff(u_), -K_*e(u_)**2), e(u_))
    en = sp.simplify(sol.rhs)                                # epsilon(u) = 1/(K u + C1)
    C1 = sp.Symbol('C1')
    ok1 = sp.simplify(en.subs(C1, C_) - 1/(K_*u_ + C_)) == 0
    # fold law: r_+ = r_h + q sqrt(epsilon)  (ell=1), kappa = ck sqrt(epsilon);  d eps/du = -K eps^2
    ee = sp.Symbol('e', positive=True); q = 2*sp.sqrt(2)*sp.Integer(3)**sp.Rational(1, 4)/3
    deps = -K_*ee**2
    drdu = sp.diff(q*sp.sqrt(ee), ee)*deps                   # = -(q K/2) eps^{3/2}
    kap = ck*sp.sqrt(ee)
    delta_over_r = sp.simplify(-drdu/(kap*rh.subs(L, 1)))     # delta/r_+ with r_+ -> r_h
    kdot = sp.diff(kap, ee)*deps
    ratio = sp.simplify(kdot/kap**2)
    ok2 = sp.simplify(drdu/ee**sp.Rational(3, 2)).free_symbols == {K_}
    ok3 = sp.simplify(delta_over_r/ee).free_symbols == {K_}
    ok4 = sp.simplify(ratio/sp.sqrt(ee)).free_symbols == {K_}
    return ok1 and ok2 and ok3 and ok4, (f"d eps/du = -K eps^2 => eps = {en.subs(C1, C_)} (never 0 at finite u); dr_+/du = {sp.simplify(drdu)} ~ eps^(3/2); "
                                       f"delta/r_+ = {delta_over_r} ~ eps; kappa'/kappa^2 = {ratio} ~ sqrt(eps)")
check("near_Mstar_scaling_chain", c_scaling_chain)


def c_far_field():
    # G=c=hbar=1: dM/du = -P(M); r_+ = 2M, kappa = 1/(4M), delta = -r_+'/kappa; A = P (2M)/M
    Pm = 1/(15360*sp.pi*Mm**2)
    rp_dot = 2*(-Pm)
    kap = 1/(4*Mm)
    delta_over_r = (-rp_dot/kap)/(2*Mm)
    Acal = Pm*(2*Mm)/Mm
    ok = sp.simplify(delta_over_r - 2*Acal) == 0 and sp.simplify(Acal - 1/(7680*sp.pi*Mm**2)) == 0
    return ok, f"delta/r_+ = {sp.simplify(delta_over_r)} = 2 A(M) with A = {sp.simplify(Acal)}"
check("far_field_identity_delta_over_r_equals_2A", c_far_field)


def c_dUdu():
    u_ = sp.Symbol('u'); kap = sp.Function('kappa')(u_); Kc = sp.Function('K')(u_)
    U_ = -sp.exp(-Kc)/kap
    dU = sp.diff(U_, u_).subs(sp.Derivative(Kc, u_), kap)
    ok = sp.simplify(dU - sp.exp(-Kc)*(1 + sp.diff(kap, u_)/kap**2)) == 0
    return ok, "U = -exp(-K)/kappa, K' = kappa: dU/du = exp(-K)[1 + kappa'/kappa^2]"
check("kruskal_like_coordinate_derivative", c_dUdu)

# ---------------------------------------------------------------- Sec 5
def c_rmin():
    t = sp.Symbol('t', positive=True)                       # t = (M/M*)^(1/3)
    Mt = Mstar*t**3
    fp = sp.factor(sp.together(sp.diff(f(r, M), r)))
    num = sp.factor(sp.numer(fp))
    sols = [q for q in sp.solve(num, r) if q.is_positive]
    ok1 = len(sols) == 1 and sp.simplify(sols[0] - (4*M*L**2)**sp.Rational(1, 3)) == 0
    rm = sp.sqrt(3)*L*t
    ok2 = sp.simplify(sp.diff(f(r, Mt), r).subs(r, rm)) == 0
    fm = sp.simplify(f(rm, Mt))
    ok3 = sp.simplify(fm - (1 - t**2)) == 0
    ok4 = sp.simplify(sols[0].subs(M, Mt) - rm) == 0 and sp.simplify((4*M*L**2)**sp.Rational(1, 3) - 2**sp.Rational(2, 3)*M**sp.Rational(1, 3)*L**sp.Rational(2, 3)) == 0
    return ok1 and ok2 and ok3 and ok4, f"f' numerator = {num}: unique positive root r_min = (4 M ell^2)^(1/3) = 2^(2/3) M^(1/3) ell^(2/3); f(r_min) = {fm} = 1-(M/M*)^(2/3) with t=(M/M*)^(1/3)"
check("horizonless_r_min_and_f_min", c_rmin)

check("area_grows_linearly", lambda: (sp.simplify(sp.diff(4*sp.pi*r**2, r) - 8*sp.pi*r) == 0, "d(4 pi r^2)/dr = 8 pi r != 0 for every r>0"))


def c_lightring():
    one = sp.Integer(1)
    fr = f(r, M, one)
    g = fr/r**2
    N = sp.numer(sp.together(sp.diff(g, r)))
    Nf = sp.factor(N)
    # the positive-r polynomial whose roots are the light rings (stationary points of f/r^2)
    core = [fac for fac, mult in sp.factor_list(N)[1] if sp.degree(fac, r) > 1][0]
    disc = sp.factor(sp.discriminant(core, r))
    Mc = [q for q in sp.solve(sp.numer(sp.together(disc)), M) if q.is_positive]      # double root <=> light ring appears/disappears
    ok1 = len(Mc) == 1
    ratio = sp.simplify(Mc[0]/Mstar.subs(L, 1))
    ok1 = ok1 and sp.simplify(ratio - sp.Rational(32, 25)*sp.sqrt(sp.Rational(2, 5))) == 0 and 0.805 < float(ratio) < 0.81
    # number of positive roots just above / below the threshold (Sturm sequence on a rational M)
    def npos(Mv):
        return sp.Poly(core.subs(M, Mv), r).count_roots(sp.Integer(0), sp.oo)
    n_hi = npos(sp.Rational(1, 1)*sp.Rational(6, 5))       # M = 1.2 ell > M_c = 1.0516 ell, M < M* = 1.299 ell
    n_lo = npos(sp.Rational(1, 1))                          # M = 1.0 ell < M_c
    ok2 = (n_hi == 2) and (n_lo == 0)
    # at M*: f/r^2 has its inner minimum (value 0) exactly at r = r_h = sqrt 3
    gs = g.subs(M, Mstar.subs(L, 1))
    ok3 = sp.simplify(gs.subs(r, sp.sqrt(3))) == 0 and sp.simplify(sp.diff(gs, r).subs(r, sp.sqrt(3))) == 0 and sp.simplify(sp.diff(gs, r, 2).subs(r, sp.sqrt(3))) > 0
    Nst = sp.Poly(sp.numer(sp.together(sp.diff(gs, r))), r)
    pr = sorted([sp.re(q) for q in Nst.nroots(n=20, maxsteps=500) if abs(sp.im(q)) < 1e-12 and sp.re(q) > 0])
    return ok1 and ok2 and ok3, (f"light rings = stationary points of f/r^2 (eikonal potential); numerator N(r;M) = {Nf} (ell=1); discriminant_r = {disc}; "
                                 f"double root at M_c = {Mc[0]} ell, i.e. M_c/M* = 32/25 sqrt(2/5) = {sp.N(ratio, 8)}: for M_c < M < M* there are 2 positive roots "
                                 f"(inner minimum + outer maximum of f/r^2), for M < M_c none (checked with Sturm counts {n_hi}, {n_lo}); "
                                 f"at M*: inner minimum f/r^2=0 at r=sqrt3=r_h; stationary points r = {[sp.N(q, 6) for q in pr]}")
check("light_ring_threshold_and_Mstar_rings", c_lightring)

# ---------------------------------------------------------------- Sec 5.4
def c_rw():
    Psi = sp.Function('Psi'); w = sp.Symbol('omega'); ll = sp.Symbol('lam', positive=True)   # lam = l(l+1)
    F = sp.Function('F')(r)
    R = Psi(r)/r
    rad = sp.diff(r**2*F*sp.diff(R, r), r)/r**2 + (w**2/F - ll/r**2)*R      # (1/r^2)(r^2 f R')' + (w^2/f - l(l+1)/r^2) R = 0
    V = F*(ll/r**2 + sp.diff(F, r)/r)
    rw = F*sp.diff(F*sp.diff(Psi(r), r), r) + (w**2 - V)*Psi(r)              # d^2 psi/dx^2 + (w^2 - V) psi,  d/dx = f d/dr
    ok = sp.simplify(sp.expand(rw - F*r*rad)) == 0
    return ok, "for metric -f dt^2 + dr^2/f + r^2 dOmega^2 and Phi = (psi/r) Y e^{-i w t}: d^2psi/dx^2 + (w^2 - V) psi = f r [radial wave operator on Phi]; V = f [l(l+1)/r^2 + f'/r]"
check("regge_wheeler_reduction_of_scalar_wave_equation", c_rw)


def c_vg():
    k, w, Vv = sp.symbols('k omega V', positive=True)
    wk = sp.sqrt(k**2 + Vv)
    vg = sp.diff(wk, k).subs(k, sp.sqrt(w**2 - Vv))
    ok = sp.simplify(vg.subs(Vv, w**2 - sp.Symbol('q', positive=True)**2) - sp.Symbol('q', positive=True)/w) == 0
    return ok, "omega^2 = k^2 + V: group velocity d omega/dk = k/omega = sqrt(omega^2 - V)/omega (round-trip time T = omega_R * 2 int dr/(f sqrt(omega^2-V)))"
check("group_velocity_round_trip_time", c_vg)

n_ok = sum(v["ok"] for v in RES.values())
print(f"\nSYMBOLIC SUMMARY: {n_ok}/{len(RES)} checks passed")
save_json("symbolic_checks.json", RES)

# Verification report - bh-paperA-evaporation-islands

Mode: **quick**.  Date (UTC): 2026-10-21.

## Summary

* checks passed: **95/95**, failed: 0, skipped: 0
* scripts: check_symbolic: OK, check_analytic: OK, wkb: OK, crosscheck_wkb_mpmath: OK, timedomain_nx6000: OK, timedomain_nx3000: OK, timedomain_nx12000: OK, make_fig_sublimation: OK, fit_timedomain: OK, qnm_shooting: OK

## Environment

* python: 3.13.16
* numpy: 2.5.3
* scipy: 1.18.1
* sympy: 1.14.0
* mpmath: 1.3.0
* matplotlib: 3.11.2
* gcc: gcc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 (not used by this repository)
* platform: Linux-6.18.44-fc-v139-x86_64-with-glibc2.39

## Checks (paper value vs value computed by the scripts in `ancillary/`)

| ID | Paper location | Quantity | Paper value | Computed | Tolerance (reason) | Status |
|---|---|---|---|---|---|---|
| S01 | Sec. 2.2 | [sympy, exact] P(r_h)=P'(r_h)=0 at M*=3sqrt3 ell/4, r_h=sqrt3 ell; P'' and dP/dM nonzero ('exact signature of two roots merging') | True | True | exact (exact symbolic identity (check_symbolic.py)) | PASS |
| S02 | Sec. 2.2 | [sympy, exact] discriminant of r^3-2Mr^2+2M ell^2 vanishes only at M*; 3 real roots above, 1 below | True | True | exact (exact symbolic identity (check_symbolic.py)) | PASS |
| S03 | Sec. 2.2 / 5.1 | [sympy, exact] min_(r>0) P = P(4M/3) > 0 iff M<M*: no positive root for M<M* | True | True | exact (exact symbolic identity (check_symbolic.py)) | PASS |
| S04 | Sec. 2.2 | [sympy, exact] A(M*) = 12 pi ell^2 | True | True | exact (exact symbolic identity (check_symbolic.py)) | PASS |
| S05 | Eq. (fold) | [sympy, exact] r_pm - r_h = +-(2 sqrt2 3^(1/4)/3) sqrt(ell eps) to leading order | True | True | exact (exact symbolic identity (check_symbolic.py)) | PASS |
| S06 | Eq. (fold) | [sympy, exact] kappa = c_kappa sqrt(eps)/ell^(3/2), c_kappa = 2 sqrt2 3^(1/4)/9 = 0.4136 | True | True | exact (exact symbolic identity (check_symbolic.py)) | PASS |
| S07 | Sec. 2.2 | [sympy, exact] f ~ (r-r_h)^2/(3 ell^2) at M*; AdS2 radius = r_h = sqrt3 ell | True | True | exact (exact symbolic identity (check_symbolic.py)) | PASS |
| S08 | Sec. 2.2 | [sympy, exact] r_+ = 2M(1-eta/4-eta^2/8-...), kappa/kappa_Schw = 1-eta/2-3eta^2/8-... | True | True | exact (exact symbolic identity (check_symbolic.py)) | PASS |
| S09 | Sec. 2.3, Fig. 1 | [sympy, exact] Gamma = sqrt(f(2M)) = ell/sqrt(4M^2+ell^2) | True | True | exact (exact symbolic identity (check_symbolic.py)) | PASS |
| S10 | Sec. 3.1 | [sympy, exact] P = sigma A_H T_H^4; A(M)=P(r_s/c)/(Mc^2)=hbar c/(7680 pi G M^2); M_break = sqrt(hbar c/(7680 pi G)) | True | True | exact (exact symbolic identity (check_symbolic.py)) | PASS |
| S11 | Sec. 3.1 | [sympy, exact] dM/du=-P gives M ~ (u_f-u)^(1/3), t_evap = 5120 pi G^2 M^3/(hbar c^4) | True | True | exact (exact symbolic identity (check_symbolic.py)) | PASS |
| S12 | Sec. 3.2-3.3 | [sympy, exact] P~eps^2 => eps~1/u (never 0 at finite u); dr_+/du ~ eps^(3/2); delta/r_+ ~ eps; kappa'/kappa^2 ~ sqrt(eps) | True | True | exact (exact symbolic identity (check_symbolic.py)) | PASS |
| S13 | Eq. (far-field-identity) | [sympy, exact] delta/r_+ = 2 A(M) far from M* | True | True | exact (exact symbolic identity (check_symbolic.py)) | PASS |
| S14 | Sec. 4.2 | [sympy, exact] dU/du = exp(-K)[1+kappa'/kappa^2] | True | True | exact (exact symbolic identity (check_symbolic.py)) | PASS |
| S15 | Sec. 5.1 | [sympy, exact] r_min = 2^(2/3) M^(1/3) ell^(2/3), f(r_min) = 1-(M/M*)^(2/3) | True | True | exact (exact symbolic identity (check_symbolic.py)) | PASS |
| S16 | Sec. 5.2 | [sympy, exact] dA/dr = 8 pi r != 0 for r>0 | True | True | exact (exact symbolic identity (check_symbolic.py)) | PASS |
| S17 | Sec. 5.3 | [sympy, exact] light-ring pair (min+max of f/r^2) exists iff M > M_c = (32/25)sqrt(2/5) M* = 0.8095 M*; at M* inner ring at r_h | True | True | exact (exact symbolic identity (check_symbolic.py)) | PASS |
| S18 | Eq. (RWpotential) | [sympy, exact] scalar wave eq. -> psi_xx+(w^2-V)psi=0, V=f[l(l+1)/r^2+f'/r] | True | True | exact (exact symbolic identity (check_symbolic.py)) | PASS |
| S19 | Sec. 5.4 | [sympy, exact] group velocity sqrt(w^2-V)/w, T = w_R * 2 int dr/(f sqrt(w^2-V)) | True | True | exact (exact symbolic identity (check_symbolic.py)) | PASS |
| A01 | Sec. 2.2 | M* for ell=l_P [kg] | 2.83e-08 | 2.82727e-08 | rel 0.003 (paper quotes 3 s.f.) | PASS |
| A02 | Sec. 2.2 | M*/m_P | 1.3 | 1.29904 | abs 0.05 (paper quotes 2 s.f.) | PASS |
| A03 | Sec. 4.4 | m_P [kg] | 2.176e-08 | 2.17643e-08 | rel 0.001 (paper quotes 4 s.f. (CODATA 2018)) | PASS |
| A04 | Eq. (fold) | (r_- - r_h)/[2sqrt2 3^(1/4)/3 sqrt(eps)] at eps=1e-8 | -1 | -0.999916 | abs 0.001 (limit eps->0; O(sqrt(eps)) corrections ~1e-4) | PASS |
| A05 | Eq. (fold) | (r_+ - r_h)/[2sqrt2 3^(1/4)/3 sqrt(eps)] at eps=1e-8 | 1 | 1.00008 | abs 0.001 (limit eps->0; O(sqrt(eps)) corrections ~1e-4) | PASS |
| A06 | Eq. (fold) | kappa/(c_kappa sqrt(eps)) at eps=1e-8 | 1 | 0.999904 | abs 0.001 (limit eps->0; O(sqrt(eps)) corrections ~1e-4) | PASS |
| A07 | Eq. (fold) | c_kappa | 0.41 | 0.413602 | abs 0.005 (paper quotes 2 s.f.) | PASS |
| A08 | Sec. 2.2 | (1/2) f''(r_h) 3 ell^2 at M* (AdS2) | 1 | 1 | abs 1e-06 (numerical second derivative, mpmath 30 digits) | PASS |
| A09 | Sec. 2.2 | r_+/2M-1: exact root vs series (M=1e3 ell) | -2.5e-07 | -2.5e-07 | rel 0.0001 (series truncated at eta^2; numerical derivative) | PASS |
| A10 | Sec. 2.2 | kappa/kappa_S-1: exact vs series (M=1e3 ell) | -5e-07 | -5e-07 | rel 0.0001 (series truncated at eta^2; numerical derivative) | PASS |
| A11 | Sec. 2.2 | r_+/2M-1: exact root vs series (M=1e4 ell) | -2.5e-09 | -2.5e-09 | rel 0.0001 (series truncated at eta^2; numerical derivative) | PASS |
| A12 | Sec. 2.2 | kappa/kappa_S-1: exact vs series (M=1e4 ell) | -5e-09 | -5e-09 | rel 0.001 (series truncated at eta^2; numerical derivative (value 5e-9: finite-difference noise ~1e-12)) | PASS |
| A13 | Sec. 2.2 | eta = (l_P/(G 10Msun/c^2))^2 | 1.2e-78 | 1.19799e-78 | rel 0.02 (paper quotes 2 s.f.) | PASS |
| A14 | Sec. 3.1 | Hawking power P(10 Msun) [W] | 9e-31 | 9.00761e-31 | rel 0.01 (paper quotes 2 s.f.) | PASS |
| A15 | Sec. 3.1 | T_H(10 Msun) [K] | 6e-09 | 6.17007e-09 | rel 0.05 (paper quotes 1 s.f. (6.17e-9)) | PASS |
| A16 | Sec. 3.1 | adiabaticity A(10 Msun) | 5e-83 | 4.96526e-83 | rel 0.01 (paper quotes 1 s.f. (4.97e-83)) | PASS |
| A17 | Eq. (adiabaticity-param) | P (r_s/c)/(M c^2) equals hbar c/(7680 pi G M^2) | 4.96526e-83 | 4.96526e-83 | rel 1e-06 (same number two ways) | PASS |
| A18 | Sec. 3.1 | M_break [kg] | 1.4e-10 | 1.40117e-10 | rel 0.01 (paper quotes 2 s.f.) | PASS |
| A19 | Sec. 3.1 | M_break/m_P | 0.0064 | 0.0064379 | rel 0.01 (paper quotes 2 s.f.) | PASS |
| A20 | Sec. 3.1 | M*/M_break ('more than two orders of magnitude') |  | 201.78 | in [100.0, None] (qualitative statement; value 202) | PASS |
| A21 | Sec. 3.1 | A(M*) (Schwarzschild expression at threshold) | 2.5e-05 | 2.45609e-05 | rel 0.02 (paper quotes 2 s.f.) | PASS |
| A22 | Eq. (far-field-identity) | delta/r_+ = |r_+'|/(kappa r_+) vs 2A (SI, 10 Msun) | 9.93052e-83 | 9.93052e-83 | rel 1e-06 (same number two ways) | PASS |
| A23 | Sec. 3.2 | eps_c = (c_kappa/2pi)^2 [Planck masses] | 0.004 | 0.00433317 | rel 0.1 (paper: 'a crude estimate', 1 s.f. (4.33e-3)) | PASS |
| A24 | Sec. 4.4 | S_BH(10 Msun) [bits] ('~1e79') | 1e+79 | 1.51332e+79 | |dlog10| 0.5 (order-of-magnitude statement) | PASS |
| A25 | Sec. 4.4 | E = M* c^2 [J] | 2.5e+09 | 2.54102e+09 | rel 0.02 (paper quotes 2 s.f.) | PASS |
| A26 | Sec. 2.3 | r_s/c for 10 Msun [s] ('order 1e-4 s') | 0.0001 | 9.85128e-05 | |dlog10| 0.5 (order-of-magnitude statement) | PASS |
| A27 | Sec. 2.2 | Hawking evaporation time of 10 Msun [yr] ('~1e70 years') | 1e+70 | 2.09568e+70 | |dlog10| 0.5 (order-of-magnitude statement) | PASS |
| A28 | Sec. 4.4 | Bekenstein bound for M* c^2 within R=l_P [bits] ('a handful') |  | 11.7754 | in [1.0, 20.0] (qualitative ('a handful'); computed value is ~12 bits) | PASS |
| A29 | Sec. 5.1 | r_min numeric vs 2^(2/3)M^(1/3) (M=0.5 M*, ell=1) | 1.37473 | 1.37473 | rel 0.0001 (grid search, spacing 1e-5) | PASS |
| A30 | Sec. 5.1 | f(r_min) numeric vs 1-(M/M*)^(2/3) (M=0.5 M*) | 0.370039 | 0.370039 | rel 1e-05 (grid search) | PASS |
| A31 | Sec. 5.1 | r_min numeric vs 2^(2/3)M^(1/3) (M=0.9 M*, ell=1) | 1.67228 | 1.67228 | rel 0.0001 (grid search, spacing 1e-5) | PASS |
| A32 | Sec. 5.1 | f(r_min) numeric vs 1-(M/M*)^(2/3) (M=0.9 M*) | 0.0678302 | 0.0678302 | rel 1e-05 (grid search) | PASS |
| A33 | Sec. 5.1 | r_min numeric vs 2^(2/3)M^(1/3) (M=0.99 M*, ell=1) | 1.72626 | 1.72625 | rel 0.0001 (grid search, spacing 1e-5) | PASS |
| A34 | Sec. 5.1 | f(r_min) numeric vs 1-(M/M*)^(2/3) (M=0.99 M*) | 0.00667783 | 0.00667783 | rel 1e-05 (grid search) | PASS |
| A35 | Sec. 5.3 | smallest M/M* with inner light ring (bisection; paper: 0.81 M* <~ M) | 0.81 | 0.809543 | in [0.805, 0.81] (paper quotes 0.81 M* (2 s.f., '<~'); exact value 0.80954 (S17)) | PASS |
| A36 | Sec. 5.3 | inner light-ring radius at 0.95 M* [ell] (paper 1.73-1.77) |  | 1.76664 | in [1.73, 1.77] (range quoted in the paper) | PASS |
| A37 | Sec. 5.3 | inner light-ring radius at M* [ell] (paper 1.73-1.77) |  | 1.73204 | in [1.73, 1.77] (range quoted in the paper) | PASS |
| A38 | Sec. 5.3 | outer maximum of f/r^2 near M* [ell^-2] ('about 0.02') |  | 0.0250836 | in [0.02, 0.03] ('about 0.02' (value 0.025)) | PASS |
| A39 | Sec. 5.3 | inner minimum of f/r^2 near M* ('nearly zero') |  | 0.00222581 | in [0.0, 0.005] ('nearly zero' (value 0.0022)) | PASS |
| F01 | Fig. 1 | plotted Gamma=sqrt(f(2M)) equals closed form ell/sqrt(4M^2+ell^2) (max abs deviation) | 0 | 2.49093e-13 | abs 1e-09 (floating point) | PASS |
| F02 | Fig. 1 | Gamma increases monotonically as M decreases (1 -> 1e-4) | True | True | exact (shape of curve) | PASS |
| F03 | Fig. 1 | regenerated PNG has the pixel size of the shipped sublimation_curve.png | True | True | exact (1000x600) | PASS |
| F04 | Fig. 1 | mean abs RGB difference to the shipped PNG | 0 | 0.0331088 | abs 0.06 (same curve; text antialiasing/font hinting differ between matplotlib versions) | PASS |
| W01 | Table 1, row 1 | LaTeX row regenerated by wkb.py equals the row in main.tex | row present in main.tex | $0.90$ & $1$ & $0.277$ & $2.3\times10^{1}$ & $1.7\times10^{2}$ \\ | exact string (exact string match of Re(w) (3 d.p.), Q and tau (2 s.f.)) | PASS |
| W02 | Table 1, row 2 | LaTeX row regenerated by wkb.py equals the row in main.tex | row present in main.tex | $0.90$ & $2$ & $0.430$ & $7.9\times10^{1}$ & $3.7\times10^{2}$ \\ | exact string (exact string match of Re(w) (3 d.p.), Q and tau (2 s.f.)) | PASS |
| W03 | Table 1, row 3 | LaTeX row regenerated by wkb.py equals the row in main.tex | row present in main.tex | $0.90$ & $4$ & $0.732$ & $1.7\times10^{3}$ & $4.6\times10^{3}$ \\ | exact string (exact string match of Re(w) (3 d.p.), Q and tau (2 s.f.)) | PASS |
| W04 | Table 1, row 4 | LaTeX row regenerated by wkb.py equals the row in main.tex | row present in main.tex | $0.95$ & $1$ & $0.200$ & $3.6\times10^{2}$ & $3.6\times10^{3}$ \\ | exact string (exact string match of Re(w) (3 d.p.), Q and tau (2 s.f.)) | PASS |
| W05 | Table 1, row 5 | LaTeX row regenerated by wkb.py equals the row in main.tex | row present in main.tex | $0.95$ & $2$ & $0.309$ & $1.0\times10^{4}$ & $6.7\times10^{4}$ \\ | exact string (exact string match of Re(w) (3 d.p.), Q and tau (2 s.f.)) | PASS |
| W06 | Table 1, row 6 | LaTeX row regenerated by wkb.py equals the row in main.tex | row present in main.tex | $0.99$ & $1$ & $0.090$ & $2.4\times10^{5}$ & $5.4\times10^{6}$ \\ | exact string (exact string match of Re(w) (3 d.p.), Q and tau (2 s.f.)) | PASS |
| W07 | Table 1 caption / Sec. 5.4 | well-to-barrier gap at l=5, M=0.90 M* [ell^-2] | 0.28 | 0.283984 | abs 0.005 (paper quotes 2 d.p.) | PASS |
| W08 | Table 1 caption | well-to-barrier gap at l=5, M=0.99 M* [ell^-2] | 0.7 | 0.700237 | abs 0.005 (paper quotes 2 d.p.) | PASS |
| W09 | Sec. 5.4 | Q (WKB) at l=1, M=0.90 M* ('Q ~ 20') | 20 | 23.0382 | rel 0.2 ('~20') | PASS |
| W10 | Sec. 5.4 | tau (WKB) at l=1, M=0.90 M* [ell/c] ('~1e2') | 100 | 166.164 | |dlog10| 0.3 ('~1e2') | PASS |
| W11 | Table 1 | scipy vs mpmath(30 digits, tanh-sinh): max rel. difference of Q,tau,Im(w),S over the six rows | 0 | 3.91476e-10 | abs 1e-06 (independent quadrature of the singular endpoint integrals) | PASS |
| W12 | Table 1 | Bohr-Sommerfeld residual of the six modes (mpmath) | 0 | 2.0192e-09 | abs 1e-07 (E solved by brentq to xtol 1e-10) | PASS |
| T01 | Sec. 5.4 | time-domain omega_R (l=2, 0.90 M*, Nx=6000) | 0.4296 | 0.429613 | abs 6e-05 (paper quotes 4 d.p.) | PASS |
| T02 | Sec. 5.4 | time-domain decay rate [1/ell] | 0.00308 | 0.00308444 | rel 0.002 (paper quotes 3 s.f.) | PASS |
| T03 | Sec. 5.4 | WKB omega_R | 0.4303 | 0.43032 | abs 6e-05 (paper quotes 4 d.p.) | PASS |
| T04 | Sec. 5.4 | WKB |Im omega| | 0.00271 | 0.00271173 | rel 0.002 (paper quotes 3 s.f.) | PASS |
| T05 | Sec. 5.4 | |omega_R(TD)-omega_R(WKB)|/omega_R(WKB) ('0.2%') | 0.002 | 0.00164458 | abs 0.0005 (paper quotes 1 s.f.) | PASS |
| T06 | Sec. 5.4 | decay(TD)/|Im omega|(WKB) ('a factor of 1.14') | 1.14 | 1.13744 | abs 0.005 (paper quotes 3 s.f.) | PASS |
| T07 | Sec. 5.4 | max change of TD decay rate between Nx=6000 and Nx=3000 / 12000 ('<0.5%') |  | 0.0019222 | in [0.0, 0.005] ('less than 0.5%') | PASS |
| T08 | Sec. 5.4 | max change of TD omega_R between Nx=6000 and Nx=3000 / 12000 |  | 0.000497759 | in [0.0, 0.005] (grid convergence (not quoted numerically)) | PASS |
| Q01 | Sec. 5.4 | direct Re(omega), l=2, 0.90 M* | 0.4296 | 0.429629 | abs 6e-05 (paper quotes 4 d.p.) | PASS |
| Q02 | Sec. 5.4 | direct -Im(omega), l=2, 0.90 M* | 0.00307 | 0.00307192 | abs 6e-06 (paper quotes 3 s.f. (sign: Im<0, compared as |Im|)) | PASS |
| Q03 | Sec. 5.4 | direct Re(omega) agrees with time domain (abs. diff) | 0 | 1.62837e-05 | abs 0.001 ('in agreement') | PASS |
| Q04 | Sec. 5.4 | direct |Im omega| agrees with time-domain decay rate (rel. diff) | 0 | 0.00405644 | abs 0.01 ('in agreement') | PASS |
| Q05 | Sec. 5.4 | direct Q (M/M*, l) = 0.90, 1 | 17 | 17.3769 | rel 0.03 (paper quotes 2 s.f.) | PASS |
| Q06 | Sec. 5.4 | direct Q (M/M*, l) = 0.90, 2 | 70 | 69.9283 | rel 0.01 (paper quotes 2 s.f.) | PASS |
| Q07 | Sec. 5.4 | direct Q (M/M*, l) = 0.90, 4 | 1400 | 1393.12 | rel 0.04 (paper quotes 2 s.f.) | PASS |
| Q08 | Sec. 5.4 | direct Q (M/M*, l) = 0.95, 1 | 200 | 200.464 | rel 0.03 (paper quotes 2 s.f.) | PASS |
| Q09 | Sec. 5.4 | direct Q (M/M*, l) = 0.95, 2 | 6800 | 6775.93 | rel 0.01 (paper quotes 2 s.f.) | PASS |
| Q10 | Sec. 5.4 | direct Q (M/M*, l) = 0.99, 1 | 160000 | 160967 | rel 0.04 (paper quotes 2 s.f.) | PASS |
| Q11 | Sec. 5.4 | min over the six Table-1 entries of |Im|_direct/|Im|_WKB ('between 1.1 and 1.8') | 1.1 | 1.13283 | abs 0.05 (paper quotes 2 s.f.) | PASS |
| Q12 | Sec. 5.4 | max over the six Table-1 entries of |Im|_direct/|Im|_WKB ('between 1.1 and 1.8') | 1.8 | 1.81574 | abs 0.05 (paper quotes 2 s.f. (value 1.82 rounds to 1.8)) | PASS |
| Q13 | Sec. 5.4 | robustness: Re/Im(omega) with matching point r_m=1.8 vs 2.0 (max abs diff) | 0 | 1.16787e-12 | abs 1e-06 (result must not depend on r_m) | PASS |

Symbolic (S01-S19) rows are exact sympy derivations, see `ancillary/output/symbolic_checks.txt` for the derived expressions.

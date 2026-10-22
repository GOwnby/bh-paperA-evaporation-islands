Code for the numerical and symbolic results of the paper.  Run everything with `bash ../verify.sh` (see ../README.md for the
result -> script -> output table); outputs go to output/ (relative to the script), checks to ../verification/.

_util.py                  shared output helper (stdout tee into output/, JSON writer); no numerical code
wkb.py                    three-turning-point WKB (Bohr-Sommerfeld + Gamow) trapped-mode spectrum; Table 1            [scipy]
                          (Gamow prefactor corrected relative to the original, see PROVENANCE.md; original: wkb_original.py)
wkb_original.py           the author's original script (v1.0.1), verbatim; prints the v1.0.1 spectrum, writes rows.pkl
timedomain_original.py    the author's original script (v1.0.1), verbatim; same signal as timedomain.py (bit-identical)
compare_wkb_original.py   original vs corrected WKB on Table 1, and both against the direct integration (added after v1.0.1)  [scipy]
crosscheck_wkb_mpmath.py  independent mpmath (30 digit) recomputation of the WKB integrals of Table 1                  [scipy, mpmath]
timedomain.py             finite-difference time-domain evolution of the same wave equation (l=2, M=0.90 M*)          [scipy]
fit_timedomain.py         omega_R and decay rate from the time-domain signals                                          [scipy]
qnm_shooting.py           direct complex-frequency integration (non-WKB cross-check)                                    [scipy]
make_fig_sublimation.py   Fig. 1                                                                                        [matplotlib]
check_analytic.py         numerical (mpmath) checks of the analytic scalings and quoted numbers                         [mpmath]
check_symbolic.py         exact (sympy) derivations of the analytic statements                                          [sympy]

Individual scripts can be run from any directory, e.g.  python3 wkb.py

Provenance (details and evidence: ../PROVENANCE.md)
  original scripts of the author (v1.0.1)       : wkb_original.py, timedomain_original.py (and Fig. 1, sublimation_curve.png)
  changed relative to the original              : wkb.py (Gamow prefactor corrected: Im(w) = -exp(-2S)/(2 w_R Ts); Table 1 and Sec. 5.4 changed),
                                                  timedomain.py (arguments, output format; bit-identical signal)
  written after v1.0.1, no earlier version      : fit_timedomain.py (the v1.0.1 fit was done by hand, code lost), qnm_shooting.py (new, non-WKB check),
                                                  compare_wkb_original.py, crosscheck_wkb_mpmath.py, check_analytic.py, check_symbolic.py, _util.py
  re-creation of an unpreserved generator       : make_fig_sublimation.py (the shipped sublimation_curve.png is the original; the regenerated one has the same
                                                  curve, axes shifted by ~2 px)

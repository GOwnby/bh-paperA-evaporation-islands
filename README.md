# Evaporation Toward a Degenerate Horizon: Islands, Asymptotics and the Horizonless Branch of the Hayward Black Hole

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23455575.svg)](https://doi.org/10.5281/zenodo.23455575)

**Author:** Gracen Ownby  
**Permanent Zenodo archive:** <https://doi.org/10.5281/zenodo.23455575>  
**License:** code MIT, manuscript/figures/data CC-BY-4.0 (see `LICENSE.md`); how to cite: `CITATION.cff`.

## Abstract

We study the late-stage evaporation of a Hayward-type regular black hole whose regularization scale is tied to the mass so that the core density is mass-independent, eps^3 = 2GM l^2/c^2. The metric has a critical mass M* = (3 sqrt(3)/4) l c^2/G at which its inner and outer horizons merge. Far from M* it reproduces Schwarzschild to fractional order (l/M)^2, so to a distant observer the object is an extremely time-dilated collapsed star. We make four points. First, the quasi-static (eternal black hole) approximation used in the island literature is justified along the evaporation: an adiabaticity parameter is of order 10^-83 for a 10 solar-mass hole, is small at M*, and, with the near-extremal mass-loss law, vanishes as the endpoint is approached. Second, within the continuous semiclassical mass-loss law the transition at M* is never completed in finite time. The horizon area at M* is a fixed nonzero constant while the surface gravity vanishes, so the emission power scales as eps^2 in the mass excess eps = M - M* and eps ~ 1/u. This is the analogue for this metric of the unattainability of extremality known for charged black holes, and infinite evaporation times are known for regular black holes in general; it is not a new structural statement. The continuous law fails in a small window above M*, where the emission of a single quantum or near-extremal quantum- gravity effects become important. Third, the island mechanism is argued, at leading adiabatic order and through a heuristic extension of the static prescription, to remain intact along the approach. The entropy plateau approaches, at leading order, the area value 3 pi l^2/G, the leading term of the extremal Hayward island entropy of Kim and Nam (2021); the extremal result is not the limit of the non-extremal one, and we do not resolve the subleading terms. Fourth, for M < M* the geometry has no horizon and supports no island, by one argument in two steps: a purity argument reduces the extremization to the area term, and the area then grows linearly near the point of maximal redshift, leaving no critical point. For M >~ 0.81 M* it also has a stable inner light ring, as a general theorem requires for smooth ultracompact objects satisfying the null energy condition; a test-scalar WKB calculation, cross-validated against a time-domain evolution, finds genuinely long-lived trapped modes there, consistent with the expectation that the horizonless configuration is not a stable endpoint. Whether general relativity coupled to standard matter can supply the core is treated in a companion paper.

## Directory map

| Path | Content |
|---|---|
| `main.tex`, `references.bib`, `main.bbl` | LaTeX source, bibliography, and the bibliography as built |
| `paperA_evaporation_islands.pdf` | compiled manuscript (17 pages) |
| `sublimation_curve.png` | Fig. 1 of the paper |
| `ancillary/` | all code (Python 3) used for the numerical and symbolic results; see below |
| `ancillary/output/` | outputs of the scripts (text logs, JSON results, regenerated Fig. 1) |
| `verify.sh`, `verify.py`, `verification/` | one-command reproduction; `verification/` holds `expected.json` (paper values), `REPORT.md`, logs |
| `requirements.txt` | Python dependencies (minimum versions) |
| `PROVENANCE.md` | which files are the author's original (Zenodo v1.0.1), which differ, which were written later, and which quoted numbers changed and why |
| `LICENSE.md`, `CITATION.cff` | licence and citation metadata |

## Build the PDF

```
pdflatex main
bibtex main
pdflatex main
pdflatex main
pdflatex main     # third pass removes the "Label(s) may have changed" notice
```

(TeX Live 2022 or newer; packages used: amsmath, amssymb, graphicx, hyperref, natbib, geometry.)

## Verify (reproduce every number, table and figure)

```
pip install -r requirements.txt      # numpy>=2, scipy>=1.11, sympy, mpmath, matplotlib (no compiler, no network)
bash verify.sh            # --quick (default): all generators, about 3 min of CPU time (measured: 2.5 min CPU, 14 min wall-clock on a heavily shared 2-core machine)
bash verify.sh --full     # additionally a finer time-domain grid and halved RK4 steps (about 25 min of CPU time; the shipped `verification/REPORT_full.md` is from such a run: 97/97 checks)
```

`verify.sh` runs every script in `ancillary/` (it uses `${PYTHON:-python3}`, Python >= 3.10; Linux, macOS, or Windows via WSL),
compares the computed values with the values quoted in `main.tex` (`verification/expected.json`: paper value, tolerance and reason),
rewrites `ancillary/output/` and writes `verification/REPORT.md` (PASS/FAIL table and environment), `verification/logs/<script>.log`,
`verification/checks.json` and `verification/timings.txt` (wall-clock per script). The last line is `RESULT bh-paperA-evaporation-islands: PASSED n/n checks (quick)`;
the exit code is non-zero if a check fails. If a Python module is missing the affected scripts are reported as SKIPPED and the run
is flagged INCOMPLETE (exit code 3). The scripts can also be run one by one from `ancillary/` (e.g. `python3 wkb.py`); they find their
output directory relative to the script. The time-domain signals `ancillary/output/td_*.npz` (3-12 MB) are regenerated, not tracked.

### Result in paper -> script -> output file

| Result in the paper | Script | Output file (`ancillary/output/`) |
|---|---|---|
| Table 1 (WKB: Re(omega), Q, tau; rows printed ready for LaTeX under "Table 1 rows"; `verify.py` compares them string-for-string with `main.tex`) | `wkb.py` | `wkb_output.txt`, `wkb.json` |
| Table 1: independent mpmath recomputation of the WKB integrals | `crosscheck_wkb_mpmath.py` | `crosscheck_wkb_mpmath.txt` |
| Sec. 5.4 / Table 1 caption: well-to-barrier gap 0.28 -> 0.70 at l=5 | `wkb.py` | `wkb_output.txt` |
| Sec. 5.4: WKB value 0.4303 / 2.71e-3 (l=2, M=0.90 M*); time-domain 0.4296 +- 0.0003 / 3.08e-3 (envelope fit and damped-sinusoid fit); grid change <0.5% | `timedomain.py` (Nx = 3000, 6000, 12000), `fit_timedomain.py` (written after v1.0.1: the original fit was done by hand) | `timedomain_fit.txt` |
| Sec. 5.4: correction of the Gamow prefactor (v1.0.1 used 1/Ts instead of 1/(2 w_R Ts)): the original script reproduces the v1.0.1 Table 1, Im(w) differs by exactly 2 w_R; direct integration supports the corrected version | `compare_wkb_original.py` (uses `wkb_original.py`, `qnm_shooting.py`) | `compare_wkb_original.txt`, `.json` |
| Sec. 5.4: direct complex-frequency values (0.4296 - 0.00307 i; Q values; factor 1.1-1.8 vs WKB); added after v1.0.1, no earlier script | `qnm_shooting.py` | `qnm_shooting.txt` |
| Fig. 1 (time-dilation factor vs remaining mass) | `make_fig_sublimation.py` | `sublimation_curve.png` (same curve as the root `sublimation_curve.png`), `sublimation_curve.json` |
| Secs. 2-3, 4.4: M*, c_kappa, fold scalings, eta expansions, eta for 10 solar masses, P, T_H, A(10 Msun), A(M*), M_break, far-field identity, epsilon_c, E = M* c^2, bits of S_BH | `check_analytic.py` (numerical, mpmath) | `analytic_checks.txt`, `analytic_checks.json` |
| Same statements, exact derivations (double root and fold, c_kappa, eta series, Hawking-law identities, scaling chain eps ~ 1/u, r_min, f(r_min), no horizon for M<M*, light-ring threshold, Regge-Wheeler reduction) | `check_symbolic.py` (sympy) | `symbolic_checks.txt`, `symbolic_checks.json` |
| Sec. 5: r_min, f(r_min); inner light ring for M > 0.81 M* (exact threshold 0.8095 M*), its radius, outer maximum of f/r^2 | `check_analytic.py`, `check_symbolic.py` | `analytic_checks.txt`, `symbolic_checks.txt` |

All other statements of the paper (the island and Page-curve discussion, the scaling arguments of Sec. 3.2-3.3, the
purity/area-growth argument of Sec. 5.2) are analytic or taken from the cited literature; the algebra behind them is in `check_symbolic.py`.

### Notes

* All numerical code requires SciPy (`quad`, `brentq`, `cumulative_trapezoid`, `interp1d` with cubic interpolation); earlier
  versions of this archive contained a NumPy/mpmath fallback (`_compat.py`) which has been removed. The reference outputs and
  `verification/` were produced with the library versions listed in `verification/REPORT.md`.
* The WKB integrals have square-root endpoint behaviour at the turning points; they are evaluated after the exact substitution
  r = (a+b)/2 - (b-a)/2 cos(theta) and cross-checked against mpmath tanh-sinh quadrature (agreement better than 1e-9).
* The Gamow decay rate in `wkb.py` uses the true round-trip time T = omega_R * 2 int dr/(f sqrt(omega^2 - V)) and
  Im(omega) = -exp(-2 S)/(2T), as in Sec. 5.4 of the paper.
* `ancillary/wkb_original.py` and `ancillary/timedomain_original.py` are the author's original scripts (v1.0.1), verbatim, for reference; they are not run by `verify.sh`
  except that `compare_wkb_original.py` loads the functions of `wkb_original.py`. `timedomain.py` gives a bit-identical signal to `timedomain_original.py`.
  `wkb.py` differs from the original in the Gamow prefactor (Table 1 and Sec. 5.4 changed accordingly); see `PROVENANCE.md`.
* No external data are used.

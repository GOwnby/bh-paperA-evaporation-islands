# PROVENANCE of this repository (Paper A: evaporation toward a degenerate horizon, islands, horizonless branch)

## 0 Baselines used in this file

Three states are compared: (1) the author's tree `old_archive/1_sources/sources/<paper>`; (2) the first Zenodo release (the initial git commit of this repository); (3) the current repository. On 2026-10-21 (1) and (2) were compared byte for byte: they are identical except for the repository metadata (`.gitignore`, `.zenodo.json`, `CITATION.cff`), the compiled PDF. The rest of the old archive (`old_archive/2_docs_and_archive_1/archive/`: earlier manuscript versions, session logs, audit directories) was searched for code and data behind this paper's results; files found there that are not in the Zenodo release are listed below where they matter. Wherever this file says "original" or "v1.0.1" it means (1) = (2).

## 1. Lineage

* This repository is a reorganization of the author's original archive of Paper A, which is the content of Zenodo v1.0.1. That archive
  contained `main.tex`, `references.bib`, `sublimation_curve.png` and two scripts, `ancillary/wkb.py` and `ancillary/timedomain.py`.
* The reorganization added `verify.sh`, `verify.py`, `verification/`, `README.md`, `ancillary/output/`, further scripts, and edited the two original scripts and the manuscript.
* On 2026-10-21 the repository was compared file by file with the original archive (and with the older versions of the papers and scripts in the author's
  archive of earlier sessions) in an audit by Claude (Anthropic). Everything below that says "run" was executed in that audit; what was only read is marked "read".
* Standing rule applied: the original scripts and results are not silently replaced. Where this repository differs from the original, it says so here.

## 2. Files that differ from the original

Classes: a = cosmetic or additive (paths, arguments, comments, output format); b = bug fix without a visible effect on results; c = changes a quoted result.

| path | class | what changed | evidence |
|---|---|---|---|
| `ancillary/wkb.py` | a | docstring, import of `_util`, reformatting; `E` also returned; more (l, M) printed; Table-1 LaTeX rows, gap block and JSON output added; `rows.pkl` no longer written | the original printed no rows for (0.95, l=1,2) and (0.99, l=1) and no l=5 gap; both obtained from its own `solve_mode`/`well_barrier` give the numbers of v1.0.1 |
| `ancillary/wkb.py` | b | f' analytic instead of central difference (h=1e-6); the three integrals evaluated after r = (a+b)/2 - (b-a)/2 cos(theta); floor 1e-14 -> 1e-300 | V_l differs by <= 4.5e-10; round-trip integral by 1.2e-11; the original emits an `IntegrationWarning` (quad roundoff) |
| `ancillary/wkb.py` | **c** | **Gamow prefactor.** Original: Im(w) = -exp(-2S)/Ts with Ts = 2 int dr/(f sqrt(E-V)). Here: T = w_R Ts, Im(w) = -exp(-2S)/(2T) | Im_new = Im_old/(2 w_R) to 1.1e-7 on all six Table-1 rows, Re(w) and S unchanged (5e-12, 1e-9) (`compare_wkb_original.py`, checks P01-P05). Ts has dimension length^2, so 1/Ts is not a rate; the group velocity is sqrt(w^2-V)/w, the round-trip time is w_R Ts. Direct complex-frequency integration (`qnm_shooting.py`), M=0.90 M*, l=1..7: |Im|_direct/|Im|_corrected = 1.33, 1.13, 1.17, 1.21, 1.22, 1.21, 1.20; /|Im|_original = 2.40, 1.32, 1.01, 0.83, 0.69, 0.59, 0.51. The corrected formula is that of Cardoso et al. (2014) with the oscillating factor replaced by its average (as stated in Sec. 5.4; the paper was read from the manuscript, not re-derived here) |
| `ancillary/timedomain.py` | a | command-line arguments (l, M/M*, Nx, T), `.npz` output instead of `signal.npy`/`tarr.npy`, docstring | default run (l=2, M=0.90 M*, Nx=6000, T=2200) of the original and of this version give a bit-identical signal and time array (max abs difference 0.0) |
| `ancillary/README.txt` | a | rewritten file list and provenance section | text only |
| `sublimation_curve.png` | none | identical to the original (md5 29fcf882...) | the curve read from its pixels agrees with Gamma = sqrt(f(2M)) for ell = 1e-3 to 0.005 (2 px) |
| `main.tex` | c | numbers and wording changed relative to v1.0.1, see section 4 | |

The originals are shipped unchanged as `ancillary/wkb_original.py` and `ancillary/timedomain_original.py`.

## 3. Files without an original

None of these exists in the original archive or in any earlier version found in the author's archives (searched by name and by content).

| path | what it generates | basis | quoted numbers that depend on it (v1.0.1 value where there was one) |
|---|---|---|---|
| `ancillary/fit_timedomain.py` | omega_R and decay rate of the time-domain signal (FFT peak; envelope fit with 10% trim; damped-sinusoid least squares) | re-creation: the v1.0.1 fit was "done by hand in the paper" (old README), code not preserved | time-domain omega_R (v1.0.1: 0.4298; now 0.4296 +- 0.0003) and decay rate 3.08e-3 (v1.0.1: 3.08e-3, reproduced), ratio to WKB (v1.0.1: 1.3; now 1.14 because the WKB value changed), "<0.5% under grid change" (new) |
| `ancillary/qnm_shooting.py` | direct complex-frequency integration (regular solution outward, outgoing solution along a rotated contour, secant iteration) | new method, no earlier version; standard complex-contour/shooting technique; written after v1.0.1 | the paragraph of Sec. 5.4 on the direct integration: 0.4296 - 0.00307 i; Q = 17, 70, 1.4e3, 2.0e2, 6.8e3, 1.6e5; "factors between 1.1 and 1.8"; none of this exists in v1.0.1 |
| `ancillary/compare_wkb_original.py` | original vs corrected WKB on the Table-1 rows; both against the direct integration | written for this audit; uses `wkb_original.py` verbatim | justifies the change of Table 1 (section 4); no number of the paper |
| `ancillary/make_fig_sublimation.py` | Fig. 1 | re-creation; the generator of the shipped (original) PNG is not preserved | Fig. 1 only (the figure file is the original) |
| `ancillary/crosscheck_wkb_mpmath.py` | mpmath (30 digit) recomputation of the WKB integrals | verification script, no original | Table 1 integrals (checks W11, W12); read, not run in the audit (mpmath not available) |
| `ancillary/check_analytic.py` | numerical check of the SI/analytic numbers (M*, c_kappa, eta, P, T_H, A, M_break, eps_c, E, bits, r_min, light rings) | verification script; v1.0.1 states these numbers without a script | A01-A39; all reproduced independently with a numpy script in the audit (run), the script itself only read |
| `ancillary/check_symbolic.py` | sympy derivations of the analytic statements | verification script, no original | S01-S19 and the sentence "exactly (32/25) sqrt(2/5) M*"; read only (sympy not available in the audit); the threshold 0.809543 was confirmed numerically (run) |
| `ancillary/_util.py` | output helper | no numerics | none |
| `ancillary/output/*`, `verification/*` | regenerated outputs and check reports | produced by `verify.sh` | |

## 4. Results that changed relative to the original (v1.0.1) and why

| quantity | v1.0.1 | this version | reason / backing |
|---|---|---|---|
| Table 1, Q and tau, all six rows | 4.2e1/3.0e2, 9.2e1/4.3e2, 1.2e3/3.1e3, 8.9e2/8.9e3, 1.7e4/1.1e5, 1.3e6/3.0e7 | 2.3e1/1.7e2, 7.9e1/3.7e2, 1.7e3/4.6e3, 3.6e2/3.6e3, 1.0e4/6.7e4, 2.4e5/5.4e6 | Gamow prefactor (section 2). Re(w) column unchanged. The v1.0.1 rows are reproduced by `wkb_original.py` (P01) |
| WKB decay rate, l=2, 0.90 M* | 2.33e-3 | 2.71e-3 | same; "factor 1.3" -> "1.14" |
| time-domain omega_R | 0.4298 | 0.4296 +- 0.0003 | not reproducible from the original signal: FFT peak 0.42961, least squares 0.42965; FFT bin width 2.5e-4, peak 0.4294-0.4297 over grids |
| omega_R mismatch to WKB | 0.1% | 0.2% | 0.4296 vs 0.4303 (0.16%) |
| Sec. 5.4 direct integration, "factor 1.1-1.8", Q after correction | absent | added | `qnm_shooting.py`, checks Q01-Q14; labelled as added after v1.0.1 in the text |
| eps_c (Sec. 3.2) | "order 1e-2" | "about 4e-3" | (c_kappa/2pi)^2 = 4.33e-3 (A23); the old number was an order of magnitude off |
| remnant energy (Sec. 4.4) | m_P c^2 = 1.95e9 J | M* c^2 = 2.5e9 J, M* = 1.3 m_P | A01, A02, A25 |
| Bekenstein bound for M* c^2 in one Planck length | "a handful of bits" | about 12 (11.8) | A28 |
| outer maximum of f/r^2 near M* (Sec. 5.3) | about 0.02 | about 0.025 | computed 0.0251 at M* (A38) |
| light-ring threshold | 0.81 M* | 0.81 M*, "exactly 0.8095 M*" | S17, A35 |
| Fig. 1 caption | no ell | ell = 1e-3 in units of the largest mass | consistent with the PNG (checked from pixels) |

## 5. What reproduces from the ORIGINAL code (runs of the audit)

* `wkb_original.py` (3.3 s): all six rows of Table 1 of v1.0.1 (Q, tau to 2 significant figures, Re(w) to 3 decimals), l=1..4 at 0.90 M* (Q = 41.5, 92.2, 302, 1152), gaps at l=5: 0.2840 and 0.7002 (paper 0.28, 0.70). Check P01.
* `timedomain_original.py` (25 s, Nx=6000): signal identical to that of `timedomain.py`. Fits of this signal: decay 3.0844e-3 (envelope), 3.0814e-3 (damped-sinusoid least squares), v1.0.1 value 3.08e-3; omega_R 0.42961 (FFT) and 0.42965 (least squares), v1.0.1 value 0.4298 (not reproduced).
* `sublimation_curve.png`: unchanged; matches the formula (see section 2). `make_fig_sublimation.py` regenerates the same curve (3096 vs 3092 curve pixels) with axes shifted by ~2 px (different matplotlib).
* Analytic and SI numbers of the paper (no script in v1.0.1): recomputed with an independent numpy script: M* = 2.827e-8 kg = 1.299 m_P, c_kappa = 0.4136, eta(10 Msun) = 1.198e-78, P = 9.01e-31 W, T_H = 6.17e-9 K, A(10 Msun) = 4.97e-83, M_break = 1.401e-10 kg, A(M*) = 2.46e-5, t_evap = 2.1e70 yr, E = M* c^2 = 2.541e9 J, 11.78 bits, r_min and f(r_min) at 0.5, 0.9, 0.99 M*, light-ring threshold 0.809543.
* Direct integration (new code) agrees with the original time-domain signal: 0.42963 - 0.003072 i against 0.42965 and 3.081e-3 (0.3% in the decay rate).

## 6. Known limits

* Not re-run in the audit (modules missing): `check_symbolic.py`, `check_analytic.py`, `crosscheck_wkb_mpmath.py`; their shipped logs come from the author's machine (97/97 checks in `verification/REPORT_full.md`, which predates the edits of this change and must be regenerated). Not re-run: the `--full` checks T09 and Q14 (about 55 min of CPU).
* Loose tolerances in `verification/expected.json`: order-of-magnitude checks A24, A26, A27 (factor 3.2); ranges for qualitative statements A20, A28, A35-A39; 3-5% for A02, A13, A15, A21, A23, A25; 20% for W09, W10 ("Q ~ 20", "tau ~ 1e2"); 25% for T05; Q05-Q10 3-4%.
* The WKB values of Table 1 are leading-order and overestimate Q by factors 1.1-1.8 (direct integration, Sec. 5.4); they are not independent of the corrected prefactor: checks W01-W06 compare the table with the same code. Only P01-P05 compare with the original.
* The envelope-fit decay rate depends on a 10% edge trim (without it 2.87e-3); the least-squares fit (3.081e-3) is the independent confirmation. The omega_R digits beyond +-0.0003 are not significant.
* The corrected prefactor is a leading-order Gamow estimate; direct integration shows a residual stable offset of about 20% at 0.90 M* (l=4..7).
* A check that a rerun of `verify.sh` regenerates `ancillary/output/` and `verification/*` identically across machines is not claimed: numbers agree to the printed digits, files differ in version header lines.

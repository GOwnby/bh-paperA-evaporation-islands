"""Extract omega_R and the amplitude decay rate from the late-time signal written by timedomain.py
(window t in [500, 2200]; omega_R from the peak of a zero-padded Hann-windowed FFT, decay rate from a
straight-line fit of log(analytic-signal envelope) with 10% trimmed at each end).
Writes output/timedomain_fit.txt and .json.   Usage: python3 fit_timedomain.py   (needs timedomain.py output first;
optional grid-convergence files: python3 timedomain.py --nx 3000 ; python3 timedomain.py --nx 12000;
an explicit list of grid sizes may be given: python3 fit_timedomain.py 3000 6000 12000 24000)
Written after the original paper (v1.0.1): the original fit was done by hand and its code is not preserved.  The estimators
"decay"/"wR" above reproduce the v1.0.1 decay rate 3.08e-3; the independent damped-sinusoid least-squares fit
("decay_lsq", "wR_lsq": y = exp(-g (t-500)) [A cos(w t) + B sin(w t)], every 10th sample) was added as a check of the 10%-trim envelope
estimator and gives omega_R to ~1e-5 resolution (the FFT peak has bin width 2.5e-4).
"""
import os, sys, numpy as np
from scipy.optimize import least_squares
from _util import OUT, tee_stdout, save_json
from wkb import solve_mode, Mstar

tee_stdout("timedomain_fit.txt")
print("l=2, M=0.90 M*, observer x=25, window t in [500,2200]")
_w = solve_mode(0.90*Mstar, 2, 0); WI = abs(_w['wI'])
print(f"WKB (wkb.py): omega_R = {_w['wR']:.4f}, |Im omega| = {abs(_w['wI']):.3e}")
found = False
J = {"wkb": {"wR": _w["wR"], "absIm": WI}, "runs": {}}
NXS = [int(a) for a in sys.argv[1:]] or [3000, 6000, 12000]     # grid sizes to analyse (default: the three shipped-by-verify runs)
for nx in NXS:
    fn = os.path.join(OUT, f"td_l2_nx{nx}.npz")
    if not os.path.exists(fn):
        print(f"Nx={nx}: (file missing; run timedomain.py --nx {nx})"); continue
    found = True
    d = np.load(fn); t, s = d["t"], d["s"]; dt = t[1]-t[0]
    m = (t >= 500) & (t <= 2200); y = s[m]; tt = t[m]
    w = np.hanning(len(y)); N = 1 << 22
    F = np.abs(np.fft.rfft(y*w, N)); fr = np.fft.rfftfreq(N, dt)*2*np.pi
    wR = fr[F.argmax()]
    Y = np.fft.fft(y); h = np.zeros(len(y)); h[0] = 1; h[1:len(y)//2] = 2; h[len(y)//2] = 1
    env = np.abs(np.fft.ifft(Y*h)); inn = slice(len(y)//10, -len(y)//10)
    p = np.polyfit(tt[inn], np.log(env[inn]), 1)
    ts, ys = tt[::10], y[::10]
    resf = lambda a: np.exp(-a[3]*(ts-500.0))*(a[0]*np.cos(a[2]*ts) + a[1]*np.sin(a[2]*ts)) - ys
    sol = least_squares(resf, [ys.max(), 0.0, wR, -p[0]], x_scale=[ys.max(), ys.max(), 1e-3, 1e-3])
    print(f"Nx={nx}: omega_R = {wR:.4f}   decay rate = {-p[0]:.3e}   (decay / WKB = {-p[0]/WI:.2f})")
    print(f"        damped-sinusoid least squares: omega_R = {sol.x[2]:.5f}   decay rate = {sol.x[3]:.4e}")
    J["runs"][str(nx)] = {"wR": wR, "decay": -p[0], "decay_over_wkb": -p[0]/WI, "wR_lsq": sol.x[2], "decay_lsq": sol.x[3]}
if not found:
    print("no time-domain output found: run `python3 timedomain.py` first")

save_json("timedomain_fit.json", J)

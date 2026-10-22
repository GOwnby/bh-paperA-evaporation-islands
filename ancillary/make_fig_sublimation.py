"""Regenerates Fig. 1 (time-dilation factor Gamma = sqrt(f(r=2M)) = ell/sqrt(4M^2+ell^2) versus remaining mass).
Units: G=c=1, ell = 1e-3 in units of the largest plotted mass (M = 1 ... 1e-4); M* = (3 sqrt3/4) ell.
Writes output/sublimation_curve.png (same curve as the figure shipped in the repo root) and output/sublimation_curve.json
(sample values, the closed form Gamma = ell/sqrt(4M^2+ell^2), and the pixel comparison with ../sublimation_curve.png).
Usage: python3 make_fig_sublimation.py [--out PATH]    (about 2 s; needs matplotlib)
"""
import argparse, os
import numpy as np
import matplotlib
import matplotlib.image
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from _util import OUT, save_json

ap = argparse.ArgumentParser(); ap.add_argument("--out", default=os.path.join(OUT, "sublimation_curve.png"))
a = ap.parse_args()
ell = 1.0e-3
M = np.logspace(0, -4, 600)
Gamma = np.sqrt(1 - 2*M*(2*M)**2/((2*M)**3 + 2*M*ell**2))   # sqrt(f) at the naive Schwarzschild radius r=2M
Mstar = 3*np.sqrt(3)/4*ell
fig, ax = plt.subplots(figsize=(10, 6))
ax.semilogx(M, Gamma, color="crimson", lw=2.5, label=r"Time Dilation Factor ($\Gamma$)")
ax.axvline(Mstar, color="k", ls="--", label=r"Critical Mass ($M^*=\frac{3\sqrt{3}}{4}\ell$)")
ax.invert_xaxis(); ax.grid(True, which="both", ls="--", alpha=0.4); ax.legend()
ax.set_xlabel("Remaining Mass (Log Scale)"); ax.set_ylabel(r"Passage of Coordinate Time ($\Gamma=\frac{d\tau}{dt}$)")
ax.set_title(r"Relativistic Sublimation: The Unfreezing Rate Near $M^*$")
fig.tight_layout(); fig.savefig(a.out, dpi=100)
print(f"saved {os.path.relpath(a.out)}; Gamma(M*) = {ell/np.sqrt(4*Mstar**2+ell**2):.4f}, Gamma(M=1e-4) = {Gamma[-1]:.4f}, Gamma(M=1) = {Gamma[0]:.1e}")

# ---- checks recorded for verify.py ----
Gam_closed = ell/np.sqrt(4*M**2+ell**2)
J = {"Gamma_at_Mstar": float(ell/np.sqrt(4*Mstar**2+ell**2)), "Gamma_M1e-4": float(Gamma[-1]), "Gamma_M1": float(Gamma[0]),
     "max_abs_dev_from_closed_form": float(np.max(np.abs(Gamma-Gam_closed))), "Mstar": float(Mstar), "ell": ell,
     "monotone_increasing_as_M_decreases": bool(np.all(np.diff(Gamma) >= 0))}
ref = os.path.join(os.path.dirname(OUT), "..", "sublimation_curve.png")
if os.path.exists(ref):
    A = matplotlib.image.imread(a.out); B = matplotlib.image.imread(ref)
    J["shipped_png_same_shape"] = bool(A.shape[:2] == B.shape[:2])
    J["shipped_png_mean_abs_pixel_diff"] = float(np.mean(np.abs(A[..., :3]-B[..., :3]))) if A.shape[:2] == B.shape[:2] else None
save_json("sublimation_curve.json", J)

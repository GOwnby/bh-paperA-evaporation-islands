Scripts used for the trapped-mode spectrum of Section 5.4.

wkb.py         Three-turning-point WKB (Bohr-Sommerfeld + Gamow tunneling)
               for the effective potential V_l(r) = f[l(l+1)/r^2 + f'/r],
               units ell=1. Produces the (Re(w), Im(w), Q, tau) values in
               Table 1. Requires numpy, scipy.

timedomain.py  Independent cross-check: finite-difference evolution of the
               1+1 wave equation d^2 psi/dx^2 + (w^2 - V_l) psi = 0 in the
               time domain, for l=2, M=0.90 M*. Regularity (Dirichlet) at
               r=0, Sommerfeld outgoing condition at the outer boundary.
               Saves signal.npy, tarr.npy for frequency/decay extraction
               (done by hand in the paper; see the session transcript for
               the exact fitting code).

Run: python3 wkb.py ; python3 timedomain.py

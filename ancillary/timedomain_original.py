import numpy as np
from scipy.integrate import quad, cumulative_trapezoid
from scipy.interpolate import interp1d
import time

l_scale=1.0; Mstar=3*np.sqrt(3)/4*l_scale
M=0.90*Mstar
l=2

def f(r): return 1 - 2*M*r**2/(r**3+2*M*l_scale**2)
def fp(r,h=1e-6): return (f(r+h)-f(r-h))/(2*h)
def V(r): return f(r)*(l*(l+1)/r**2 + fp(r)/r)

# build r(x) via tortoise coordinate: x(r) = int_0^r dr'/f(r')
rgrid = np.linspace(1e-5, 80, 400000)
integrand = 1.0/f(rgrid)
xgrid = cumulative_trapezoid(integrand, rgrid, initial=0.0)
r_of_x = interp1d(xgrid, rgrid, kind='cubic', bounds_error=False, fill_value=(0.0, rgrid[-1]))
xmax = xgrid[-1]
print(f"x range: [0, {xmax:.2f}]  (r up to {rgrid[-1]:.1f})")

# spatial grid in x
Nx = 6000
xarr = np.linspace(0, min(xmax, 90), Nx)
dx = xarr[1]-xarr[0]
rarr = r_of_x(xarr)
Varr = V(np.clip(rarr,1e-6,None))
Varr[0] = Varr[1]*3 - Varr[2]*2   # extrapolate near the (regular, steep) origin rather than trust interp there

dt = 0.4*dx
Tfinal = 2200.0
Nt = int(Tfinal/dt)
print(f"Nx={Nx}, dx={dx:.4f}, dt={dt:.4f}, Nt={Nt}")

# initial data: Gaussian pulse in psi outside the barrier, psi_t=0
x0 = 20.0; sigma=2.0
psi = np.exp(-(xarr-x0)**2/(2*sigma**2))
psi[0]=0.0  # regularity/Dirichlet at the origin (x=0 <-> r=0), required for l>=1
psi_prev = psi.copy()  # psi_t=0 => psi(-dt)=psi(0) to O(dt^2) using the wave eq for the first step properly below

# proper first step using Taylor expansion with psi_tt = psi_xx - V psi, psi_t(0)=0
lap = np.zeros_like(psi)
lap[1:-1] = (psi[2:]-2*psi[1:-1]+psi[:-2])/dx**2
psi_tt0 = lap - Varr*psi
psi_next_est = psi + 0.5*dt**2*psi_tt0
psi_prev = psi.copy()
psi_cur  = psi_next_est.copy()
psi_cur[0]=0.0

# observer location: just outside the barrier
obs_x = 25.0
obs_i = np.argmin(np.abs(xarr-obs_x))
signal = []

t0=time.time()
for n in range(Nt):
    lap = np.zeros_like(psi_cur)
    lap[1:-1] = (psi_cur[2:]-2*psi_cur[1:-1]+psi_cur[:-2])/dx**2
    psi_new = 2*psi_cur - psi_prev + dt**2*(lap - Varr*psi_cur)
    psi_new[0]=0.0
    # simple outgoing (Sommerfeld) condition at outer boundary: d/dt psi + d/dx psi = 0
    psi_new[-1] = psi_cur[-1] - (dt/dx)*(psi_cur[-1]-psi_cur[-2])
    signal.append(psi_cur[obs_i])
    psi_prev, psi_cur = psi_cur, psi_new
print("evolution time:", time.time()-t0, "s")

signal=np.array(signal)
tarr = np.arange(Nt)*dt
np.save('signal.npy', signal); np.save('tarr.npy', tarr)
print("done, signal length", len(signal))

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

l_scale = 1.0
Mstar = 3*np.sqrt(3)/4*l_scale

def f_of(M):
    return lambda r: 1 - 2*M*r**2/(r**3+2*M*l_scale**2)
def fp_of(M,h=1e-6):
    f=f_of(M); return lambda r:(f(r+h)-f(r-h))/(2*h)
def V_of(M,l):
    f=f_of(M); fp=fp_of(M)
    return lambda r: f(r)*(l*(l+1)/r**2 + fp(r)/r)

def well_barrier(M,l):
    V=V_of(M,l)
    r=np.linspace(1e-4,40,300000); Vv=V(r); dV=np.gradient(Vv,r)
    idx=np.where(dV[:-1]*dV[1:]<0)[0]
    pts=[brentq(lambda x:(V(x+1e-5)-V(x-1e-5))/2e-5, r[i], r[i+1]) for i in idx]
    assert len(pts)==2
    rmin,rmax=pts
    return rmin,V(rmin),rmax,V(rmax)

def turning_points(V,E,rmin,rmax,rhi=60.0):
    # r1 in (0,rmin): V(r1)=E on the falling branch from +inf
    r1=brentq(lambda x:V(x)-E, 1e-6, rmin)
    r2=brentq(lambda x:V(x)-E, rmin, rmax)
    r3=brentq(lambda x:V(x)-E, rmax, rhi)
    return r1,r2,r3

def bohr_sommerfeld_residual(E,V,f,rmin,rmax,n):
    r1,r2,_=turning_points(V,E,rmin,rmax)
    integrand=lambda r: np.sqrt(max(E-V(r),0))/f(r)
    val,_=quad(integrand,r1,r2,limit=200,points=[r1,r2],epsabs=1e-9)
    return val-(n+0.5)*np.pi

def solve_mode(M,l,n=0):
    f=f_of(M); V=V_of(M,l)
    rmin,Vmin,rmax,Vmax=well_barrier(M,l)
    lo,hi=Vmin*1.0001,Vmax*0.9999
    # bracket: residual negative near Vmin (well too shallow for n+1/2 quantum), positive near Vmax
    Elo=bohr_sommerfeld_residual(lo,V,f,rmin,rmax,n)
    Ehi=bohr_sommerfeld_residual(hi,V,f,rmin,rmax,n)
    if Elo*Ehi>0:
        return None  # this (n,l) doesn't fit in the well
    E=brentq(lambda E: bohr_sommerfeld_residual(E,V,f,rmin,rmax,n), lo,hi, xtol=1e-10)
    wR=np.sqrt(E)
    r1,r2,r3=turning_points(V,E,rmin,rmax)
    # tunneling integral through the barrier
    integrand=lambda r: np.sqrt(max(V(r)-E,0))/f(r)
    barrier,_=quad(integrand,r2,r3,limit=200,points=[r2,r3],epsabs=1e-9)
    # period (round trip) in the well
    per_integrand=lambda r: 1.0/(f(r)*np.sqrt(max(E-V(r),1e-14)))
    T,_=quad(per_integrand,r1,r2,limit=200,points=[r1,r2],epsabs=1e-9)
    T*=2
    wI = -(1.0/T)*np.exp(-2*barrier)
    return dict(wR=wR, wI=wI, T=T, barrier_action=barrier, r1=r1,r2=r2,r3=r3, Vmin=Vmin,Vmax=Vmax)

print(f"{'M/M*':>6} {'l':>3} {'n':>2}  {'Re(w) l':>10} {'Im(w) l':>12} {'Q=Re/2|Im|':>12} {'tau=1/|Im|':>12} {'barrier action':>15}")
rows=[]
for Mfrac in [0.90,0.95,0.99]:
    M=Mfrac*Mstar
    for l in [5,10,20,30]:
        for n in [0,1]:
            res=solve_mode(M,l,n)
            if res is None:
                print(f"{Mfrac:6.2f} {l:3d} {n:2d}   (no mode fits: well too shallow for this n)")
                continue
            Q = res['wR']/(2*abs(res['wI']))
            tau = 1.0/abs(res['wI'])
            print(f"{Mfrac:6.2f} {l:3d} {n:2d}  {res['wR']:10.5f} {res['wI']:12.3e} {Q:12.4e} {tau:12.4e} {res['barrier_action']:15.4f}")
            rows.append((Mfrac,l,n,res))
import pickle
pickle.dump(rows, open('rows.pkl','wb'))

print("\n--- low-l check (l=1,2,3,4), M=0.90 M* ---")
M=0.90*Mstar
for l in [1,2,3,4]:
    res=solve_mode(M,l,0)
    if res is None:
        print(f"l={l}: no n=0 mode fits in the well")
    else:
        Q=res['wR']/(2*abs(res['wI'])); tau=1/abs(res['wI'])
        print(f"l={l}: Re(w)={res['wR']:.5f}  Im(w)={res['wI']:.3e}  Q={Q:.3e}  tau={tau:.3e} (units of l_Planck/c)  barrier_action={res['barrier_action']:.3f}")

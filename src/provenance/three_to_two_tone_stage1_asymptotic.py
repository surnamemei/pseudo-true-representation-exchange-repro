"""Finite-N small-spacing tangent model for the 3-to-2 tone crossing."""

import numpy as np
from scipy.optimize import brentq, minimize_scalar


def tangent_cost(u, lam, n=21, b=10.0, phase=np.pi):
    s=(np.arange(n)-(n-1)/2)/n
    z=np.exp(1j*u*s)
    target=-s*s+lam*np.exp(1j*phase)*np.exp(1j*b*s)
    # Base amplitude is complex; central-frequency displacement is real;
    # satellite amplitude is complex. Frequency shift is 2*kappa*d^2/N.
    design=np.column_stack((np.ones(n),1j*np.ones(n),
                            2j*s,z,1j*z))
    real=np.vstack((design.real,design.imag))
    yy=np.r_[target.real,target.imag]
    coef=np.linalg.lstsq(real,yy,rcond=None)[0]
    return float(np.linalg.norm(yy-real@coef)**2)


def tangent_branch(lam, seed, n=21, b=10.0, phase=np.pi):
    bounds=(-8.0,-1.0) if seed<0 else (8.0,17.0)
    res=minimize_scalar(lambda u:tangent_cost(u,lam,n,b,phase),
                        bounds=bounds,method='bounded',
                        options={'xatol':1e-10})
    return res.x,res.fun


def tangent_crossing(n=21,b=10.0,phase=np.pi):
    def diff(lam):
        return tangent_branch(lam,-4.17,n,b,phase)[1]-tangent_branch(lam,12.47,n,b,phase)[1]
    ll=np.linspace(.04,.1,50)
    vv=[diff(v) for v in ll]
    for i in range(len(ll)-1):
        if vv[i]*vv[i+1]<0:
            root=brentq(diff,ll[i],ll[i+1],xtol=1e-12)
            return root,tangent_branch(root,-4.17,n,b,phase),tangent_branch(root,12.47,n,b,phase)
    return None


if __name__=='__main__':
    for n in (17,21,25,31,41,61):
        print(n,tangent_crossing(n),flush=True)

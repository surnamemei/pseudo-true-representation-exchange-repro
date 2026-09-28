"""Explanatory N=21 tangent landscape; plot is visual, not a certificate."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

ROOT=Path(__file__).resolve().parent
N=21
s=np.arange(-10,11,dtype=float)/N
S2=np.sum(s*s);S4=np.sum(s**4)
def D(v):return np.cos(np.multiply.outer(np.asarray(v),s)).sum(axis=-1)
def D1(v):return -(np.sin(np.multiply.outer(np.asarray(v),s))*s).sum(axis=-1)
def D2(v):return -(np.cos(np.multiply.outer(np.asarray(v),s))*s*s).sum(axis=-1)
def R(v,lam):
    v=np.asarray(v)
    q0=-S2-lam*D(10)
    q1=2*lam*D1(10)
    C=S4-2*lam*D2(10)+N*lam*lam-q0*q0/N-q1*q1/(4*S2)
    a=N-D(v)**2/N-D1(v)**2/S2
    b=D2(v)-lam*D(10-v)-D(v)*q0/N+D1(v)*q1/(2*S2)
    with np.errstate(divide='ignore',invalid='ignore'):
        out=C-b*b/a
    return np.where(np.abs(v)<0.15,np.nan,out)

lc=float(json.loads((ROOT/'stage3_tangent_reference.json').read_text())['lambda_N'])
lamvals=[lc-.002,lc,lc+.002]
grid=np.linspace(-7,16,3000)
fig,ax=plt.subplots(1,3,figsize=(9.4,2.75),sharex=True,sharey=True)
for k,(a,lam) in enumerate(zip(ax,lamvals)):
    mA=minimize_scalar(lambda v:float(R(v,lam)),bounds=(-5,-3),method='bounded')
    mB=minimize_scalar(lambda v:float(R(v,lam)),bounds=(11,14),method='bounded')
    baseline=min(mA.fun,mB.fun)
    y=R(grid,lam)-baseline
    a.plot(grid,y,color='#155e75',lw=1.7)
    a.plot([mA.x,mB.x],[mA.fun-baseline,mB.fun-baseline],'o',ms=5,
           color='#be123c')
    a.text(mA.x-.35,mA.fun-baseline+.005,'A',fontsize=9)
    a.text(mB.x-.15,mB.fun-baseline+.005,'B',fontsize=9)
    a.set_title([r'$\lambda_c-0.002$',r'$\lambda_c$',r'$\lambda_c+0.002$'][k],fontsize=10)
    a.set_xlabel(r'Satellite coordinate $v$',fontsize=9)
    a.grid(alpha=.22)
    a.set_xlim(-7,16);a.set_ylim(-.001,.07)
ax[0].set_ylabel(r'$R_{21}(v,\lambda)-\min_v R_{21}$',fontsize=9)
fig.tight_layout(pad=.6)
path=ROOT/'paper'/'figures'/'fig_stage17_tangent_landscape.pdf'
fig.savefig(path)
print(path)

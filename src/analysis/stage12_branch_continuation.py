import json,csv
from pathlib import Path
import numpy as np
from scipy.optimize import root,minimize_scalar
from mpmath import mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from three_to_two_tone_stage3_full_crossings import value_grad
R=Path(__file__).resolve().parent
ref=json.loads((R/'stage2_reference.json').read_text());ec=float(ref['epsilon_cross']);s=np.arange(-10,11)/21

def D(q,k=0):return np.sum((1j*s)**k*np.exp(1j*q*s)).real

def calc(u,e):
 u1,u2=u;g=D(u2-u1);gp=D(u2-u1,1);gpp=D(u2-u1,2)
 h=np.array([D(v-2)+D(v+2)-e*D(v-10) for v in u]);hp=np.array([D(v-2,1)+D(v+2,1)-e*D(v-10,1) for v in u]);hpp=np.array([D(v-2,2)+D(v+2,2)-e*D(v-10,2) for v in u])
 G=np.array([[21,g],[g,21]]);a=np.linalg.solve(G,h);E=42+2*D(4)+21*e*e-2*e*(D(8)+D(12));J=E-h@a
 gpairs=[np.array([[0,-gp],[-gp,0]]),np.array([[0,gp],[gp,0]])]
 grad=np.array([-2*hp[k]*a[k]+a@gpairs[k]@a for k in range(2)])
 deriv=[]
 for k in range(2):
  q=np.zeros(2);q[k]=hp[k];deriv.append(np.linalg.solve(G,q-gpairs[k]@a))
 H=np.zeros((2,2))
 for i in range(2):
  for j in range(2):
   gg=np.array([[0,gpp*(1 if i==j else -1)],[gpp*(1 if i==j else -1),0]])
   H[i,j]=-2*((hpp[i]*a[i] if i==j else 0)+hp[i]*deriv[j][i])+2*deriv[j]@gpairs[i]@a+a@gg@a
 return J,grad,H,a

def track(branch,es,seed):
 rows=[];stop=None
 for e in es:
  sol=root(lambda u:calc(u,e)[1],seed,jac=lambda u:calc(u,e)[2],tol=1e-11)
  u=sol.x;J,g,H,a=calc(u,e);eig=np.linalg.eigvalsh(H)
  if np.linalg.norm(g)>1e-8 or eig[0]<=1e-9 or np.linalg.norm(u-seed)>.5:
   stop=dict(epsilon=float(e),u=list(u),gradient=list(g),hessian_eigs=list(eig),previous=list(seed));break
  row=dict(branch=branch,epsilon=float(e),u1=u[0],u2=u[1],amplitude1_re=a[0],amplitude1_im=0.,amplitude2_re=a[1],amplitude2_im=0.,J=J,Hessian_min=eig[0],Hessian_max=eig[1],gradient_norm=np.linalg.norm(g))
  rows.append(row);seed=u
 return rows,stop

def main():
 A,sa=track('A',np.unique(np.r_[np.linspace(0,.30,601),.05,.1,ec]),np.array([-2.,2.]))
 B,sb=track('B',np.unique(np.r_[np.arange(0,ec,.0005),ec])[::-1],np.array(list(map(float,ref['B_at_crossing']))))
 Bup,_=track('B',np.arange(ec+.0005,.3001,.0005),np.array(list(map(float,ref['B_at_crossing']))));B+=Bup
 rows=sorted(A+B,key=lambda r:(r['branch'],r['epsilon']))
 with (R/'branch_continuation.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 mp.dps=65;high=[]
 for br in ('A','B'):
  for e in ('0','0.05','0.1',ref['epsilon_cross']):
   candidates=[r for r in rows if r['branch']==br];near=min(candidates,key=lambda r:abs(r['epsilon']-float(e)))
   if abs(near['epsilon']-float(e))>.00051:continue
   try:
    u=mp.findroot(lambda a,b:value_grad(a,b,mp.mpf(e),mp.mpf(2))[1:],(mp.mpf(str(near['u1'])),mp.mpf(str(near['u2']))),tol=mp.mpf('1e-55'))
    high.append(dict(branch=br,epsilon=e,u1=mp.nstr(u[0],55),u2=mp.nstr(u[1],55),J=mp.nstr(value_grad(*u,mp.mpf(e),mp.mpf(2))[0],55)))
   except Exception as ex:high.append(dict(branch=br,epsilon=e,error=str(ex)))
 fixed=minimize_scalar(lambda v:calc([v,10],ec)[0],bounds=(-2,2),method='bounded',options={'xatol':1e-13})
 report=dict(epsilon_cross=ref['epsilon_cross'],A_stop=sa,B_backward_stop=sb,A_points=len(A),B_points=len(B),high_precision_checks=high,satellite_fixed_b=dict(u1=float(fixed.x),u2=10,J=float(fixed.fun)),crossing_cost=float(next(r for r in high if r['branch']=='A' and r['epsilon']==ref['epsilon_cross'])['J']))
 (R/'stage12_continuation_checks.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
 fig,ax=plt.subplots(figsize=(7.1,2.8))
 for br,color in [('A','#185fa5'),('B','#b65318')]:
  rs=sorted([r for r in rows if r['branch']==br],key=lambda r:r['epsilon'])
  for k,style in [(1,'-'),(2,'--')]:ax.plot([r['epsilon'] for r in rs],[r['u'+str(k)] for r in rs],style,color=color,lw=1.5,label=br+str(k))
 for val in [-2,2,10]:ax.axhline(val,color='.65',lw=.75,ls=':');ax.text(.303,val,f'{val:g}',va='center',color='.4',fontsize=8)
 ax.axvline(ec,color='black',ls='-.',lw=1);ax.text(ec-.004,14.4,r'$\epsilon_c$',ha='right',fontsize=9)
 ax.set(xlabel=r'Omitted-tone amplitude $\epsilon$',ylabel=r'Fitted frequency $u=N\nu$',xlim=(0,.30));ax.grid(alpha=.15);ax.legend(ncol=4,fontsize=8,loc='center',bbox_to_anchor=(.5,.5));fig.tight_layout()
 fig.savefig(R/'branch_continuation_figure.pdf');fig.savefig(R/'branch_continuation_figure.png',dpi=180);fig.savefig(R/'paper/figures/fig1_continuation.pdf');plt.close(fig)
if __name__=='__main__':main()

"""Numerical weighted two-tone LS landscape; all frequencies remain off grid."""
import csv,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize,brentq
from scipy.ndimage import minimum_filter
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
N=21;S=(np.arange(N)-10)/N;Z=2.;B=10.
EC=.2481906301722774
SEEDS=[[-4.772447353918614,.605170320276862],[-.057653561926680,12.46341935680166],[-2,2]]
def signal(e):return 2*np.cos(Z*S)-e*np.exp(1j*B*S)
def canonical(u):return np.sort((np.asarray(u)+np.pi*N)%(2*np.pi*N)-np.pi*N)
def objective(u,y,w,detail=False):
    v=np.exp(1j*np.outer(S,u));sw=np.sqrt(w)
    a=np.linalg.lstsq(sw[:,None]*v,sw*y,rcond=None)[0];r=y-v@a
    val=np.vdot(r,w*r).real
    g=-2*np.real(np.sum(np.conj((w*r)[:,None])*(1j*S[:,None]*v)*a,axis=0))
    return (val,g,a) if detail else (val,g)
def fit(seed,y,w):
    r=minimize(lambda u:objective(u,y,w),seed,jac=True,method='BFGS',options={'gtol':1e-10,'maxiter':220})
    u=canonical(r.x);j,g,a=objective(u,y,w,True)
    h=np.column_stack([(objective(u+np.eye(2)[k]*1e-4,y,w)[1]-objective(u-np.eye(2)[k]*1e-4,y,w)[1])/2e-4 for k in range(2)])
    return dict(u=u,J=float(j),gradient=float(np.linalg.norm(g)),curvature=float(np.linalg.eigvalsh((h+h.T)/2)[0]),amplitudes=a)
class Grid:
    def __init__(self,w,m=384):
        self.w=w;self.g=np.linspace(-np.pi*N,np.pi*N,m,endpoint=False)
        self.v=np.exp(1j*np.outer(S,self.g));self.gram=(self.v.conj().T@(w[:,None]*self.v)).real
        self.den=np.sum(w)**2-self.gram**2
    def seeds(self,y,top=45):
        h=self.v.conj().T@(self.w*y);nn=np.sum(self.w)
        cap=(nn*(abs(h[:,None])**2+abs(h[None,:])**2)-2*self.gram*np.real(np.conj(h[:,None])*h[None,:]))/np.maximum(self.den,1e-15)
        vals=np.vdot(y,self.w*y).real-cap;vals[np.tril_indices(len(h),0)]=np.inf
        ii=np.argwhere(np.isfinite(vals)&(vals<=minimum_filter(vals,size=3,mode='nearest')+1e-10))
        ii=sorted(ii,key=lambda ij:vals[tuple(ij)])[:top]
        return [self.g[ij] for ij in ii]
def search(e,w,grid,extra=()):
    y=signal(e);rng=np.random.default_rng(1207)
    seeds=grid.seeds(y)+SEEDS+list(extra)+list(rng.uniform(-np.pi*N,np.pi*N,(12,2)))
    roots=[]
    for seed in seeds:
        f=fit(seed,y,w)
        if f['curvature']>1e-9 and f['gradient']<1e-6 and np.diff(f['u'])[0]>.001 and not any(np.linalg.norm(f['u']-q['u'])<.02 for q in roots):roots.append(f)
    return sorted(roots,key=lambda q:q['J'])
def write(name,rows):
    with (ROOT/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
def run():
    start=time.perf_counter();out=[];summary=[];plots=[]
    for name,raw in [('rectangular',np.ones(N)),('Hann',np.hanning(N)),('Hamming',np.hamming(N))]:
        # Window denotes LS precision weight, J_w=sum w_n |r_n|^2.
        w=N*raw/raw.sum();grid=Grid(w);etas=np.linspace(.05,.8,31);prev=[];best=[]
        for e in etas:
            rr=search(e,w,grid,prev);prev=[q['u'] for q in rr[:2]];best.append((e,rr))
            for rank,q in enumerate(rr[:3]):out.append(dict(window=name,epsilon=e,rank=rank+1,u1=q['u'][0],u2=q['u'][1],J=q['J'],curvature=q['curvature'],gradient=q['gradient'],n_minima=len(rr)))
        candidates=[]
        for (e0,r0),(e1,r1) in zip(best,best[1:]):
            if np.linalg.norm(r0[0]['u']-r1[0]['u'])>1 and len(r0)>1 and len(r1)>1:
                a,b=r0[0]['u'],r1[0]['u']
                def gap(e):return fit(a,signal(e),w)['J']-fit(b,signal(e),w)['J']
                if gap(e0)*gap(e1)<0:
                    ec=brentq(gap,e0,e1,xtol=1e-12);rr=search(ec,w,Grid(w,768),[a,b]);aa=fit(a,signal(ec),w);bb=fit(b,signal(ec),w)
                    distance=np.linalg.norm(aa['u']-bb['u'])
                    if distance<=1 or min(aa['curvature'],bb['curvature'])<=0:continue
                    if max(aa['J'],bb['J'])>rr[0]['J']+1e-7:continue
                    slope=(gap(ec+1e-5)-gap(ec-1e-5))/2e-5
                    if abs(slope)<1e-4:continue
                    candidates.append(dict(window=name,crossing_found=True,epsilon_crossing=ec,uA1=aa['u'][0],uA2=aa['u'][1],uB1=bb['u'][0],uB2=bb['u'][1],cost=(aa['J']+bb['J'])/2,branch_distance=distance,min_curvature=min(aa['curvature'],bb['curvature']),third_margin=rr[2]['J']-rr[0]['J'] if len(rr)>2 else float('nan'),n_minima=len(rr),slope=slope,evidence='numerical_only'))
        if not candidates:candidates=[dict(window=name,crossing_found=False,epsilon_crossing='',uA1='',uA2='',uB1='',uB2='',cost='',branch_distance='',min_curvature='',third_margin='',n_minima='',slope='',evidence='none_detected_on_0.05_to_0.8')]
        summary+=candidates;plots.append((name,best));print(name,candidates,flush=True)
    write('window_robustness_minima.csv',out);write('window_robustness_summary.csv',summary)
    fig,axes=plt.subplots(1,3,figsize=(9,2.8),sharey=True)
    for ax,(name,points) in zip(axes,plots):
        for k in range(2):ax.plot([e for e,_ in points],[r[0]['u'][k] for _,r in points],'.-',ms=3,lw=1,label=rf'$u_{k+1}$')
        for c in summary:
            if c['window']==name and c['crossing_found']:ax.axvline(c['epsilon_crossing'],ls='--',color='grey',lw=1)
        ax.set_title(name);ax.set_xlabel(r'Omitted amplitude $\epsilon$');ax.grid(alpha=.2)
    axes[0].set_ylabel('Best numerical fitted frequencies');axes[0].legend(frameon=False,fontsize=8)
    fig.tight_layout();fig.savefig(ROOT/'window_robustness.pdf');plt.close(fig)
    (ROOT/'window_robustness_run.json').write_text(json.dumps(dict(seconds=time.perf_counter()-start,grid=384,crossing_grid=768,random_starts=12,grid_minimum_seeds=45,weight_convention='sum(w)=N; J=sum w|residual|^2; symmetric Hann/Hamming as weights'),indent=2))
if __name__=='__main__':run()

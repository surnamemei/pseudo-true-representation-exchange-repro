import time
import numpy as np
from three_to_two_tone_stage5_core import *

n=21; s=times(n); e=.2481906301722774; x=signal(n,2,10,e)
rg=np.random.default_rng(20260926)
for snr in (20,30,40):
    sig=(np.vdot(x,x).real/n*10**(-snr/10))**.5
    st=time.time(); mism=0; other=0; ngood=0; maxgap=0
    for i in range(40):
        y=x+sig*(rg.normal(size=n)+1j*rg.normal(size=n))/2**.5
        a=refine([-4.772447,.60517],y,s)
        b=refine([-.057654,12.463419],y,s)
        gr=grid_minima(y,s,192,24)
        best=gr[0]
        g=min(a['J'],b['J'])
        mism+=int(best['J']<g-1e-6)
        other+=int(all(np.linalg.norm(best['u']-f['u'])>.05 for f in (a,b)))
        ngood+=int(a['converged'] and b['converged'])
        maxgap=max(maxgap,g-best['J'])
    print(snr,mism,other,ngood,maxgap,time.time()-st,flush=True)

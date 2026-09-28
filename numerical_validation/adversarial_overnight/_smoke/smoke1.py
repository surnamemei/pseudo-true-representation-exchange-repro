import os, sys, time
os.environ.setdefault("MKL_NUM_THREADS","1"); os.environ.setdefault("OMP_NUM_THREADS","1")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from advrun import run_configs, shutdown_pool
if __name__=='__main__':
    cfgs=[dict(cid='smoke_N15_rect',stage=99,N=15,z=2.0,b=10.0,phi=3.141592653589793,delta=0.0,amp_minus=1.0,amp_plus=1.0,window='rect',scale='abs',lo=0.0,hi=1.0)]
    t=time.time()
    out=run_configs(99,os.path.join(os.path.dirname(os.path.abspath(__file__)),'out'),cfgs,desc='smoke')
    print('seconds',time.time()-t)
    a=out['smoke_N15_rect']
    print(a['verdict'], a['seconds'])
    for s in a['switches']:
        print({k:v for k,v in s.items() if k not in ('P','Q','third','P_lo','Q_hi','validation_lower_other')})
    shutdown_pool()

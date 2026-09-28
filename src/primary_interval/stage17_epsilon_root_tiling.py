"""Krawczyk continuation of A/B roots across epsilon_c +/-0.0025."""
import json
import sys
from pathlib import Path
from mpmath import mp,iv
from three_to_two_tone_stage3_full_crossings import value_grad
import stage7_n31_certificate as lc
import stage17_parametric_jet as p
import three_to_two_tone_stage2_interval as it

ROOT=Path(__file__).resolve().parent
mp.dps=85;iv.dps=80
it.N=21
lc.interval.N=21
lc.interval.objective_jet=p.objective
ref=json.loads((ROOT/'stage2_reference.json').read_text())
bracket=json.loads((ROOT/'stage2_crossing_global.json').read_text())['epsilon_bracket']
ec=(mp.mpf(bracket[0])+mp.mpf(bracket[1]))/2
radius=mp.mpf(sys.argv[1]) if len(sys.argv)>1 else mp.mpf('.0025')
lo=ec-radius;hi=ec+radius
bins=int(sys.argv[2]) if len(sys.argv)>2 else 256
records=[]
seed={w:tuple(map(mp.mpf,ref[w+'_at_crossing'])) for w in ('A','B')}
all_pass=True
for i in range(bins):
    a=lo+(hi-lo)*i/bins;b=lo+(hi-lo)*(i+1)/bins;m=(a+b)/2
    row=dict(index=i,epsilon_interval=[str(a),str(b)],branches={})
    for w,rad in [('A','0.0003'),('B','0.01')]:
        x=mp.findroot(lambda u,v:value_grad(u,v,m,2)[1:],seed[w],tol=mp.mpf('1e-60'),maxsteps=40)
        seed[w]=tuple(x)
        cert,_,X=lc.local_check([str(v) for v in x],[str(a),str(b)],rad)
        g=p.coefficients(*X)
        ep=iv.mpf([str(a),str(b)])
        slope_part=2*g[0].v*ep+g[1].v
        row['branches'][w]=dict(center=[str(v) for v in x],
                                X=[[str(v.a),str(v.b)] for v in X],
                                inclusion=cert['inclusion'],
                                positive_Hessian=cert['positive_Hessian'],
                                contraction_upper=cert['contraction_upper'],
                                dJ_de=[str(slope_part.a),str(slope_part.b)])
        if not (cert['inclusion'] and cert['positive_Hessian'] and cert['contraction_lt_one']):
            all_pass=False
    row['slope_lower']=str((iv.mpf(row['branches']['A']['dJ_de'][0])-
                            iv.mpf(row['branches']['B']['dJ_de'][1])).a)
    if not iv.mpf(row['slope_lower']).a>0:all_pass=False
    records.append(row)
    if (i+1)%32==0:print('root slabs',i+1,'pass',all_pass,flush=True)
out=dict(interval=[str(lo),str(hi)],slabs=bins,all_pass=all_pass,
         minimum_slope_lower=str(min((iv.mpf(r['slope_lower']).a for r in records))),
         records=records)
suffix='' if radius==mp.mpf('.0025') else '_'+str(radius).replace('.','p')
(ROOT/f'stage17_epsilon_root_tiling{suffix}.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='records'},indent=2))

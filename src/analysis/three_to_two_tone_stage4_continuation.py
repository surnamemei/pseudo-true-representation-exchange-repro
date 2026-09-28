"""High-precision numerical continuation of the two crossing branches.

Rows are numerical until independently interval validated. No global claim is
inferred from a smooth curve or from a finite set of rows.
"""
import csv, json
from pathlib import Path
from mpmath import mp
from three_to_two_tone_stage3_full_crossings import crossing, value_grad, D

ROOT=Path(__file__).resolve().parent
N=21;B=10
mp.dps=75

def statistics(u1,u2,e,d):
    J=value_grad(u1,u2,e,d)[0]
    def gg(a,c):return value_grad(a,c,e,d)
    h11=mp.diff(lambda a:gg(a,u2)[1],u1)
    h12=mp.diff(lambda c:gg(u1,c)[1],u2)
    h22=mp.diff(lambda c:gg(u1,c)[2],u2)
    rad=mp.sqrt((h11-h22)**2+4*h12*h12)
    emin=(h11+h22-rad)/2;emax=(h11+h22+rad)/2
    g=D(u2-u1,N)
    h1=D(u1-d,N)+D(u1+d,N)-e*D(u1-B,N)
    h2=D(u2-d,N)+D(u2+d,N)-e*D(u2-B,N)
    den=N*N-g*g
    amp1=(N*h1-g*h2)/den; amp2=(N*h2-g*h1)/den
    de=mp.diff(lambda ee:value_grad(u1,u2,ee,d)[0],e)
    return dict(J=J,hmin=emin,hmax=emax,sep_u=u2-u1,
                amp1=amp1,amp2=amp2,dJ_de=de)

def row(d,seed=None):
    z=crossing(str(d),seed=seed)
    e=mp.mpf(z['epsilon_cross'])
    a1,a2,b1,b2=[mp.mpf(z[k]) for k in ('u_A1','u_A2','u_B1','u_B2')]
    A=statistics(a1,a2,e,mp.mpf(d));C=statistics(b1,b2,e,mp.mpf(d))
    out=dict(z_NDelta=mp.nstr(d,12),epsilon_cross=mp.nstr(e,60),
             u_A1=z['u_A1'],u_A2=z['u_A2'],u_B1=z['u_B1'],u_B2=z['u_B2'],
             J_A=mp.nstr(A['J'],45),J_B=mp.nstr(C['J'],45),
             slope_deltaJ_epsilon=mp.nstr(A['dJ_de']-C['dJ_de'],45),
             A_hessian_min=mp.nstr(A['hmin'],35),A_hessian_max=mp.nstr(A['hmax'],35),
             B_hessian_min=mp.nstr(C['hmin'],35),B_hessian_max=mp.nstr(C['hmax'],35),
             A_sep_u=mp.nstr(A['sep_u'],35),B_sep_u=mp.nstr(C['sep_u'],35),
             A_amp1=mp.nstr(A['amp1'],35),A_amp2=mp.nstr(A['amp2'],35),
             B_amp1=mp.nstr(C['amp1'],35),B_amp2=mp.nstr(C['amp2'],35),
             status='numerical_branch_crossing')
    return out,(a1,a2,b1,b2,e)

if __name__=='__main__':
    ds=[mp.mpf(k)/20 for k in range(2,41)]
    rows=[];history=[];fallbacks=0
    for d in ds:
        predictor=None
        if len(history)>=2:
            d0,x0=history[-2];d1,x1=history[-1]
            predictor=tuple(x1[i]+(x1[i]-x0[i])*(d-d1)/(d1-d0)
                            for i in range(5))
        used_fallback=False
        try:z,x=row(d,predictor)
        except (ValueError,ZeroDivisionError):
            z,x=row(d,None);fallbacks+=1;used_fallback=True
        z['predictor_corrector_used']=bool(predictor is not None)
        z['tangent_fallback']=used_fallback
        history.append((d,x))
        rows.append(z)
        print(str(d),z['epsilon_cross'][:22],z['A_hessian_min'][:13],
              z['B_hessian_min'][:13],flush=True)
    with (ROOT/'continuation_branches.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print('predictor-corrector fallbacks',fallbacks,flush=True)

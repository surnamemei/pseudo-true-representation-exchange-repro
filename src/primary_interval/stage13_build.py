"""Stage-13 evidence and figures. Tangent continuation is numerical; endpoint identities are exact."""
from __future__ import annotations
import csv
import json
from pathlib import Path
import sys
from mpmath import mp, iv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent
mp.dps=75
iv.dps=65
N=21
s=[mp.mpf(t)/N for t in range(-10,11)]

def D(x,k=0):
    return mp.fsum((1j*t)**k*mp.exp(1j*x*t) for t in s).real

S2=mp.fsum(t*t for t in s)
S4=mp.fsum(t**4 for t in s)
S6=mp.fsum(t**6 for t in s)
V2=S4-S2*S2/N
V3=S6-S4*S4/S2
M=D(10,3)+S4/S2*D(10,1)
Q=D(10,2)+S2*D(10)/N
LAM_G=V2/Q
LAM_AB=mp.mpf('0.0659804111269913661716')

def Cbase(l):
    q0=-S2-l*D(10);q1=2*l*D(10,1)
    return S4-2*l*D(10,2)+N*l*l-q0*q0/N-q1*q1/(4*S2)

def cost(v,l):
    if abs(v)<mp.mpf('1e-18'):
        return coal(l)
    d=D(v);dp=D(v,1)
    A=N-d*d/N-dp*dp/S2
    q0=-S2-l*D(10);q1=2*l*D(10,1)
    B=D(v,2)-l*D(10-v)-d*q0/N+dp*q1/(2*S2)
    return Cbase(l)-B*B/A

def coal(l):
    q0=-S2-l*D(10)
    q2=-S4+l*D(10,2)-S2*q0/N
    return Cbase(l)-q2*q2/V2

def strong(l):
    return l*l*(N-D(10)**2/N-D(10,1)**2/S2)

def root(l,seed):
    f=lambda x:mp.diff(lambda y:cost(y,l),x)
    x=mp.findroot(f,(mp.mpf(seed)*mp.mpf('.999'),mp.mpf(seed)*mp.mpf('1.001')),
                  solver='secant',tol=mp.mpf('1e-48'),maxsteps=45)
    return x

def ivD(x,k=0):
    ss=[iv.mpf(t)/N for t in range(-10,11)]
    if k%2==0:
        return sum((-1)**(k//2)*t**k*iv.cos(x*t) for t in ss)
    return sum((-1)**((k+1)//2)*t**k*iv.sin(x*t) for t in ss)

def ibounds(x):
    return [str(mp.mpf(x.a)),str(mp.mpf(x.b))]

def endpoint_certificate():
    si=[iv.mpf(t)/N for t in range(-10,11)]
    iS2=sum(t*t for t in si);iS4=sum(t**4 for t in si);iS6=sum(t**6 for t in si)
    iV2=iS4-iS2*iS2/N;iV3=iS6-iS4*iS4/iS2
    iM=ivD(10,3)+iS4/iS2*ivD(10,1)
    iQ=ivD(10,2)+iS2*ivD(10)/N
    ilam=iV2/iQ
    iDb=ivD(10);iDpb=ivD(10,1);iD2b=ivD(10,2)
    iKcoal=N-iDb*iDb/N-iDpb*iDpb/iS2-iQ*iQ/iV2
    iA=N-iDb*iDb/N-iDpb*iDpb/iS2
    iB=iD2b+iDb*iS2/N
    iRb=iS4-iS2*iS2/N-iB*iB/iA
    ithresh=iv.sqrt(iRb/iKcoal)
    cc=coal(LAM_AB);ab=cost(root(LAM_AB,mp.mpf('-4.15344048664877')),LAM_AB)
    cert={
      'scope':'analytic coalescent exclusion for all positive lambda; local near-zero branch emergence at lambda=0',
      'N':21,'b':10,'phase':'pi',
      'V2_interval':ibounds(iV2),'V3_interval':ibounds(iV3),
      'M_interval':ibounds(iM),'Q_interval':ibounds(iQ),
      'lambda_gamma_interval':ibounds(ilam),
      'Ccoal_lambda_squared_coefficient_interval':ibounds(iKcoal),
      'R_at_b_constant_interval':ibounds(iRb),
      'lambda_fixed_b_beats_coalescent_above_interval':ibounds(ithresh),
      'vA_linear_coefficient_interval':ibounds(3*iM/iV3),
      'coalescent_penalty_quadratic_coefficient_interval':ibounds(iM*iM/iV3),
      'lambda_AB_reference':mp.nstr(LAM_AB,25),
      'coalescent_at_AB':mp.nstr(cc,22),
      'AB_cost_at_AB_approx':mp.nstr(ab,22),
      'endpoint_gap_at_AB_approx':mp.nstr(cc-ab,22),
      'lambda_gamma_is_above_AB':bool(ilam.a>iv.mpf(str(LAM_AB)).b),
      'lower_positive_coalescent_crossing_excluded_on_0_to_AB':bool(ilam.a>iv.mpf(str(LAM_AB)).b and iM.b<0 and iV3.a>0),
      'coalescent_excluded_for_all_positive_lambda':bool(ithresh.b<ilam.a and iM.b<0 and iV3.a>0),
      'full_regular_global_continuation_certified':False,
      'note':'Directed interval constants prove nonzero endpoint derivative below lambda_gamma. They do not certify global identity of regular branch A on the whole interval.'
    }
    (ROOT/'lower_transition_certificate.json').write_text(json.dumps(cert,indent=2),encoding='utf-8')
    return cert

def main():
    cert=endpoint_certificate()
    vals=['0','0.00001','0.0001','0.001','0.002','0.005','0.01','0.015','0.02','0.025',
          '0.03','0.035','0.04','0.045','0.05','0.055','0.06','0.065',str(LAM_AB),
          '0.07','0.075','0.08','0.09','0.1','0.12','0.15','0.2']
    rows=[];a=mp.mpf('-0.1');b=mp.mpf('14.59')
    for z in vals:
        l=mp.mpf(z)
        if l==0:
            va=mp.mpf(0);ra=mp.mpf(0)
            vb=root(l,b);rb=cost(vb,l)
        else:
            if l<mp.mpf('.001'):a=3*M*l/V3
            va=root(l,a);ra=cost(va,l)
            vb=root(l,b);rb=cost(vb,l)
            a=va
        b=vb
        winner='endpoint' if l==0 else ('A' if ra<rb else 'B')
        rows.append({
          'lambda':mp.nstr(l,20),'v_A':mp.nstr(va,20),'R_A':mp.nstr(ra,20),
          'v_B':mp.nstr(vb,20),'R_B':mp.nstr(rb,20),
          'C_coal':mp.nstr(coal(l),20),'C_strong':mp.nstr(strong(l),20),
          'R_at_b':mp.nstr(cost(mp.mpf(10),l),20),
          'sampled_best_identity':winner,
          'global_status':'certified_crossing_at_lambda_AB' if abs(l-LAM_AB)<mp.mpf('1e-20') else 'numerical_regular_comparison',
          'Ccoal_minus_A':mp.nstr(coal(l)-ra,20)
        })
    with (ROOT/'tangent_cost_vs_lambda.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    with (ROOT/'smallz_global_regimes.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f)
        w.writerow(['lambda_range','observed_best','proof_status','notes'])
        w.writerow(['0','coalescent endpoint','exact at tangent level','A emerges from this endpoint'])
        w.writerow(['(0, lambda_AB)','near-zero regular branch connects numerically to A','numerical over plotted interval; independently certified locally near lambda_AB and lambda=0','no positive coalescent global regime'])
        w.writerow(['lambda_AB','A=B','interval certified by Stage 11/12','transverse separated minima'])
        w.writerow(['(lambda_AB, 0.2]','regular B','numerical over plotted interval; certified locally near lambda_AB','no statement beyond 0.2'])
    x=np.array([float(r['lambda']) for r in rows]);a_cost=np.array([float(r['R_A']) for r in rows]);b_cost=np.array([float(r['R_B']) for r in rows]);c=np.array([float(r['C_coal']) for r in rows]);cs=np.array([float(r['C_strong']) for r in rows])
    fig,(ax,strip)=plt.subplots(2,1,figsize=(8.5,5.8),gridspec_kw={'height_ratios':[5,0.7]},sharex=True)
    ax.plot(x,a_cost,label='regular A (numerical continuation)',color='#1464a0',lw=2)
    ax.plot(x,b_cost,label='regular B (numerical continuation)',color='#d06322',lw=2)
    ax.plot(x,c,label='coalescent endpoint',color='#555555',lw=1.5,ls='--')
    ax.plot(x,cs,label='fixed strong-pair control',color='#888888',lw=1.1,ls=':')
    ax.scatter([float(LAM_AB)],[float(cost(mp.mpf('-4.15344048664877'),LAM_AB))],color='black',zorder=5,s=25)
    ax.axvline(float(LAM_AB),color='black',alpha=.4,lw=.8)
    ax.annotate(r'certified A/B crossing, $\lambda_{21}$',xy=(float(LAM_AB),float(cost(mp.mpf('-4.15344048664877'),LAM_AB))),xytext=(.085,.16),arrowprops=dict(arrowstyle='->',lw=.8),fontsize=9)
    ax.set_ylabel('tangent cost');ax.set_ylim(-.005,.34);ax.set_xlim(0,.2);ax.grid(alpha=.18);ax.legend(fontsize=8,loc='upper left')
    strip.axvspan(0,float(LAM_AB),color='#1464a0',alpha=.35)
    strip.axvspan(float(LAM_AB),.2,color='#d06322',alpha=.35)
    strip.scatter([0],[.5],color='black',s=30,zorder=5)
    strip.text(.03,.5,'A: sampled best',ha='center',va='center',fontsize=8)
    strip.text(.14,.5,'B: sampled best',ha='center',va='center',fontsize=8)
    strip.set_yticks([]);strip.set_ylim(0,1);strip.set_xlabel(r'$\lambda=\epsilon/z^2$');strip.spines[['left','right','top']].set_visible(False)
    fig.text(.12,.005,'Full-interval global A/B ordering is numerical; endpoint exclusion and the A/B crossing have proofs.',fontsize=7.7)
    fig.tight_layout(rect=(0,.025,1,1))
    fig.savefig(ROOT/'smallz_global_cost_figure.pdf')
    fig.savefig(ROOT/'paper'/'figures'/'fig_stage13_global_cost.pdf')
    plt.close(fig)
    print(json.dumps({'rows':len(rows),'certificate':cert,'R_11_at_lambda_gamma':mp.nstr(cost(mp.mpf(11),LAM_G),25),'Ccoal_at_lambda_gamma':mp.nstr(coal(LAM_G),25)},indent=2))

if __name__=='__main__':main()

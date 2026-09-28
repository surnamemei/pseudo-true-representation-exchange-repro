"""Directed interval endpoint constants for the continuum tangent functional."""
import json
from pathlib import Path
from mpmath import iv,mp

ROOT=Path(__file__).resolve().parent
iv.dps=75;mp.dps=85
m2=iv.mpf(1)/12;m4=iv.mpf(1)/80;m6=iv.mpf(1)/448
V2=m4-m2*m2;V3=m6-m4*m4/m2

def D(x,k=0):
    y=2*iv.sin(x/2)/x
    for j in range(1,k+1):
        y=(2*iv.sin(x/2+j*iv.pi/2)/2**j-j*y)/x
    return y

def bounds(x):return [str(mp.mpf(x.a)),str(mp.mpf(x.b))]

def main():
    d=D(iv.mpf(10));dp=D(iv.mpf(10),1);d2=D(iv.mpf(10),2)
    Q=d2+m2*d
    M=D(iv.mpf(10),3)+m4/m2*dp
    lamg=V2/Q
    Kcoal=1-d*d-dp*dp/m2-Q*Q/V2
    A=1-d*d-dp*dp/m2
    Rb=V2-Q*Q/A
    threshold=iv.sqrt(Rb/Kcoal)
    out={
      'M':bounds(M),'Q':bounds(Q),'V2':bounds(V2),'V3':bounds(V3),
      'lambda_gamma':bounds(lamg),
      'v0_linear_coefficient':bounds(3*M/V3),
      'endpoint_cost_gain_coefficient':bounds(M*M/V3),
      'coalescent_quadratic_coefficient':bounds(Kcoal),
      'R_at_b_constant':bounds(Rb),
      'fixed_b_witness_threshold':bounds(threshold),
      'no_positive_coalescent_phase':bool(M.b<0 and Q.a>0 and V3.a>0 and threshold.b<lamg.a)
    }
    (ROOT/'stage14_continuum_endpoint_interval.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()

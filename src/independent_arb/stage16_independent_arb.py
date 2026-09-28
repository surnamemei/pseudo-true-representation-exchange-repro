"""Independent Arb replay of the frozen continuum/all-large-N certificate.

This file does not import the Stage 14/15 evaluators.  It uses python-flint
balls, a separate second-order dual implementation, exact rational kernel
coefficients, and the frozen partition only as a list of obligations.
"""
from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

from flint import arb, ctx

ctx.dps = 85
ROOT = Path(__file__).resolve().parent
N0 = 35377
TM = arb(0).union(arb(1) / (N0 * N0))
ZERO = arb(0)
ONE = arb(1)
ERR = arb(0, "1e-30")
NEAR_ERR = arb(0, "1e-25")
REF = json.loads((ROOT / "stage14_continuum_candidate.json").read_text())


def B(x):
    if isinstance(x, arb):
        return x
    if isinstance(x, Fraction):
        return arb(x.numerator) / x.denominator
    return arb(x)


def hull(a, b):
    return B(a).union(B(b))


def lower(x):
    return x.lower()


def upper(x):
    return x.upper()


def positive(x):
    return lower(x) > 0


def negative(x):
    return upper(x) < 0


def away(x):
    return positive(x) or negative(x)


def minimum_ball(x, y):
    return x if float(lower(x)) < float(lower(y)) else y


# H(w)=sqrt(w)/sin(sqrt(w)); recurrence derived from H(w)*sinc(sqrt(w))=1.
H = [Fraction(1)]
for m in range(1, 9):
    H.append(-sum(H[m-j] * Fraction((-1)**j, math.factorial(2*j+1))
                  for j in range(1, m+1)))


def sinc_deriv(q, k):
    """Derivatives of 2 sin(q/2)/q, with a separate analytic zero chart."""
    # The zero chart is valid much farther out; using |q|<=2 also avoids
    # an unnecessary wide rational enclosure when a ball grazes |q|=1.
    if (k <= 4 and float(abs(q).upper()) <= 2) or float(abs(q).upper()) <= 1:
        y = ZERO
        for m in range((k+1)//2, 38):
            p = 2*m-k
            y += B(Fraction((-1)**m, (2*m+1) * math.factorial(p) * 2**(2*m))) * q**p
        return y + arb(0, "1e-50")
    s, c = (q/2).sin(), (q/2).cos()
    y = 2*s/q
    for j in range(1, k+1):
        trig = (s, c, -s, -c)[j % 4]
        y = (2*trig / (2**j) - j*y) / q
    return y


def kernel(q, k=0):
    """Uniform t-ball kernel derivatives; Cauchy remainder checked below."""
    ds = [sinc_deriv(q, j) for j in range(k+1)]
    result = ZERO
    for j in range(k+1):
        hp = ZERO
        for m in range((j+1)//2, 9):
            p = 2*m
            hp += B(H[m]) * TM**m * q**(p-j) * (math.factorial(p)//math.factorial(p-j)) / 4**m
        result += math.comb(k, j) * ds[k-j] * hp
    return result + ERR


M2 = (ONE-TM)/12
M4 = ONE/80-TM/24+7*TM**2/240


@dataclass
class Dual:
    x: arb
    g: tuple
    h: tuple

    @classmethod
    def const(cls, x):
        return cls(B(x), (ZERO, ZERO), ((ZERO, ZERO), (ZERO, ZERO)))

    @classmethod
    def variable(cls, x, i):
        return cls(B(x), tuple(ONE if j == i else ZERO for j in range(2)),
                   ((ZERO, ZERO), (ZERO, ZERO)))

    def __add__(self, other):
        other = other if isinstance(other, Dual) else Dual.const(other)
        return Dual(self.x+other.x,
                    tuple(self.g[i]+other.g[i] for i in range(2)),
                    tuple(tuple(self.h[i][j]+other.h[i][j] for j in range(2)) for i in range(2)))
    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.x, tuple(-a for a in self.g),
                    tuple(tuple(-a for a in row) for row in self.h))

    def __sub__(self, other):
        return self + (-other if isinstance(other, Dual) else -B(other))

    def __rsub__(self, other):
        return other + (-self)

    def __mul__(self, other):
        other = other if isinstance(other, Dual) else Dual.const(other)
        return Dual(self.x*other.x,
                    tuple(self.g[i]*other.x+self.x*other.g[i] for i in range(2)),
                    tuple(tuple(self.h[i][j]*other.x+self.x*other.h[i][j]
                                +self.g[i]*other.g[j]+self.g[j]*other.g[i]
                                for j in range(2)) for i in range(2)))
    __rmul__ = __mul__

    def reciprocal(self):
        x1 = ONE/self.x
        return Dual(x1, tuple(-self.g[i]*x1*x1 for i in range(2)),
                    tuple(tuple(2*self.g[i]*self.g[j]*x1**3-self.h[i][j]*x1*x1
                                for j in range(2)) for i in range(2)))

    def __truediv__(self, other):
        other = other if isinstance(other, Dual) else Dual.const(other)
        return self*other.reciprocal()

    def __rtruediv__(self, other):
        return Dual.const(other)*self.reciprocal()


def kd(q, order=0):
    """Composition of the scalar kernel derivative with a dual variable."""
    d0, d1, d2 = (kernel(q.x, order+j) for j in range(3))
    return Dual(d0, tuple(d1*q.g[i] for i in range(2)),
                tuple(tuple(d2*q.g[i]*q.g[j]+d1*q.h[i][j] for j in range(2)) for i in range(2)))


def tangent(v, lam):
    v, lam = Dual.variable(v, 0), Dual.variable(lam, 1)
    d, dp = kd(v), kd(v, 1)
    q0 = -M2-lam*kernel(arb(10))
    q1 = 2*lam*kernel(arb(10), 1)
    aa = 1-d*d-dp*dp/M2
    bb = kd(v, 2)-lam*kd(10-v)-d*q0+dp*q1/(2*M2)
    cc = M4-2*lam*kernel(arb(10), 2)+lam*lam-q0*q0-q1*q1/(4*M2)
    return cc-bb*bb/aa, aa, bb/aa


def system(x):
    a, aa, beta_a = tangent(x[0], x[2])
    b, ab, beta_b = tangent(x[1], x[2])
    f = [a.g[0], b.g[0], a.x-b.x]
    j = [[a.h[0][0], ZERO, a.h[0][1]],
         [ZERO, b.h[0][0], b.h[0][1]],
         [a.g[0], -b.g[0], a.g[1]-b.g[1]]]
    return f, j, a, b, aa, ab, beta_a, beta_b


def inverse3(a):
    t = [[float(a[i][j].mid()) for j in range(3)]+[float(i==j) for j in range(3)] for i in range(3)]
    for p in range(3):
        k = max(range(p,3), key=lambda i: abs(t[i][p]))
        t[p], t[k] = t[k], t[p]
        scale=t[p][p]
        t[p]=[v/scale for v in t[p]]
        for i in range(3):
            if i != p:
                z=t[i][p]
                t[i]=[t[i][j]-z*t[p][j] for j in range(6)]
    return [[B(repr(t[i][j+3])) for j in range(3)] for i in range(3)]


def local_replay():
    x0 = [B(REF[k]) for k in ("v_A", "v_B", "lambda")]
    tmax = 1/N0**2
    radii = [max(1e-8,10000*tmax),max(1e-8,10000*tmax),max(1e-9,2*tmax)]
    x = [x0[i] + arb(0,repr(radii[i])) for i in range(3)]
    f0, j0, *_ = system(x0)
    _, jx, a, b, aa, ab, beta_a, beta_b = system(x)
    c = inverse3(j0)
    e = [[B(i==j)-sum(c[i][k]*jx[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    kbox = [x0[i]-sum(c[i][k]*f0[k] for k in range(3))
            +sum(e[i][j]*(x[j]-x0[j]) for j in range(3)) for i in range(3)]
    inside = all(lower(kbox[i])>upper(lower(x[i])) and
                 upper(kbox[i])<lower(upper(x[i])) for i in range(3))
    contraction = max(sum(float(abs(e[i][j]).upper()) for j in range(3)) for i in range(3))
    signs = (positive(a.h[0][0]) and positive(b.h[0][0]) and
             positive(a.g[1]-b.g[1]) and positive(aa.x) and positive(ab.x)
             and positive(beta_a.x) and negative(beta_b.x))
    return x, {"krawczyk_inside":inside,"contraction_upper":contraction,
               "signs":signs,"curvature_A_lower":str(lower(a.h[0][0])),
               "curvature_B_lower":str(lower(b.h[0][0])),
               "slope_lower":str(lower(a.g[1]-b.g[1])),
               "beta_A_lower":str(lower(beta_a.x)),
               "abs_beta_B_lower":str(-upper(beta_b.x)),
               "K":[str(v) for v in kbox],"X":[str(v) for v in x]}


def make_eval(lam):
    q0 = -M2-lam*kernel(arb(10))
    q1 = 2*lam*kernel(arb(10),1)
    base = M4-2*lam*kernel(arb(10),2)+lam*lam-q0*q0-q1*q1/(4*M2)
    # Zero chart: direct Maclaurin coefficients, independent Arb operations.
    K = 46
    dc = [ZERO for _ in range(K+3)]
    for p in range(0,K+3,2):
        m=p//2
        dc[p]=sum(B(Fraction((-1)**(m-j),math.factorial(2*(m-j)+1)*2**(2*(m-j))))
                  *B(H[j])*TM**j/4**j for j in range(min(m,8)+1))
    dp = [(i+1)*dc[i+1] for i in range(K+1)]
    db = [(-1)**i*kernel(arb(10),i)/math.factorial(i) for i in range(K+1)]
    ac = [(ONE if i==0 else ZERO)-sum(dc[j]*dc[i-j]+dp[j]*dp[i-j]/M2 for j in range(i+1))
          for i in range(K+1)][4:45]
    bc = [(i+2)*(i+1)*dc[i+2]-lam*db[i]-dc[i]*q0+dp[i]*q1/(2*M2)
          for i in range(K+1)][2:43]

    def poly(coeff, q, order):
        y=ZERO
        for i in range(len(coeff)-1,order-1,-1):
            y=y*q+coeff[i]*(math.factorial(i)//math.factorial(i-order))
        return y+NEAR_ERR

    def evaluate(lo,hi):
        v=hull(lo,hi)
        if max(abs(float(lo)),abs(float(hi)))<=1:
            aa,ap,app=(poly(ac,v,i) for i in range(3))
            bb,bp,bpp=(poly(bc,v,i) for i in range(3))
        else:
            d=[kernel(v,i) for i in range(5)]
            aa=1-d[0]*d[0]-d[1]*d[1]/M2
            ap=-2*d[0]*d[1]-2*d[1]*d[2]/M2
            app=-2*(d[1]*d[1]+d[0]*d[2])-2*(d[2]*d[2]+d[1]*d[3])/M2
            bb=d[2]-lam*kernel(10-v)-d[0]*q0+d[1]*q1/(2*M2)
            bp=d[3]+lam*kernel(10-v,1)-d[1]*q0+d[2]*q1/(2*M2)
            bpp=d[4]-lam*kernel(10-v,2)-d[2]*q0+d[3]*q1/(2*M2)
        if not positive(aa):
            return None
        cost=base-bb*bb/aa
        grad=-2*bb*bp/aa+bb*bb*ap/aa**2
        curv=-2*(bp*bp+bb*bpp)/aa+4*bb*bp*ap/aa**2+bb*bb*app/aa**2-2*bb*bb*ap*ap/aa**3
        return cost,grad,curv,aa
    return evaluate,base,q0,q1


def independent_remainder():
    # q-Cauchy radius 55 about any real |q|<=110; t-Cauchy ratio z=(165/N)^2.
    half=B(Fraction(1,2))
    eta=half.sinh()/half-1
    def rem(n):
        z=(B(165)/n)**2
        return (B(Fraction(55,2))).exp()/(1-eta)*z**9/(1-z)
    n,prev=rem(N0),rem(N0-2)
    return {"bound_N0":str(n),"bound_previous_odd":str(prev),
            "under_1e-30":upper(n)<B("1e-30"),
            "previous_over_1e-30":lower(prev)>B("1e-30")}


def run():
    x,local=local_replay()
    print('Arb local:',local['krawczyk_inside'],local['signs'],flush=True)
    evaluate,base,q0,q1=make_eval(x[2])
    trial=evaluate(REF['v_A'],REF['v_A'])
    incumbent=upper(trial[0])
    rows=list(csv.DictReader((ROOT/'stage14_continuum_global_partition.csv').open()))
    ordered=sorted(rows,key=lambda r:Fraction(r['lo']))
    covered=(Fraction(ordered[0]['lo']) == -100 and Fraction(ordered[-1]['hi']) == 100
             and all(Fraction(a['hi'])==Fraction(b['lo']) for a,b in zip(ordered,ordered[1:])))
    counts={};failed=[];refined=0;initial_unpassed=[];min_cost=None;min_grad=None;min_curv=None
    output=[]
    def certify_refined(lo,hi,root=False,depth=0):
        nonlocal refined
        refined += 1
        z=evaluate(lo,hi)
        if z is not None:
            cost,grad,curv,_=z
            if root:
                if positive(curv): return True
            else:
                if positive(lower(cost)-incumbent) or away(grad): return True
                ga,gb=evaluate(lo,lo)[1],evaluate(hi,hi)[1]
                if away(curv) and ((positive(ga) and positive(gb)) or (negative(ga) and negative(gb))):
                    return True
        if depth >= 16:
            return False
        mid=(lo+hi)/2
        return certify_refined(lo,mid,root,depth+1) and certify_refined(mid,hi,root,depth+1)

    for i,r in enumerate(rows):
        lo,hi=r['lo'],r['hi'];kind=r['predicate']
        z=evaluate(lo,hi)
        ok=False; margin=None
        if z is not None:
            cost,grad,curv,_=z
            if kind=='cost_excluded':
                margin=lower(cost)-incumbent
                ok=positive(margin)
                if ok:min_cost=margin if min_cost is None else minimum_ball(margin,min_cost)
            elif kind=='gradient_excluded':
                ok=away(grad)
                margin=grad.abs_lower()
                if ok:min_grad=margin if min_grad is None else minimum_ball(margin,min_grad)
            elif kind=='monotone_derivative_no_root':
                ga,gb=evaluate(lo,lo)[1],evaluate(hi,hi)[1]
                ok=away(curv) and ((positive(ga) and positive(gb)) or (negative(ga) and negative(gb)))
                margin=curv.abs_lower()
                if ok:min_curv=margin if min_curv is None else minimum_ball(margin,min_curv)
            elif kind.startswith('unique_root_'):
                ga,gb=evaluate(lo,lo)[1],evaluate(hi,hi)[1]
                ok=positive(curv) and negative(ga) and positive(gb)
                margin=lower(curv)
        counts[kind]=counts.get(kind,0)+1
        if not ok:
            initial_unpassed.append(i)
            aa,bb=Fraction(lo),Fraction(hi)
            if kind.startswith('unique_root_'):
                ga,gb=evaluate(lo,lo)[1],evaluate(hi,hi)[1]
                ok=negative(ga) and positive(gb) and certify_refined(aa,bb,root=True)
            else:
                ok=certify_refined(aa,bb)
            if not ok: failed.append((i,lo,hi,kind))
        output.append({'index':i,'lo':lo,'hi':hi,'predicate':kind,'pass':ok,
                       'refined':i in initial_unpassed,
                       'initial_margin':str(margin) if margin is not None else ''})
        if (i+1)%5000==0:print('Arb cells:',i+1,'failures:',len(failed),flush=True)
    # Endpoint and tail use independently evaluated Arb moment and sinc envelopes.
    q2=-M4+x[2]*kernel(arb(10),2)-M2*q0
    coal=base-q2*q2/(M4-M2*M2)
    V=B(100)
    d0=2*arb.pi()/V;d1=arb.pi()/V;d2=arb.pi()/(2*V)
    a_tail=1-d0*d0-d1*d1/M2
    b_tail=d2+upper(x[2])*2*arb.pi()/(V-10)+d0*upper(abs(q0))+d1*upper(abs(q1))/(2*lower(M2))
    tail=base-b_tail*b_tail/a_tail
    roots=[r for r in rows if r['predicate'].startswith('unique_root_')]
    linked=(len(roots)==2 and all(
        B(Fraction(r['lo'])) < lower(lower(x[0 if '_A_' in r['predicate'] else 1]))
        and upper(upper(x[0 if '_A_' in r['predicate'] else 1])) < B(Fraction(r['hi']))
        for r in roots))
    bounds=independent_remainder()
    summary={'backend':'python-flint/Arb 0.9.0','precision_decimal_digits':ctx.dps,
             'N0':N0,'t_max':str(upper(TM)),'partition_cells':len(rows),'partition_cover':covered,
             'predicate_counts':counts,'passed':len(rows)-len(failed),
             'original_cells_refined':len(initial_unpassed),'refinement_nodes':refined,
             'refined_original_indices':initial_unpassed,'failed':len(failed),
             'first_failures':failed[:20],'root_boxes_linked':linked,'local':local,
             'cost_margin_lower':str(lower(min_cost)) if min_cost is not None else None,
             'gradient_margin_lower':str(lower(min_grad)) if min_grad is not None else None,
             'curvature_sign_margin_lower':str(lower(min_curv)) if min_curv is not None else None,
             'coalescent_gap_lower':str(lower(coal-incumbent)),
             'tail_gap_lower':str(lower(tail-incumbent)),
             'kernel_remainder':bounds}
    summary['verified']=(covered and linked and len(failed)==0 and local['krawczyk_inside']
                         and local['signs'] and positive(coal-incumbent) and positive(tail-incumbent)
                         and bounds['under_1e-30'] and bounds['previous_over_1e-30'])
    with (ROOT/'independent_N0_replay.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=output[0]);w.writeheader();w.writerows(output)
    (ROOT/'independent_N0_replay.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps({k:v for k,v in summary.items() if k not in ('local','first_failures')},indent=2),flush=True)
    return summary


if __name__=='__main__':
    run()

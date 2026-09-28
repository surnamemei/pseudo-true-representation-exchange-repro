"""Independent, limited MPFR directed-rounding replay of frozen certificates.

This script uses gmpy2/MPFR, not mpmath.iv. It deliberately checks selected
root and outer-cell predicates, not the complete global partition.
"""
from __future__ import annotations

import csv
import json
import math
import os
import re
import sys
from decimal import Decimal, getcontext
from pathlib import Path

if os.name == "nt":
    candidate = Path(sys.executable).parent / "Library" / "bin"
    if candidate.is_dir():
        _mpfr_dll_dir = os.add_dll_directory(str(candidate))
import gmpy2 as G

ROOT = Path(__file__).resolve().parent
getcontext().prec = 125
BASE = G.get_context().copy()
BASE.precision = 360
G.get_context().precision = 360
DOWN = BASE.copy()
DOWN.round = G.RoundDown
UP = BASE.copy()
UP.round = G.RoundUp
NEAR = BASE.copy()
NEAR.round = G.RoundToNearest


def rd(fn, *args):
    with G.local_context(DOWN):
        return fn(*args)


def ru(fn, *args):
    with G.local_context(UP):
        return fn(*args)


def rn(fn, *args):
    with G.local_context(NEAR):
        return fn(*args)


class I:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        if isinstance(lo, I) and hi is None:
            self.lo, self.hi = lo.lo, lo.hi
            return
        if hi is None:
            hi = lo
        self.lo = rd(G.mpfr, str(lo))
        self.hi = ru(G.mpfr, str(hi))
        assert self.lo <= self.hi

    @classmethod
    def raw(cls, lo, hi):
        x = object.__new__(cls)
        x.lo, x.hi = lo, hi
        assert lo <= hi
        return x

    def __str__(self):
        return "[" + str(self.lo) + ", " + str(self.hi) + "]"

    def mid(self):
        return rn(lambda a, b: (a + b) / 2, self.lo, self.hi)

    def absmax(self):
        return max(abs(self.lo), abs(self.hi))

    def __add__(a, b):
        b = I(b)
        return I.raw(rd(lambda x, y: x + y, a.lo, b.lo),
                     ru(lambda x, y: x + y, a.hi, b.hi))

    __radd__ = __add__

    def __neg__(a):
        return I.raw(-a.hi, -a.lo)

    def __sub__(a, b):
        return a + -I(b)

    def __rsub__(a, b):
        return I(b) - a

    def __mul__(a, b):
        b = I(b)
        pairs = [(x, y) for x in (a.lo, a.hi) for y in (b.lo, b.hi)]
        return I.raw(min(rd(lambda x, y: x * y, x, y) for x, y in pairs),
                     max(ru(lambda x, y: x * y, x, y) for x, y in pairs))

    __rmul__ = __mul__

    def reciprocal(a):
        assert a.lo > 0 or a.hi < 0, "division interval crosses zero"
        return I.raw(rd(lambda x: 1 / x, a.hi),
                     ru(lambda x: 1 / x, a.lo))

    def __truediv__(a, b):
        return a * I(b).reciprocal()

    def __rtruediv__(a, b):
        return I(b) * a.reciprocal()

    def __pow__(a, n):
        assert isinstance(n, int) and n >= 0
        if n == 0:
            return I(1)
        if n == 1:
            return a
        if n == 2:
            upper = max(ru(lambda x: x * x, a.lo),
                        ru(lambda x: x * x, a.hi))
            if a.lo <= 0 <= a.hi:
                return I.raw(G.mpfr(0), upper)
            lower = min(rd(lambda x: x * x, a.lo),
                        rd(lambda x: x * x, a.hi))
            return I.raw(lower, upper)
        out = I(1)
        for _ in range(n):
            out = out * a
        return out

    def sqrt(a):
        assert a.lo >= 0
        return I.raw(rd(G.sqrt, a.lo), ru(G.sqrt, a.hi))

    @staticmethod
    def pi():
        return I.raw(rd(G.const_pi), ru(G.const_pi))

    def _trig(a, cosine=False):
        fun = G.cos if cosine else G.sin
        lows = [rd(fun, a.lo), rd(fun, a.hi)]
        highs = [ru(fun, a.lo), ru(fun, a.hi)]
        step = I.pi() / 2
        start = math.floor(float(a.lo) / (math.pi / 2)) - 4
        stop = math.ceil(float(a.hi) / (math.pi / 2)) + 4
        if stop - start > 100:
            return I(-1, 1)
        for k in range(start, stop + 1):
            peak = step * k
            if peak.lo <= a.hi and peak.hi >= a.lo:
                if cosine and k % 2 == 0:
                    highs.append(G.mpfr(1) if k % 4 == 0 else G.mpfr(-1))
                    lows.append(G.mpfr(1) if k % 4 == 0 else G.mpfr(-1))
                if not cosine and k % 2:
                    highs.append(G.mpfr(1) if k % 4 == 1 else G.mpfr(-1))
                    lows.append(G.mpfr(1) if k % 4 == 1 else G.mpfr(-1))
        return I.raw(min(lows), max(highs))

    def sin(a):
        return a._trig(False)

    def cos(a):
        return a._trig(True)


Z = I(0)
O = I(1)


class Jet:
    __slots__ = ("v", "g", "H")

    def __init__(self, v, g=None, H=None):
        self.v = I(v)
        self.g = g if g is not None else (Z, Z)
        self.H = H if H is not None else ((Z, Z), (Z, Z))

    @staticmethod
    def var(v, k):
        return Jet(v, (O if k == 0 else Z, O if k == 1 else Z))

    def __add__(a, b):
        b = b if isinstance(b, Jet) else Jet(b)
        return Jet(a.v + b.v,
                   tuple(a.g[i] + b.g[i] for i in range(2)),
                   tuple(tuple(a.H[i][j] + b.H[i][j] for j in range(2))
                         for i in range(2)))

    __radd__ = __add__

    def __neg__(a):
        return Jet(-a.v, tuple(-x for x in a.g),
                   tuple(tuple(-x for x in row) for row in a.H))

    def __sub__(a, b):
        return a + -(b if isinstance(b, Jet) else Jet(b))

    def __rsub__(a, b):
        return Jet(b) - a

    def __mul__(a, b):
        b = b if isinstance(b, Jet) else Jet(b)
        return Jet(a.v * b.v,
                   tuple(a.g[i] * b.v + a.v * b.g[i] for i in range(2)),
                   tuple(tuple(a.H[i][j] * b.v + a.g[i] * b.g[j]
                               + a.g[j] * b.g[i] + a.v * b.H[i][j]
                               for j in range(2)) for i in range(2)))

    __rmul__ = __mul__

    def reciprocal(a):
        v = a.v.reciprocal()
        return Jet(v,
                   tuple(-a.g[i] * v**2 for i in range(2)),
                   tuple(tuple(2*a.g[i]*a.g[j]*v**3-a.H[i][j]*v**2
                               for j in range(2)) for i in range(2)))

    def __truediv__(a, b):
        return a * (b if isinstance(b, Jet) else Jet(b)).reciprocal()

    def __rtruediv__(a, b):
        return Jet(b) * a.reciprocal()

    def cos(a):
        c, s = a.v.cos(), a.v.sin()
        return Jet(c, tuple(-s*a.g[i] for i in range(2)),
                   tuple(tuple(-s*a.H[i][j]-c*a.g[i]*a.g[j]
                               for j in range(2)) for i in range(2)))


def D(z, N):
    out = Jet(1)
    for k in range(1, (N-1)//2+1):
        out += 2*(z*k/N).cos()
    return out


def d0(z, N):
    out = I(1)
    for k in range(1, (N-1)//2+1):
        out += 2*(I(z)*k/N).cos()
    return out


def energy(e, N):
    out = I(0)
    for t in range(-(N-1)//2, (N-1)//2+1):
        c = (I(2)*t/N).cos()
        out += 4*c**2 + e**2 - 4*e*c*(I(10)*t/N).cos()
    return out


def objective(u1, u2, e, N):
    a, b = Jet.var(u1, 0), Jet.var(u2, 1)
    h1 = D(a+2,N)+D(a-2,N)-D(a-10,N)*e
    h2 = D(b+2,N)+D(b-2,N)-D(b-10,N)*e
    g = D(b-a,N)
    den = N*N-g*g
    cap = N*(h1*h1+h2*h2)-2*g*h1*h2
    return Jet(energy(e,N))-cap/den


def slope_branch(u1, u2, e, N):
    g = d0(u2-u1, N)
    den = N*N-g*g
    h1 = d0(u1+2,N)+d0(u1-2,N)-e*d0(u1-10,N)
    h2 = d0(u2+2,N)+d0(u2-2,N)-e*d0(u2-10,N)
    z1, z2 = d0(u1-10,N), d0(u2-10,N)
    Eder = 2*N*e
    for t in range(-(N-1)//2, (N-1)//2+1):
        Eder -= 4*(I(2)*t/N).cos()*(I(10)*t/N).cos()
    capder = -2*N*(h1*z1+h2*z2)+2*g*(z1*h2+h1*z2)
    return Eder-capder/den


def parse_original(s):
    nums = re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", str(s))
    if not nums:
        return None
    return I(nums[0], nums[-1])


def root_check(N, center, eps, radius="1e-10"):
    x0 = [I(v) for v in center]
    r = I(radius)
    X = [I.raw(rd(lambda a,b:a-b,x.lo,r.hi),
               ru(lambda a,b:a+b,x.hi,r.hi)) for x in x0]
    j0 = objective(*x0,eps,N)
    jx = objective(*X,eps,N)
    H = jx.H
    h0 = objective(*x0,eps,N).H
    det0 = rn(lambda a,b,c,d:a*d-b*c,
              h0[0][0].mid(),h0[0][1].mid(),
              h0[1][0].mid(),h0[1][1].mid())
    c00 = rn(lambda a,b:b/a, det0, h0[1][1].mid())
    c01 = rn(lambda a,b:-b/a, det0, h0[0][1].mid())
    c10 = rn(lambda a,b:-b/a, det0, h0[1][0].mid())
    c11 = rn(lambda a,b:b/a, det0, h0[0][0].mid())
    C = [[I(c00),I(c01)],[I(c10),I(c11)]]
    M = [[I(int(i==j))-sum(C[i][k]*H[k][j] for k in range(2))
          for j in range(2)] for i in range(2)]
    Y = [X[i]-x0[i] for i in range(2)]
    K = [x0[i]-sum(C[i][k]*j0.g[k] for k in range(2))
         +sum(M[i][k]*Y[k] for k in range(2)) for i in range(2)]
    included = all(K[i].lo>X[i].lo and K[i].hi<X[i].hi for i in range(2))
    det = H[0][0]*H[1][1]-H[0][1]*H[1][0]
    positive = H[0][0].lo>0 and det.lo>0
    J = j0.v
    for i in range(2):
        J += j0.g[i]*Y[i]
    for i in range(2):
        for k in range(2):
            J += I("0.5")*H[i][k]*Y[i]*Y[k]
    return dict(X=X,K=K,J=J,H11=H[0][0],det=det,
                included=included,positive=positive,
                contraction=max(sum(M[i][j].absmax() for j in range(2))
                                for i in range(2)))


def center_projected(m,h,e,E,N):
    tmax=(N-1)//2
    h1=I(2)-e;h2=I(0)
    c1sq=I(1);c2sq=I(0);w1sq=I(0);w2sq=I(0)
    for t in range(1,tmax+1):
        z=h*t
        co,si=z.cos(),z.sin()
        sc=si/z
        sp=(z*co-si)/(z*z)
        p=(I(2)*t/N-m*t).cos()+(-I(2)*t/N-m*t).cos()
        p-=e*(I(10)*t/N-m*t).cos()
        q=(I(2)*t/N-m*t).sin()+(-I(2)*t/N-m*t).sin()
        q-=e*(I(10)*t/N-m*t).sin()
        h1+=2*co*p;h2+=2*t*sc*q
        c1sq+=2*co**2;c2sq+=2*t*t*sc**2
        w1sq+=2*t*t*si**2;w2sq+=2*t**4*sp**2
    J=E-h1**2/c1sq-h2**2/c2sq
    return J,c1sq.sqrt(),c2sq.sqrt(),w1sq.sqrt(),w2sq.sqrt()


def outer_bound(row,e,N):
    a,b,c,d=[Decimal(row[k]) for k in ("m_lo","m_hi","h_lo","h_hi")]
    inflate=Decimal("1e-12")
    m=I(str((a+b)/2));h=I(str((c+d)/2))
    rm=I(str((b-a)/2+inflate));rh=I(str((d-c)/2+inflate))
    E=energy(e,N)
    Jc,c1,c2,w1,w2=center_projected(m,h,e,E,N)
    tmax=(N-1)//2
    t2=I(sum(t**4 for t in range(-tmax,tmax+1))).sqrt()
    t3=I(sum(t**6 for t in range(-tmax,tmax+1))).sqrt()
    D1=w1*rh+I("0.5")*t2*rh**2
    D2=w2*rh+t3*rh**2/6
    assert c1.lo>D1.hi and c2.lo>D2.hi
    theta=I(2).sqrt()*tmax*rm+((D1/(c1-D1))**2+
                                   (D2/(c2-D2))**2).sqrt()
    assert Jc.lo>0
    residual=Jc.sqrt()-E.sqrt()*theta
    return I(0) if residual.lo<=0 else I.raw(
        rd(lambda x:x*x,residual.lo),ru(lambda x:x*x,residual.hi))


def concise(x, digits=24):
    return "["+format(x.lo,f".{digits}g")+","+format(x.hi,f".{digits}g")+"]"


def main():
    c21=json.loads((ROOT/"validated_globality_certificate.json").read_text())
    c31=json.loads((ROOT/"N31_global_certificate.json").read_text())
    local21=json.loads((ROOT/"stage2_local_interval_checks.json").read_text())
    centers={
        21:{k:c21["cases"][1][k+"_box"] for k in ("A","B")},
        31:{"A":[c31["high_precision_numerical_root"]["u_A1"],
                 c31["high_precision_numerical_root"]["u_A2"]],
            "B":[c31["high_precision_numerical_root"]["u_B1"],
                 c31["high_precision_numerical_root"]["u_B2"]]},
    }
    eps={
        21:[c21["crossing"]["epsilon_lower"],c21["crossing"]["epsilon_upper"]],
        31:c31["epsilon_bracket"],
    }
    rows=[]
    root_data={}
    for N,cert in ((21,c21),(31,c31)):
        er=I(*eps[N])
        root_data[N]={}
        for branch in ("A","B"):
            q=root_check(N,centers[N][branch],er)
            root_data[N][branch]=q
            orig=(cert["crossing"] if N==21 else cert["crossing_checks"][2])
            if N==21:
                original_root=next(v for v in local21
                                   if v["label"]=="crossing" and v["branch"]==branch)
                original_K=original_root["K1"]+"; "+original_root["K2"]
                original_H11=original_root["H11_lower"]
                original_det=original_root["determinant_lower"]
            else:
                original_K=str(orig[branch]["K"])
                original_H11=orig[branch]["H11_lower"]
                original_det=orig[branch]["Hessian_det_lower"]
            for typ,orig_value,ind_interval,result in (
                ("root_inclusion",original_K,
                 "; ".join(concise(v) for v in q["K"]),q["included"]),
                ("Hessian_H11",original_H11,concise(q["H11"]),q["H11"].lo>0),
                ("Hessian_det",original_det,concise(q["det"]),q["positive"]),
            ):
                rows.append(dict(N=N,check=typ,branch=branch,case="full_bracket",
                    source=("stage2_local_interval_checks.json (point) and certificate JSON (bracket)"
                            if N==21 else "N31_global_certificate.json (bracket)"),
                    original_interval=orig_value,
                    independent_interval=ind_interval,
                    agreement="PASS" if result else "FAIL",notes="radius 1e-10 in normalized frequency"))
        Xa=root_data[N]["A"]["X"];Xb=root_data[N]["B"]["X"]
        slope=slope_branch(*Xa,er,N)-slope_branch(*Xb,er,N)
        original=(cert["crossing"]["DeltaJ_derivative_interval"] if N==21
                  else cert["slope_interval"])
        rows.append(dict(N=N,check="dDeltaJ_depsilon",branch="A-B",
            case="full_bracket",source="certificate JSON",
            original_interval=original,independent_interval=concise(slope),
            agreement="PASS" if slope.lo>0 else "FAIL",
            notes="positive interval; enclosure widths may differ"))
        for label,eval_eps in (("lower",eps[N][0]),("upper",eps[N][1])):
            qa=root_check(N,centers[N]["A"],I(eval_eps))
            qb=root_check(N,centers[N]["B"],I(eval_eps))
            diff=qa["J"]-qb["J"]
            if N==21:
                orig=next(v["DeltaJ_interval"] for v in cert["crossing"]["endpoints"]
                          if v["endpoint"]==label)
            else:
                x=next(v for v in cert["crossing_checks"] if v["label"]==label)
                orig=str(x["DeltaJ"])
            ok=diff.hi<0 if label=="lower" else diff.lo>0
            rows.append(dict(N=N,check="crossing_endpoint_DeltaJ",
                branch="A-B",case=label,source="certificate JSON",
                original_interval=orig,independent_interval=concise(diff,31),
                agreement="PASS" if ok else "FAIL",
                notes="Taylor enclosure around root box"))
        with (ROOT/("interval_boxes.csv" if N==21 else "N31_interval_cells.csv")).open(newline="",encoding="utf-8") as f:
            cells=[r for r in csv.DictReader(f)
                   if r["status"]=="excluded" and
                   (N==31 or r["label"]=="crossing")]
        cells.sort(key=lambda r:Decimal(str(parse_original(r["interval_J_lower"]).lo)))
        hardest=cells[:10]
        coalescent=[r for r in cells if Decimal(r["h_lo"])<=Decimal("0.03")]
        coalescent.sort(key=lambda r:Decimal(str(parse_original(r["interval_J_lower"]).lo)))
        chosen=[("outer_min_margin",r) for r in hardest]
        chosen += [("coalescent",r) for r in coalescent[:3]]
        incumbent=min(root_data[N][k]["J"].hi for k in ("A","B"))
        for kind,r in chosen:
            bound=outer_bound(r,er,N)
            original=parse_original(r["interval_J_lower"])
            safety=bound.lo-incumbent
            difference=abs(bound.lo-original.lo)
            close=difference<safety/1000
            ok=bound.lo>incumbent and original.lo>incumbent and close
            idx=r.get("index",str(cells.index(r)))
            rows.append(dict(N=N,check=kind,branch="-",
                case=str(idx),source="outer-cell CSV",
                original_interval=concise(original,29),
                independent_interval=concise(bound,29),
                agreement="PASS" if ok else "FAIL",
                notes="lower-bound difference="+str(difference)+
                      "; difference<0.001*MPFR_margin="+str(close)+
                      "; m=["+r["m_lo"]+","+r["m_hi"]+
                      "]; h=["+r["h_lo"]+","+r["h_hi"]+
                      "]; independent incumbent upper="+str(incumbent)))
    out=ROOT/"independent_interval_spotcheck.csv"
    with out.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)
    print("rows",len(rows),"PASS",sum(r["agreement"]=="PASS" for r in rows),
          "FAIL",sum(r["agreement"]=="FAIL" for r in rows))
    for r in rows:
        if r["agreement"]!="PASS":
            print("FAILED",r["N"],r["check"],r["branch"],r["case"],
                  r["independent_interval"])
    if any(r["agreement"]!="PASS" for r in rows):
        raise SystemExit(1)


if __name__=="__main__":
    main()

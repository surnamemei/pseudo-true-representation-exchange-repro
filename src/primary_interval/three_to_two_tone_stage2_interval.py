"""Interval automatic differentiation for the N=21 real-phase Dirichlet objective.

Uses mpmath.iv (libmpi directed rounding). This module provides interval
values, gradients, and Hessians in normalized u=N*frequency coordinates.
"""

from __future__ import annotations

from dataclasses import dataclass
from mpmath import iv

N=21
D_STRONG=2
B_WEAK=10


def I(x):
    return x if isinstance(x,iv.mpf) else iv.mpf(x)


@dataclass
class Jet:
    v: object
    g: tuple
    H: tuple

    @staticmethod
    def c(x):
        z=iv.mpf(0)
        return Jet(I(x),(z,z),((z,z),(z,z)))

    @staticmethod
    def var(x,k):
        z=iv.mpf(0);o=iv.mpf(1)
        return Jet(I(x),(o if k==0 else z,o if k==1 else z),((z,z),(z,z)))

    def __add__(a,b):
        b=b if isinstance(b,Jet) else Jet.c(b)
        return Jet(a.v+b.v,tuple(a.g[i]+b.g[i] for i in range(2)),
                   tuple(tuple(a.H[i][j]+b.H[i][j] for j in range(2)) for i in range(2)))
    __radd__=__add__

    def __neg__(a):
        return Jet(-a.v,tuple(-v for v in a.g),tuple(tuple(-v for v in row) for row in a.H))

    def __sub__(a,b):return a+(-b if isinstance(b,Jet) else -Jet.c(b))
    def __rsub__(a,b):return (-a)+b

    def __mul__(a,b):
        b=b if isinstance(b,Jet) else Jet.c(b)
        return Jet(a.v*b.v,
                   tuple(a.g[i]*b.v+a.v*b.g[i] for i in range(2)),
                   tuple(tuple(a.H[i][j]*b.v+a.v*b.H[i][j]
                               +a.g[i]*b.g[j]+b.g[i]*a.g[j]
                               for j in range(2)) for i in range(2)))
    __rmul__=__mul__

    def reciprocal(a):
        v=1/a.v
        return Jet(v,tuple(-a.g[i]*v*v for i in range(2)),
                   tuple(tuple(2*a.g[i]*a.g[j]*v*v*v-a.H[i][j]*v*v
                               for j in range(2)) for i in range(2)))

    def __truediv__(a,b):
        b=b if isinstance(b,Jet) else Jet.c(b)
        return a*b.reciprocal()
    def __rtruediv__(a,b):return Jet.c(b)*a.reciprocal()

    def cos(a):
        c=iv.cos(a.v);s=iv.sin(a.v)
        return Jet(c,tuple(-s*a.g[i] for i in range(2)),
                   tuple(tuple(-s*a.H[i][j]-c*a.g[i]*a.g[j]
                               for j in range(2)) for i in range(2)))

    def sin(a):
        c=iv.cos(a.v);s=iv.sin(a.v)
        return Jet(s,tuple(c*a.g[i] for i in range(2)),
                   tuple(tuple(c*a.H[i][j]-s*a.g[i]*a.g[j]
                               for j in range(2)) for i in range(2)))


def dirichlet(z):
    # Exact finite trigonometric sum; no removable-singularity division.
    out=Jet.c(1)
    for k in range(1,(N-1)//2+1):
        out=out+2*(z*iv.mpf(k)/N).cos()
    return out


def data_energy(e):
    d=iv.mpf(D_STRONG);b=iv.mpf(B_WEAK)
    out=iv.mpf(0)
    for t in range(-(N-1)//2,(N-1)//2+1):
        s=iv.mpf(t)/N
        c=iv.cos(d*s)
        out+=4*c*c+e*e-4*e*c*iv.cos(b*s)
    return out


def objective_jet(u1,u2,e):
    u1=Jet.var(u1,0);u2=Jet.var(u2,1)
    e=I(e)
    h1=dirichlet(u1+D_STRONG)+dirichlet(u1-D_STRONG)-dirichlet(u1-B_WEAK)*e
    h2=dirichlet(u2+D_STRONG)+dirichlet(u2-D_STRONG)-dirichlet(u2-B_WEAK)*e
    cross=dirichlet(u2-u1)
    den=N*N-cross*cross
    numer=N*(h1*h1+h2*h2)-2*cross*h1*h2
    return Jet.c(data_energy(e))-numer/den


def objective_gradient_interval(u1,u2,e):
    """Cheaper direct interval J and first derivatives (no Hessian AD)."""
    e=I(e)
    def D0(z):
        return iv.mpf(1)+sum(2*iv.cos(iv.mpf(k)*z/N) for k in range(1,(N-1)//2+1))
    def D1(z):
        return sum(-2*iv.mpf(k)*iv.sin(iv.mpf(k)*z/N)/N
                   for k in range(1,(N-1)//2+1))
    g=D0(u2-u1);gp=D1(u2-u1)
    h1=D0(u1+D_STRONG)+D0(u1-D_STRONG)-e*D0(u1-B_WEAK)
    h2=D0(u2+D_STRONG)+D0(u2-D_STRONG)-e*D0(u2-B_WEAK)
    hp1=D1(u1+D_STRONG)+D1(u1-D_STRONG)-e*D1(u1-B_WEAK)
    hp2=D1(u2+D_STRONG)+D1(u2-D_STRONG)-e*D1(u2-B_WEAK)
    den=N*N-g*g
    cap=N*(h1*h1+h2*h2)-2*g*h1*h2
    cap1=2*N*h1*hp1+2*gp*h1*h2-2*g*hp1*h2
    cap2=2*N*h2*hp2-2*gp*h1*h2-2*g*h1*hp2
    den1=2*g*gp
    den2=-2*g*gp
    J=data_energy(e)-cap/den
    grad1=-(cap1*den-cap*den1)/(den*den)
    grad2=-(cap2*den-cap*den2)/(den*den)
    return J,(grad1,grad2)


def lo(x):return x.a
def hi(x):return x.b
def width(x):return x.b-x.a
def mid(x):return x.mid


def interval_str(x,digits=35):
    return f'[{str(x.a)},{str(x.b)}]'

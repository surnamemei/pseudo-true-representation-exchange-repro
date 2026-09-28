"""Independent Arb ball arithmetic for finite-record certificate replay.
Analytic formulas are implemented separately from the mpmath proof backend.
"""
import sys,json,csv,math,time
from pathlib import Path
from decimal import Decimal,getcontext
sys.path.insert(0,str(Path(__file__).resolve().parent/'stage12_deps'))
import flint
from flint import arb,ctx
ctx.prec=320
getcontext().prec=125
class I:
    __slots__=('x',)
    def __init__(self,lo=0,hi=None):
        if isinstance(lo,I): self.x=lo.x;return
        a=lo if isinstance(lo,arb) else arb(str(lo))
        self.x=a if hi is None else a.union(hi.x if isinstance(hi,I) else arb(str(hi)))
    @property
    def lo(self):return self.x.lower()
    @property
    def hi(self):return self.x.upper()
    def mid(self):return self.x.mid()
    def absmax(self):return self.x.abs_upper()
    def __str__(self):return str(self.x)
    def __add__(a,b):return I(a.x+I(b).x)
    __radd__=__add__
    def __neg__(a):return I(-a.x)
    def __sub__(a,b):return I(a.x-I(b).x)
    def __rsub__(a,b):return I(I(b).x-a.x)
    def __mul__(a,b):return I(a.x*I(b).x)
    __rmul__=__mul__
    def __truediv__(a,b):return I(a.x/I(b).x)
    def __rtruediv__(a,b):return I(I(b).x/a.x)
    def __pow__(a,n):return I(a.x**n)
    def reciprocal(a):return I(1/a.x)
    def cos(a):return I(a.x.cos())
    def sin(a):return I(a.x.sin())
    def sqrt(a):return I(a.x.sqrt())
Z=I(0);O=I(1)
def box(mid,r):return I(I(mid).x+arb(0,str(r)))
def bounds(x):return x.x.str(95,more=True)
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



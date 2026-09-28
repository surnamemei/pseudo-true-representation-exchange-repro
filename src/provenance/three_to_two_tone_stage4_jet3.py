"""Third-order bivariate interval Taylor coefficients for the finite-record cost.

Used only for centered-cell remainder bounds. Coefficients are derivative / alpha!.
"""
from mpmath import iv
import three_to_two_tone_stage2_interval as si

MON=[(i,j) for i in range(4) for j in range(4-i)]
ZERO=iv.mpf(0);ONE=iv.mpf(1)

class P:
    def __init__(self,d=None):self.d={} if d is None else d
    @staticmethod
    def c(v):return P({(0,0):iv.mpf(v)})
    @staticmethod
    def var(v,k):
        q=P.c(v);q.d[(1,0) if k==0 else (0,1)]=ONE;return q
    def get(self,k):return self.d.get(k,ZERO)
    def __add__(self,b):
        b=asP(b);return P({k:self.get(k)+b.get(k) for k in MON})
    __radd__=__add__
    def __neg__(self):return P({k:-v for k,v in self.d.items()})
    def __sub__(self,b):return self+(-asP(b))
    def __rsub__(self,b):return asP(b)+(-self)
    def __mul__(self,b):
        b=asP(b);out={}
        for (i,j),v in self.d.items():
            for (k,l),w in b.d.items():
                if i+j+k+l<=3:
                    key=(i+k,j+l);out[key]=out.get(key,ZERO)+v*w
        return P(out)
    __rmul__=__mul__
    def recip(self):
        a=self.get((0,0));r=ONE/a
        v=self-P.c(a)
        t=v*(-r)
        return P.c(r)*(P.c(1)+t+t*t+t*t*t)
    def __truediv__(self,b):return self*asP(b).recip()
    def __rtruediv__(self,b):return asP(b)*self.recip()
    def cos(self):
        a=self.get((0,0));v=self-P.c(a)
        return P.c(iv.cos(a))*(P.c(1)-v*v/2)+P.c(-iv.sin(a))*(v-v*v*v/6)

def asP(x):return x if isinstance(x,P) else P.c(x)

def D(x):return P.c(1)+sum(2*(x*iv.mpf(k)/21).cos() for k in range(1,11))

def objective_coeffs(u1,u2,e,z):
    a=P.var(u1,0);b=P.var(u2,1)
    h1=D(a+z)+D(a-z)-D(a-10)*e
    h2=D(b+z)+D(b-z)-D(b-10)*e
    g=D(b-a);den=21*21-g*g
    cap=21*(h1*h1+h2*h2)-2*g*h1*h2
    return P.c(si.data_energy(e))-cap/den

def max_third_remainder(coeffs,w1,w2):
    # Sum of |c_alpha|*w^alpha directly avoids a crude isotropic M3.
    return sum(max(abs(coeffs.get((i,j)).a),abs(coeffs.get((i,j)).b))*w1**i*w2**j
               for i,j in ((3,0),(2,1),(1,2),(0,3)))

"""Exact-rational majorant for the centered-kernel H-series tail."""
from fractions import Fraction
import sys

def upper_error(n,qmax=110,kmax=46):
    # H(w)=sqrt(w)/sin(sqrt(w)). On |w|=4, |H(w)|<6, since
    # sinh(2)/2<11/6. Cauchy gives |h_m|<6/4^m.
    # |d^(j)(q)|<=2^-j for real q from the integral representation.
    n=int(n);qmax=int(qmax)
    r=Fraction(qmax*qmax,16*n*n)
    assert r<1
    total=Fraction(0)
    for m in range(9,101):
        base=Fraction(1,2)+Fraction(2*m,qmax)
        total+=r**m*max(Fraction(1),base**kmax)
    m=101
    base=Fraction(1,2)+Fraction(2*m,qmax)
    ratio=r*(Fraction(m+1,m)**kmax)
    assert ratio<1
    total+=r**m*max(Fraction(1),base**kmax)/(1-ratio)
    return 6*total

def decimal_power_ceiling(x):
    """Return 10^e at least x, with exact integer comparisons."""
    e=0
    while Fraction(10)**e<x:e+=1
    while Fraction(10)**(e-1)>=x:e-=1
    return e

if __name__=="__main__":
    for s in sys.argv[1:]:
        n=int(s);v=upper_error(n)
        print(n,float(v),decimal_power_ceiling(v))

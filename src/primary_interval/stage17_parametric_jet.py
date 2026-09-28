"""Exact polynomial-in-epsilon Jet for the N=21 projected loss."""
from mpmath import iv
import three_to_two_tone_stage2_interval as it

it.N=21
iv.dps=80
J=it.Jet
E0=it.data_energy(iv.mpf(0))
E1=(it.data_energy(iv.mpf(1))-it.data_energy(iv.mpf(-1)))/2
E2=iv.mpf(21)

def coefficients(u1,u2):
    u1=J.var(u1,0);u2=J.var(u2,1)
    p1=it.dirichlet(u1+2)+it.dirichlet(u1-2)
    p2=it.dirichlet(u2+2)+it.dirichlet(u2-2)
    z1=it.dirichlet(u1-10)
    z2=it.dirichlet(u2-10)
    g=it.dirichlet(u2-u1)
    den=21*21-g*g
    c=J.c(E0)-(21*(p1*p1+p2*p2)-2*g*p1*p2)/den
    b=J.c(E1)-(-2*21*(p1*z1+p2*z2)+2*g*(p1*z2+p2*z1))/den
    a=J.c(E2)-(21*(z1*z1+z2*z2)-2*g*z1*z2)/den
    return a,b,c

def objective(u1,u2,e):
    a,b,c=coefficients(u1,u2)
    return c+b*e+a*e*e

"""Independent 50-digit Dirichlet-kernel crossing check for the baseline."""

from pathlib import Path
import mpmath as mp

mp.mp.dps=50
n=21
d=mp.mpf(2)
b=mp.mpf(10)
ss=[mp.mpf(k)/n for k in range(-10,11)]


def D(u):
    return mp.fsum(mp.cos(u*s) for s in ss)


def objective(u1,u2,e):
    g=D(u2-u1)
    h1=D(u1+d)+D(u1-d)-e*D(u1-b)
    h2=D(u2+d)+D(u2-d)-e*D(u2-b)
    norm=mp.fsum(4*mp.cos(d*s)**2+e**2-4*e*mp.cos(d*s)*mp.cos(b*s)
                 for s in ss)
    return norm-(n*(h1*h1+h2*h2)-2*g*h1*h2)/(n*n-g*g)


def equations(a1,a2,c1,c2,e):
    return (mp.diff(lambda q:objective(q,a2,e),a1),
            mp.diff(lambda q:objective(a1,q,e),a2),
            mp.diff(lambda q:objective(q,c2,e),c1),
            mp.diff(lambda q:objective(c1,q,e),c2),
            objective(a1,a2,e)-objective(c1,c2,e))


sol=mp.findroot(equations,(-4.77244734,.60517033,-.05765356,12.46341936,.24819063),
                tol=mp.mpf('1e-43'),maxsteps=35,solver='mdnewton')
a1,a2,c1,c2,e=sol
lines=[f'a1={mp.nstr(a1,48)}',f'a2={mp.nstr(a2,48)}',
       f'b1={mp.nstr(c1,48)}',f'b2={mp.nstr(c2,48)}',
       f'epsilon_cross={mp.nstr(e,48)}',
       f'J_A={mp.nstr(objective(a1,a2,e),48)}',
       f'J_B={mp.nstr(objective(c1,c2,e),48)}',
       f'J_A_minus_J_B={mp.nstr(objective(a1,a2,e)-objective(c1,c2,e),10)}']
Path(__file__).with_name('high_precision_crossing.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print('\n'.join(lines))

"""Independent formula cross-checks for epsilon-quadratic interval coefficients."""
import json
import sys
from decimal import Decimal
from pathlib import Path
from mpmath import iv

sys.argv=['stage17_epsilon_quadratic_trial.py','0.0025','0']
import stage17_epsilon_quadratic_trial as q
import stage17_parametric_jet as p
import three_to_two_tone_stage2_interval as it
import three_to_two_tone_stage2_outer_global as outer

iv.dps=80;it.N=21
indices=[0,1,2,10,100,500,2000,4000,8000,12000,14286]
eps=['0.23','0.2481906301722774','0.27']
def overlaps(a,b):return not (a.b<b.a or b.b<a.a)
count_outer=0
for k in indices:
    row=q.rows[k]
    c,_=q.cell_coeff(row)
    aa,bb,cc,dd=(Decimal(row[n]) for n in ('m_lo','m_hi','h_lo','h_hi'))
    m=iv.mpf(str((aa+bb)/2));h=iv.mpf(str((cc+dd)/2))
    for e0 in eps:
        e=iv.mpf(e0)
        a=q.val(c,e)
        b=outer.center_projected_J(m,h,e,it.data_energy(e))[0]
        assert overlaps(a,b),(k,e0,a,b)
        count_outer+=1
count_jet=0
for u1,u2 in [('-4.7724473539','0.6051703203'),('-0.0576535619','12.4634193568'),
              ('-2','2'),('-3.2','1.1')]:
    for e0 in eps:
        a=p.objective(iv.mpf(u1),iv.mpf(u2),iv.mpf(e0))
        b=it.objective_jet(iv.mpf(u1),iv.mpf(u2),iv.mpf(e0))
        vals=[(a.v,b.v)]+[(a.g[i],b.g[i]) for i in range(2)]
        vals += [(a.H[i][j],b.H[i][j]) for i in range(2) for j in range(2)]
        assert all(overlaps(x,y) for x,y in vals),(u1,u2,e0)
        count_jet+=1
out=dict(outer_center_value_pairs=count_outer,parametric_jet_pairs=count_jet,
         all_interval_formulas_overlap=True)
(Path(__file__).resolve().parent/'stage17_epsilon_coeff_crosscheck.json').write_text(json.dumps(out,indent=2))
print(out)

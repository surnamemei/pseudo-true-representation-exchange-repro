"""One bounded adaptive-subdivision pilot for local crossing validation.

This does not prove globality. It diagnoses raw finite-frequency interval
dependency before committing to a full parameterized full-domain cover.
"""
import csv
from mpmath import mp
from three_to_two_tone_stage4_interval_curve import validate

mp.dps=75
rows=[]
for center in ('0.1','0.5','1.0','1.5','1.9'):
    c=mp.mpf(center)
    for width in ('0.01','0.001','0.0001','0.00001'):
        w=mp.mpf(width)
        q=validate(c-w/2,c+w/2,inflate=2)
        rows.append(dict(z_mid=center,z_width=width,
                         p_lo=mp.nstr((c-w/2)**2,18),p_hi=mp.nstr((c+w/2)**2,18),
                         verified=q['verified'],included=q['included'],
                         contraction_upper=q.get('contraction_upper'),
                         A_Hdet_lower=q.get('A_Hdet_lower'),
                         B_Hdet_lower=q.get('B_Hdet_lower'),
                         slope_lower=q.get('slope_lower'),
                         A_gram_lower=q.get('A_gram_lower'),
                         B_gram_lower=q.get('B_gram_lower'),
                         max_K_radius_ratio=max(q.get('K_radius_ratio',[float('nan')])),
                         status='local_only' if q['verified'] else 'failed_local_enclosure'))
        print(center,width,q['verified'],q.get('contraction_upper'),flush=True)
for center,widths in (('0.1',('0.0000001','0.00000001','0.000000001')),
                      ('0.5',('0.000001','0.0000001')),
                      ('1.0',('0.000001','0.0000001'))):
    c=mp.mpf(center)
    for width in widths:
        w=mp.mpf(width)
        q=validate(c-w/2,c+w/2,inflate=2)
        rows.append(dict(z_mid=center,z_width=width,
                         p_lo=mp.nstr((c-w/2)**2,18),p_hi=mp.nstr((c+w/2)**2,18),
                         verified=q['verified'],included=q['included'],
                         contraction_upper=q.get('contraction_upper'),
                         A_Hdet_lower=q.get('A_Hdet_lower'),
                         B_Hdet_lower=q.get('B_Hdet_lower'),
                         slope_lower=q.get('slope_lower'),
                         A_gram_lower=q.get('A_gram_lower'),
                         B_gram_lower=q.get('B_gram_lower'),
                         max_K_radius_ratio=max(q.get('K_radius_ratio',[float('nan')])),
                         status='local_only' if q['verified'] else 'failed_local_enclosure'))
        print(center,width,q['verified'],q.get('contraction_upper'),flush=True)
with open('stage18_local_slab_pilot.csv','w',newline='') as f:
    wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)

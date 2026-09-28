"""Refine only archived z=2 outer cells limiting a wider parameter transfer."""
import csv,json,math
from pathlib import Path
from mpmath import iv,mp
from three_to_two_tone_stage2_outer_global import verify_cell
from three_to_two_tone_stage2_interval import data_energy

ROOT=Path(__file__).resolve().parent
iv.dps=65;mp.dps=75
DZ=mp.mpf('1e-4');DE=mp.mpf('5e-5')
ETA=mp.mpf('5e-4')
assert ETA>2*mp.sqrt(mp.mpf(440)/252)*DZ+mp.sqrt(21)*DE
ref=json.loads((ROOT/'stage2_reference.json').read_text())
e0=mp.mpf(ref['epsilon_cross'])
e=iv.mpf(str(e0))+iv.mpf(['-1e-50','1e-50'])
E=data_energy(e)
outrec=next(r for r in json.loads((ROOT/'stage2_outer_global_summary.json').read_text())
            if r['label']=='crossing')
U=mp.mpf(outrec['incumbent_upper'].strip('[]').split(',')[-1])
Utransfer=(mp.sqrt(U)+ETA)**2

def transferred_margin(L):
    return max(mp.mpf(0),mp.sqrt(max(L,mp.mpf(0)))-ETA)**2-Utransfer

def parseL(r):return mp.mpf(r['interval_J_lower'].strip('[]').split(',')[0])

def refine(r,maxdepth=8):
    orig=(float(r['m_lo']),float(r['m_hi']),float(r['h_lo']),float(r['h_hi']))
    stack=[(orig,0)];leaf=[]
    while stack:
        (m0,m1,h0,h1),dep=stack.pop()
        cell=(int(r['depth'])+2*dep,m0,m1,h0,h1,0.,0.,'excluded')
        ok,B,why=verify_cell(cell,e,E,U,strong_d=iv.mpf(2))
        L=mp.mpf(str(B.a).strip('[]').split(',')[0]) if B is not None else mp.mpf(0)
        margin=transferred_margin(L)
        if margin>0:
            leaf.append(dict(depth=dep,m_lo=m0,m_hi=m1,h_lo=h0,h_hi=h1,
                             lower=str(L),transferred_margin=str(margin),status='excluded'))
        elif dep>=maxdepth:
            leaf.append(dict(depth=dep,m_lo=m0,m_hi=m1,h_lo=h0,h_hi=h1,
                             lower=str(L),transferred_margin=str(margin),status='unresolved',reason=why))
        else:
            mm=(m0+m1)/2;hh=(h0+h1)/2
            stack.extend([((m0,mm,h0,hh),dep+1),((mm,m1,h0,hh),dep+1),
                          ((m0,mm,hh,h1),dep+1),((mm,m1,hh,h1),dep+1)])
    return leaf

def main():
    rows=[r for r in csv.DictReader((ROOT/'interval_boxes.csv').open())
          if r['label']=='crossing' and r['status']=='excluded']
    bad=[(i,r) for i,r in enumerate(rows) if transferred_margin(parseL(r))<=0]
    refined=[];fail=[]
    for n,(i,r) in enumerate(bad):
        leaves=refine(r)
        refined.append(dict(original_cell_index=i,original_lower=str(parseL(r)),
                            original_margin=str(transferred_margin(parseL(r))),
                            children=len(leaves),unresolved=sum(q['status']=='unresolved' for q in leaves),
                            min_margin=min(float(q['transferred_margin']) for q in leaves),
                            leaves=leaves))
        if refined[-1]['unresolved']:fail.append(i)
        if n%10==0:print(n,'of',len(bad),'children',len(leaves),'unresolved',refined[-1]['unresolved'],flush=True)
    untouched=[transferred_margin(parseL(r)) for r in rows
               if transferred_margin(parseL(r))>0]
    result=dict(z_interval=['1.9999','2.0'],eta=str(ETA),epsilon_halfwidth=str(DE),
                original_outer_cells=len(rows),original_failing_cells=len(bad),
                refined_children=sum(r['children'] for r in refined),
                remaining_original_failures=fail,
                min_untouched_margin=str(min(untouched)),
                min_refined_margin=min((r['min_margin'] for r in refined),default=None),
                outer_uniform_exclusion=not fail,
                note='outer full-domain cells only; inner cells and local roots are separate obligations',
                refined=refined)
    (ROOT/'stage4_targeted_transfer_outer.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:result[k] for k in result if k!='refined'},indent=2),flush=True)

if __name__=='__main__':main()

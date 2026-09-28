"""Recheck archived inner cells over the ten existing local parameter boxes."""
import csv,json,math
from pathlib import Path
from mpmath import iv,mp
import three_to_two_tone_stage2_interval as si
from three_to_two_tone_stage4_interval_curve import root
from three_to_two_tone_stage4_global_interval import recheck_nonobjective

ROOT=Path(__file__).resolve().parent
iv.dps=55;mp.dps=75

def end(x,side=0):return mp.mpf(str(x).strip('[]').split(',')[side])
def run():
    boxes=json.loads((ROOT/'stage4_interval_boxes.json').read_text())
    cells=[r for r in csv.DictReader((ROOT/'stage2_inner_cells.csv').open())
           if r['label']=='crossing']
    ans=[]
    for i,q in enumerate(boxes):
        zlo=mp.mpf(q['z_lo']);zhi=mp.mpf(q['z_hi']);zm=(zlo+zhi)/2
        assert q['verified']
        er=q['X_bounds'][4]
        elo=end(er[0],0)-mp.mpf('1e-18')
        ehi=end(er[1],1)+mp.mpf('1e-18')
        em=root(zm)[4]
        dz=(zhi-zlo)/2;de=max(abs(elo-em),abs(ehi-em))
        eta=mp.mpf('2.644')*dz+mp.mpf('4.583')*de
        si.D_STRONG=iv.mpf(mp.nstr(zm,60))
        e0=iv.mpf(mp.nstr(em,60))
        a,b,c,d,_=root(zm)
        ja=si.objective_jet(iv.mpf(mp.nstr(a,65)),iv.mpf(mp.nstr(b,65)),e0).v.b
        jb=si.objective_jet(iv.mpf(mp.nstr(c,65)),iv.mpf(mp.nstr(d,65)),e0).v.b
        U0=min(end(ja),end(jb))
        U=(mp.sqrt(U0)+eta)**2
        eb=iv.mpf([mp.nstr(elo,65),mp.nstr(ehi,65)])
        db=iv.mpf([mp.nstr(zlo,65),mp.nstr(zhi,65)])
        _,fails,_=recheck_nonobjective(cells,eb,db,U,mp.mpf('5e-4'))
        ans.append(dict(i=i,z_lo=str(zlo),z_hi=str(zhi),e_lo=str(elo),e_hi=str(ehi),
                        candidate_upper=str(U),eta_local=str(eta),
                        failures=len(fails),failure_cells=fails))
        print(i,str(zlo),str(zhi),'failures',len(fails),flush=True)
    (ROOT/'stage4_targeted_transfer_inner_boxes.json').write_text(json.dumps(ans,indent=2))
    return ans

if __name__=='__main__':run()

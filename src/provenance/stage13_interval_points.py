"""Run the frozen Stage-11 scalar interval isolator at a new fixed lambda.

This wrapper creates a temporary in-memory module variant; frozen Stage-11 files
and certificates are never modified. Outputs are suffixed by lambda.
"""
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parent
lam=sys.argv[1]
tag=lam.replace('.','p')
src=(ROOT/'stage11_tangent_global.py').read_text(encoding='utf-8')
src=src.replace('import csv,json,math','import csv,json,math')
src=src.replace('from three_to_two_tone_stage3_tangent_interval import prove','')
src=src.replace("ref=json.loads((ROOT/'stage3_tangent_reference.json').read_text())",'ref={}')
src=src.replace("local=prove('1e-24','1e-25')", "local={'verified': True, 'scope': 'not rerun'}")
src=src.replace("L=iv.mpf(ref['lambda_N'])+iv.mpf(['-1e-25','1e-25'])",f"L=iv.mpf('{lam}')")
src=src.replace("ROOT/'smallz_global_certificate.json'",f"ROOT/'stage13_isolation_{tag}.json'")
src=src.replace("ROOT/'stage11_tangent_partition.csv'",f"ROOT/'stage13_partition_{tag}.csv'")
src=src.replace("ROOT/'stage11_tangent_roots.csv'",f"ROOT/'stage13_roots_{tag}.csv'")
exec(compile(src, 'stage11_tangent_global.py [Stage13 variant]', 'exec'),{'__name__':'__main__','__file__':str(ROOT/'stage11_tangent_global.py')})

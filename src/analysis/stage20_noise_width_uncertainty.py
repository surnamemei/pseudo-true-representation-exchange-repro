"""Delta-method uncertainty for archived 10--90% probit transition widths."""
import csv
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm

root=Path(__file__).resolve().parent
counts=pd.read_csv(root/'branch_probability.csv')
widths=pd.read_csv(root/'transition_width.csv')
out=[]
for row in widths.itertuples():
    g=counts[counts.snr_db==row.snr_db]
    eta=g.eta.to_numpy()
    linear=row.fitted_probit_intercept+row.fitted_probit_slope*eta
    p=np.clip(norm.cdf(linear),1e-12,1-1e-12)
    weight=g.n_trials.to_numpy()*norm.pdf(linear)**2/(p*(1-p))
    design=np.column_stack((np.ones(len(eta)),eta))
    cov=np.linalg.inv(design.T@(weight[:,None]*design))
    width=row.transition_width_eta_10_90
    se=width*np.sqrt(cov[1,1])/abs(row.fitted_probit_slope)
    out.append(dict(snr_db=int(row.snr_db),n_settings=len(g),trials_per_setting=int(g.n_trials.iloc[0]),
                    width_eta_10_90=width,delta_method_se=se,
                    approximate_95_low=width-1.96*se,approximate_95_high=width+1.96*se,
                    method='inverse expected Fisher information of binomial probit fit; delta method'))
with (root/'noise_transition_width_uncertainty.csv').open('w',newline='',encoding='utf-8') as stream:
    writer=csv.DictWriter(stream,fieldnames=list(out[0]));writer.writeheader();writer.writerows(out)
print(out)

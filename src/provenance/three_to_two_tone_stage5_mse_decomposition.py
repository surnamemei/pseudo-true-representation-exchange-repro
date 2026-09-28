"""Add the reference-dependent bias term to the exact mixture MSE identity."""
from pathlib import Path
import numpy as np
import pandas as pd

root=Path(__file__).resolve().parent
mix=pd.read_csv(root/'mixture_statistics.csv')
trials=pd.read_csv(root/'noise_branch_trials.csv')
target=np.array([-2.,2.])
for ix,row in mix.iterrows():
    d=trials[(trials.snr_db==row.snr_db)&(trials.eta==0)]
    u=d[['u1','u2']].to_numpy()
    bias_sq=np.sum((np.mean(u,axis=0)-target)**2)
    mse=np.mean(np.sum((u-target)**2,axis=1))
    mix.loc[ix,'bias_to_true_strong_sq']=bias_sq
    mix.loc[ix,'MSE_to_true_strong']=mse
    mix.loc[ix,'within_plus_between_plus_bias']=row.within_variance_trace+row.between_selection_variance_trace+bias_sq
    mix.loc[ix,'MSE_identity_error']=mse-mix.loc[ix,'within_plus_between_plus_bias']
mix.to_csv(root/'mixture_statistics.csv',index=False)
print(mix[['snr_db','within_variance_trace','between_selection_variance_trace',
           'bias_to_true_strong_sq','MSE_to_true_strong','MSE_identity_error']].to_string(index=False))

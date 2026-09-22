"""Replay all trained networks from bundled inputs; no raw GWAS/genotypes needed."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import torch
HERE=Path(__file__).resolve().parent; sys.path.insert(0,str(HERE))
from revision_guardnet_validation import m,operations
out=HERE.parent/'revision_validation'
meta=pd.read_csv(out/'challenge_metadata.tsv',sep='\t')
z1=np.loadtxt(out/'challenge_z1.tsv',delimiter='\t'); z2=np.loadtxt(out/'challenge_z2.tsv',delimiter='\t')
matrices={p.stem:np.loadtxt(p,delimiter='\t').astype('float32') for p in sorted((out/'ld_inputs').glob('*.tsv'))}
m.PANELS=tuple(matrices); names=list(matrices)
split=dict(x=np.stack([m.node_features(a,b) for a,b in zip(z1,z2)]),
    panel=np.array([names.index(name) for name in meta.matrix]),y=(meta.fault!='valid').to_numpy().astype(int))
saved=pd.read_csv(out/'all_challenge_predictions.tsv',sep='\t')
torch.set_num_threads(4)
for mode in ('signed_graph','pooled_mlp','unsigned_graph'):
    pos,neg=operations(matrices,mode); prob=[]
    for seed in (11,29,47):
        path=(HERE.parent/f'ld_guardnet/ld_guardnet_seed_{seed}.pt') if mode=='signed_graph' else out/f'{mode}_{seed}.pt'
        state=torch.load(path,weights_only=True); model=m.SignedLDGuardNet(5,48,.2); model.load_state_dict(state['state_dict'])
        prob.append(m.predict_network(model,split,pos,neg,128))
    actual=1-np.mean(prob,axis=0)[:,0]
    difference=np.max(abs(actual-saved[mode+'_anomaly_score'].to_numpy()))
    assert difference<2e-6,(mode,difference)
    print(f'{mode}: {len(actual)} scores reproduced; max_abs_difference={difference:.9g}')

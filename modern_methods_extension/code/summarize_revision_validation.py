"""Complete non-selective challenge summary, including the named diagnostic."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
out=Path(__file__).resolve().parent.parent/'revision_validation'
df=pd.read_csv(out/'challenge_predictions.tsv',sep='\t')
valid=df.fault.eq('valid').to_numpy(); y=(~valid).astype(int)
val=pd.read_csv(out/'validation_kriging_scores.tsv',sep='\t')
test=pd.read_csv(out/'challenge_kriging_scores.tsv',sep='\t')
assert len(test)==len(df) and np.array_equal(test.row_index.to_numpy(),np.arange(len(df))+1)
thresholds=json.loads((out/'validation_thresholds.json').read_text())
thresholds['susie_kriging']=float(np.quantile(val.score.dropna(),.95))
df['susie_kriging_anomaly_score']=test.score
models=list(thresholds); summaries=[]; strata=[]
groups=list(df.groupby('region',sort=True).indices.values())
rng=np.random.default_rng(20260921)
indices=[np.concatenate([groups[j] for j in rng.integers(0,len(groups),len(groups))]) for _ in range(1000)]
for model in models:
    score=df[model+'_anomaly_score'].to_numpy(); finite=np.isfinite(score); pred=score>thresholds[model]
    auc=roc_auc_score(y[finite],score[finite])
    boot=[roc_auc_score(y[ix][finite[ix]],score[ix][finite[ix]]) for ix in indices]
    ci=np.quantile(boot,[.025,.975])
    summaries.append(dict(model=model,anomaly_auc=auc,region_bootstrap_low=ci[0],region_bootstrap_high=ci[1],
        valid_false_positives=int(pred[valid].sum()),valid_n=int(valid.sum()),valid_fpr=float(pred[valid].mean()),
        threshold=thresholds[model],failed_scores=int((~finite).sum())))
    for key,ix in df.groupby(['region','fault']).indices.items():
        strata.append(dict(model=model,region=key[0],fault=key[1],n=len(ix),flagged=int(pred[ix].sum()),flag_rate=float(pred[ix].mean())))
df.to_csv(out/'all_challenge_predictions.tsv',sep='\t',index=False)
pd.DataFrame(summaries).to_csv(out/'all_model_challenge_summary.tsv',sep='\t',index=False)
pd.DataFrame(strata).to_csv(out/'all_model_region_fault_metrics.tsv',sep='\t',index=False)
(out/'all_validation_thresholds.json').write_text(json.dumps(thresholds,indent=2),encoding='utf-8')
print(pd.DataFrame(summaries).to_string(index=False))
print('kriging mean flag rates by fault:')
print(pd.DataFrame(strata).query("model=='susie_kriging'").groupby('fault').flag_rate.mean().to_string())

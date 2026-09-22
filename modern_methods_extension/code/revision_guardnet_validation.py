"""Review-triggered, non-selective challenge and controlled ablations.

The original benchmark and frozen checkpoints are unchanged. All new results,
including negative findings, are written to a separate revision directory.
"""
import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('original_guardnet', HERE/'run_ld_guardnet_benchmark.py')
m = importlib.util.module_from_spec(spec); sys.modules[spec.name] = m; spec.loader.exec_module(m)


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def bootstrap(frame, probability, baseline, repeats=2000):
    groups = list(frame.groupby('replicate', sort=True).indices.values())
    y = frame.true_class.to_numpy()
    rng = np.random.default_rng(7312026)
    rows = []
    for _ in range(repeats):
        ix = np.concatenate([groups[j] for j in rng.integers(0,len(groups),len(groups))])
        a = m.probability_metrics(y[ix], probability[ix])
        b = m.probability_metrics(y[ix], baseline[ix])
        rows.append([a['macro_ovr_auc'],a['macro_f1'],a['accuracy'],a['anomaly_auc'],
                     a['macro_ovr_auc']-b['macro_ovr_auc']])
    ci = np.quantile(rows,[.025,.975],axis=0)
    names = ['macro_ovr_auc','macro_f1','accuracy','anomaly_auc','paired_auc_gain']
    return {name:ci[:,i].tolist() for i,name in enumerate(names)}


def operations(matrices, mode):
    pos, neg = m.signed_operators(matrices)
    if mode == 'pooled_mlp':
        pos = torch.eye(pos.shape[1]).repeat(len(matrices),1,1); neg = pos.clone()
    elif mode == 'unsigned_graph':
        mats = []
        for r in matrices.values():
            a = np.abs(r).copy(); np.fill_diagonal(a,1)
            inv = 1/np.sqrt(a.sum(1)); mats.append((inv[:,None]*a*inv[None,:]).astype('float32'))
        pos = torch.from_numpy(np.stack(mats)); neg = pos.clone()
    return pos, neg


def split_with_raw(**kwargs):
    collector = []
    original = m.node_features
    def record(z1,z2):
        collector.append(np.stack([z1,z2])); return original(z1,z2)
    m.node_features = record
    try: result = m.build_split(**kwargs)
    finally: m.node_features = original
    result['raw'] = np.stack(collector)
    return result


def challenge(matrices, region_names, ids, out):
    faults = ['valid','random_permutation','block_permutation','sign_1pct','sign_5pct',
              'sign_10pct','sign_20pct','ld_scramble']
    xs, ys, panels, raw, records, baseline = [],[],[],[],[],[]
    # Lock all LD-content corruptions before creating any summary vectors.
    for index,name in enumerate(region_names):
        perm = np.random.default_rng(9941+index).permutation(100)
        matrices[name+'_scrambled'] = matrices[name][np.ix_(perm,perm)]
    panel_names = list(matrices)
    for ri,name in enumerate(region_names):
        for rep in range(40):
            arch = m.ARCHITECTURES[rep%4]; n = (2000,20000)[(rep//4)%2]
            rng = np.random.default_rng(990001+ri*100003+rep*107)
            c1,c2 = m.causal_configuration(40001+rep,arch,100)
            z1 = m.simulate_z(rng,matrices[name],n,c1,.085)
            z2 = m.simulate_z(rng,matrices[name],n,c2,.075)
            for fault in faults:
                b = z2.copy(); fit = name
                if fault=='random_permutation': b = b[rng.permutation(100)]
                elif fault=='block_permutation':
                    start = int(rng.integers(0,81)); b[start:start+20] = b[start:start+20][rng.permutation(20)]
                elif fault.startswith('sign_'):
                    k = int(fault.split('_')[1].replace('pct','')); b[rng.choice(100,k,replace=False)] *= -1
                elif fault=='ld_scramble': fit = name+'_scrambled'
                feat = m.node_features(z1,b)
                xs.append(feat); ys.append(0 if fault=='valid' else 1); panels.append(panel_names.index(fit))
                raw.append(np.stack([z1,b])); baseline.append(m.aggregate_baseline_features(feat,matrices[fit]))
                records.append(dict(region=name,replicate=rep,architecture=arch,sample_size=n,
                                    fault=fault,matrix=fit,base_id=f'{name}:{rep}'))
    return dict(x=np.stack(xs),y=np.array(ys),panel=np.array(panels),raw=np.stack(raw),
                baseline=np.stack(baseline),metadata=records)


def export_diagnostic(split, name, names, out):
    frame = pd.DataFrame(split['metadata'])
    frame['matrix'] = [names[i] for i in split['panel']]
    frame['is_valid'] = split['y']==0
    if 'fault' not in frame: frame['fault'] = frame['label']
    frame.to_csv(out/f'{name}_metadata.tsv',sep='\t',index=False)
    for trait in (0,1): np.savetxt(out/f'{name}_z{trait+1}.tsv',split['raw'][:,trait,:],delimiter='\t')


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--workspace',type=Path,required=True)
    args=parser.parse_args(); root=args.workspace
    out=HERE.parent/'revision_validation'; out.mkdir(exist_ok=True)
    torch.set_num_threads(4)
    config=m.Config()
    ids0, matrices, hashes=m.load_ld_matrices(root/'paper/revision_20260831_strengthened/ancestry_ld')
    # Cache identical eigen decompositions; draws and arithmetic are unchanged.
    cache={}
    def draw(rng,matrix):
        key=matrix.tobytes()
        if key not in cache: cache[key]=np.linalg.eigh(matrix.astype(np.float64))
        values,vectors=cache[key]
        return vectors@(np.sqrt(np.maximum(values,0))*rng.standard_normal(len(values)))
    m.draw_correlated=draw
    train=split_with_raw(matrices=matrices,replicate_start=1,replicates=220,data_panels=('EUR','EAS','AFR'),
        fit_panels=('EUR','EAS','AFR'),architectures=('one_shared','two_shared'),sample_sizes=(5000,),seed=config.seed)
    val=split_with_raw(matrices=matrices,replicate_start=10001,replicates=60,data_panels=('EUR','EAS','AFR'),
        fit_panels=('EUR','EAS','AFR'),architectures=('one_shared','two_shared'),sample_sizes=(5000,),seed=config.seed)
    oldtest=split_with_raw(matrices=matrices,replicate_start=20001,replicates=100,data_panels=('SAS','AMR'),
        fit_panels=('SAS','AMR'),architectures=('shared_plus_specific','distinct_only'),sample_sizes=(2000,20000),seed=config.seed)
    region_names=[]; regions=[]; ids={}
    for source in sorted((root/'paper/reanalysis_outputs/new_manuscript_20260828/multiregion_transfer/executed/regions').glob('R*/ld_mean_impute.npz')):
        a=np.load(source,allow_pickle=False); variants=a['variant_ids'].astype(str)
        stop = len(variants); ix=np.linspace(0,stop-1,100,dtype=int)
        assert len(set(variants[ix]))==100 and not (set(ids0)&set(variants[ix]))
        positions=[int(v.split('_')[1]) for v in variants[ix]]
        assert min(positions)>max(int(v.split('_')[1]) for v in ids0)
        r=m.nearest_correlation(a['ld'][np.ix_(ix,ix)]).astype('float32')
        name=source.parent.name; matrices[name]=r; region_names.append(name); ids[name]=variants[ix].tolist()
        regions.append(dict(region=name,source_relative=str(source.relative_to(root)),source_sha256=digest(source),
                            variants=100,min_position=min(positions),max_position=max(positions),source_variants=stop))
    newtest=challenge(matrices,region_names,ids,out)
    m.PANELS=tuple(matrices)
    names=list(matrices)
    pd.DataFrame(regions).to_csv(out/'challenge_region_manifest.tsv',sep='\t',index=False)
    (out/'challenge_variant_ids.json').write_text(json.dumps(ids,indent=2),encoding='utf-8')
    matrix_dir=out/'ld_inputs'; matrix_dir.mkdir(exist_ok=True)
    for name,r in matrices.items(): np.savetxt(matrix_dir/f'{name}.tsv',r,delimiter='\t')
    export_diagnostic(val,'validation',names,out); export_diagnostic(newtest,'challenge',names,out)
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    baseline=make_pipeline(StandardScaler(),LogisticRegression(max_iter=3000,random_state=config.seed))
    baseline.fit(train['baseline'],train['y'])
    predictions={}; old_predictions={}; thresholds={}; training=[]
    for mode in ('signed_graph','pooled_mlp','unsigned_graph'):
        pos,neg=operations(matrices,mode); vps=[]; ops=[]; nps=[]
        for seed in config.training_seeds:
            if mode=='signed_graph':
                model=m.SignedLDGuardNet(5,48,.2)
                state=torch.load(HERE.parent/f'ld_guardnet/ld_guardnet_seed_{seed}.pt',weights_only=True)
                model.load_state_dict(state['state_dict'])
            else:
                model,history=m.train_network(train,val,pos,neg,config,seed)
                torch.save({'state_dict':model.state_dict(),'mode':mode,'seed':seed},out/f'{mode}_{seed}.pt')
                training.append({'mode':mode,'seed':seed,'best_epoch':history['best_epoch'],'best_validation_auc':history['best_validation_auc']})
            vps.append(m.predict_network(model,val,pos,neg,128)); ops.append(m.predict_network(model,oldtest,pos,neg,128))
            nps.append(m.predict_network(model,newtest,pos,neg,128))
        vp=np.mean(vps,0); op=np.mean(ops,0); npred=np.mean(nps,0)
        thresholds[mode]=float(np.quantile(1-vp[val['y']==0,0],.95))
        predictions[mode]=npred; old_predictions[mode]=op
        print(f'{mode} trained/evaluated',flush=True)
    vp=baseline.predict_proba(val['baseline']); thresholds['logistic']=float(np.quantile(1-vp[val['y']==0,0],.95))
    predictions['logistic']=baseline.predict_proba(newtest['baseline']); old_predictions['logistic']=baseline.predict_proba(oldtest['baseline'])
    original=pd.read_csv(HERE.parent/'ld_guardnet/ood_predictions.tsv',sep='\t')
    stored=original[['prob_'+c for c in m.CLASSES]].to_numpy()
    maxdiff=float(np.max(np.abs(stored-old_predictions['signed_graph'])))
    assert maxdiff<2e-6, maxdiff
    ci=bootstrap(original,old_predictions['signed_graph'],old_predictions['logistic'])
    original_rows=[]; challenge_rows=[]
    frame=pd.DataFrame(newtest['metadata']); all_predictions=frame.copy()
    for mode,prob in predictions.items():
        scores=1-prob[:,0]; pred=scores>thresholds[mode]
        original_rows.append(dict(model=mode,**m.probability_metrics(oldtest['y'],old_predictions[mode])))
        all_predictions[mode+'_anomaly_score']=scores
        for fault,idx in frame.groupby('fault').indices.items():
            challenge_rows.append(dict(model=mode,fault=fault,n=len(idx),flagged=int(pred[idx].sum()),
                                       flag_rate=float(pred[idx].mean()),threshold=thresholds[mode]))
        challenge_rows.append(dict(model=mode,fault='ALL_ANOMALY_AUC',n=len(scores),flagged='',
                                   flag_rate=roc_auc_score(newtest['y'],scores),threshold=thresholds[mode]))
    all_predictions.to_csv(out/'challenge_predictions.tsv',sep='\t',index=False)
    pd.DataFrame(original_rows).to_csv(out/'original_regime_ablation_metrics.tsv',sep='\t',index=False)
    pd.DataFrame(challenge_rows).to_csv(out/'challenge_metrics.tsv',sep='\t',index=False)
    pd.DataFrame(training).to_csv(out/'ablation_training.tsv',sep='\t',index=False)
    (out/'bootstrap_intervals.json').write_text(json.dumps(ci,indent=2),encoding='utf-8')
    (out/'validation_thresholds.json').write_text(json.dumps(thresholds,indent=2),encoding='utf-8')
    (out/'audit.json').write_text(json.dumps(dict(frozen_checkpoint_max_abs_difference=maxdiff,
        original_input_sha256=hashes,regions=len(region_names),challenge_graphs=len(frame),
        bootstrap_unit='outer replicate; 100 clusters of 32 graphs',bootstrap_repeats=2000,
        scope='review-triggered exploratory evaluation; same EUR donors, different loci; simulated z',
        python=sys.version,torch=torch.__version__,numpy=np.__version__),indent=2),encoding='utf-8')
    # Test the actual gate, with z-only faults and explicit ID-order mismatch separated.
    sys.path.insert(0,str(root/'scripts/hcsmr_code'))
    from hcsmr.reanalysis.preflight_gate import run_preflight
    gate_rows=[]
    for i,record in frame.iterrows():
        variant_ids=ids[record.region]
        r=matrices[record.matrix]
        rows=[dict(target_variant=v,eqtl_se='.1',gwas_se='.1',maf='.25',an='1000',gwas_n=str(record.sample_size),
                   alignment_status='aligned',vcf_status='exact_ref_alt_match',eqtl_beta=str(newtest['raw'][i,0,j]/10),
                   gwas_beta=str(newtest['raw'][i,1,j]/10)) for j,v in enumerate(variant_ids)]
        provenance=[dict(label='synthetic_review_challenge',url='generated://revision-challenge',version_or_build='GRCh37',
            license_or_access='project-generated simulations; source LD terms unchanged',retrieval_date='2026-09-20',
            sha256=hashlib.sha256(r.tobytes()+newtest['raw'][i].tobytes()).hexdigest(),byte_count=r.nbytes+newtest['raw'][i].nbytes,structural_check='passed')]
        result=run_preflight(rows=rows,ld=r,ld_variant_ids=variant_ids,provenance=provenance,
            eqtl_trait_type='quant',gwas_trait_type='quant',case_fraction=.2,ld_rank_policy='diagnostic',max_condition_number=1e12)
        gate_rows.append(dict(region=record.region,replicate=record.replicate,fault=record.fault,status=result['status'],codes=';'.join(result['status_codes'])))
        if record.fault=='valid':
            result=run_preflight(rows=rows[::-1],ld=r,ld_variant_ids=variant_ids,provenance=provenance,
                eqtl_trait_type='quant',gwas_trait_type='quant',case_fraction=.2,ld_rank_policy='diagnostic',max_condition_number=1e12)
            gate_rows.append(dict(region=record.region,replicate=record.replicate,fault='explicit_ID_order',status=result['status'],codes=';'.join(result['status_codes'])))
    pd.DataFrame(gate_rows).to_csv(out/'actual_gate_challenge.tsv',sep='\t',index=False)
    print(json.dumps({'regions':len(region_names),'graphs':len(frame),'ci':ci}),flush=True)


if __name__=='__main__': main()

"""Audit source precision independently; never overwrite original float32 gate audit."""
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from revision_guardnet_validation import m

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--workspace',type=Path,required=True); args=parser.parse_args()
    out=HERE.parent/'revision_validation'; regions=pd.read_csv(out/'challenge_region_manifest.tsv',sep='\t')
    meta=pd.read_csv(out/'challenge_metadata.tsv',sep='\t'); ids=json.loads((out/'challenge_variant_ids.json').read_text())
    z=[np.loadtxt(out/f'challenge_z{i}.tsv',delimiter='\t') for i in (1,2)]
    matrices={}; numeric=[]
    target=out/'source_float64_ld'; target.mkdir(exist_ok=True)
    for ri,row in regions.iterrows():
        a=np.load(args.workspace/row.source_relative,allow_pickle=False)
        ix=np.linspace(0,len(a['variant_ids'])-1,100,dtype=int)
        r=m.nearest_correlation(a['ld'][np.ix_(ix,ix)])
        perm=np.random.default_rng(9941+ri).permutation(100)
        for key,value in ((row.region,r),(row.region+'_scrambled',r[np.ix_(perm,perm)])):
            matrices[key]=value
            original=np.loadtxt(out/'ld_inputs'/f'{key}.tsv',delimiter='\t')
            path=target/f'{key}.tsv'; np.savetxt(path,value,delimiter='\t')
            numeric.append(dict(matrix=key,max_abs_difference=float(np.max(abs(value-original))),
                source64_min_eigenvalue=float(np.linalg.eigvalsh(value)[0]),
                float32_min_eigenvalue=float(np.linalg.eigvalsh(original)[0]),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    sys.path.insert(0,str(args.workspace/'scripts/hcsmr_code'))
    from hcsmr.reanalysis.preflight_gate import run_preflight
    results=[]
    for i,row in meta.iterrows():
        variants=ids[row.region]; r=matrices[row.matrix]
        rows=[dict(target_variant=v,eqtl_se='.1',gwas_se='.1',maf='.25',an='1000',gwas_n=str(row.sample_size),
                   alignment_status='aligned',vcf_status='exact_ref_alt_match',eqtl_beta=str(z[0][i,j]/10),
                   gwas_beta=str(z[1][i,j]/10)) for j,v in enumerate(variants)]
        prov=[dict(label='synthetic_review_challenge_float64',url='generated://revision-challenge-float64',version_or_build='GRCh37',
            license_or_access='project-generated simulations; source LD terms unchanged',retrieval_date='2026-09-20',
            sha256=hashlib.sha256(r.tobytes()+z[0][i].tobytes()+z[1][i].tobytes()).hexdigest(),
            byte_count=r.nbytes+z[0][i].nbytes+z[1][i].nbytes,structural_check='passed')]
        for fault,active in [(row.fault,rows)]+([('explicit_ID_order',rows[::-1])] if row.fault=='valid' else []):
            for representation,matrix in [('source64',r),('float32',r.astype('float32').astype('float64'))]:
                result=run_preflight(rows=active,ld=matrix,ld_variant_ids=variants,provenance=prov,eqtl_trait_type='quant',gwas_trait_type='cc',
                    case_fraction=.2,ld_rank_policy='diagnostic',max_condition_number=1e12)
                results.append(dict(row_index=i,region=row.region,replicate=row.replicate,fault=fault,representation=representation,status=result['status'],codes=';'.join(result['status_codes'])))
    pd.DataFrame(results).to_csv(out/'source_float64_gate_challenge.tsv',sep='\t',index=False)
    pd.DataFrame(numeric).to_csv(out/'ld_precision_comparison.tsv',sep='\t',index=False)
    print(pd.DataFrame(results).groupby(['fault','status']).size().to_string())

if __name__=='__main__': main()

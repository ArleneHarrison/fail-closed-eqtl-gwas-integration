#!/usr/bin/env python3
"""Prepare one frozen GRCh37 region x four GTEx tissues x CAD/HF units."""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,subprocess,os
from pathlib import Path
import numpy as np
import pysam,certifi
from pyliftover import LiftOver

ROOT=Path('PROJECT_ROOT')
EQ={x:ROOT/f'data/eqtl/{x}.all.tsv.gz' for x in ['QTD000131','QTD000136','QTD000251','QTD000256']}
CAD=ROOT/'data/cad/GCST005195_CAD_UKBIOBANK.gz'; HF=ROOT/'data/hf/primary/FORMAT-METAL_Pheno1_EUR.tsv.gz'
PANEL=ROOT/'data/ld/integrated_call_samples_v3.20130502.ALL.panel'
FWD=ROOT/'data/liftover/hg19ToHg38.over.chain.gz'; REV=ROOT/'data/liftover/hg38ToHg19.over.chain.gz'
URL='https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz'
SRC=json.load(open(ROOT/'audit/data_manifest.json',encoding='utf-8'))['assets']; SRC={x['asset_id']:x for x in SRC}
def hbytes(xs): return hashlib.sha256(('\n'.join(xs)+'\n').encode()).hexdigest()
def hfile(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
def primary_hits(lo,c,p):
 hs=[h for h in lo.convert_coordinate('chr'+str(c).removeprefix('chr'),p-1) if h[0] in {f'chr{i}' for i in range(1,23)}|{'chrX','chrY'}]
 return hs
def map_point(lo,back,c,p):
 hs=primary_hits(lo,c,p)
 if len(hs)!=1 or hs[0][2]!='+': return None
 tc,tp=hs[0][0],int(hs[0][1])+1; bs=primary_hits(back,tc,tp)
 if len(bs)!=1 or bs[0][0]!='chr'+str(c).removeprefix('chr') or int(bs[0][1])+1!=p:return None
 return tc.removeprefix('chr'),tp
def tabix_rows(path,region):
 with gzip.open(path,'rt') as f: head=f.readline().rstrip().split('\t')
 z=subprocess.run(['tabix',str(path),region],check=True,capture_output=True,text=True).stdout.splitlines()
 return [dict(zip(head,x.split('\t'))) for x in z]
def complete(r):
 try:return all(r[k] not in {'','NA','nan'} for k in ['gene_id','chromosome','position','variant','ref','alt','beta','se']) and int(r['an'])>0 and int(r['an'])%2==0
 except:return False
def scan_cad(ch,start,end):
 out=[]
 with gzip.open(CAD,'rt') as f:
  rd=csv.DictReader(f,delimiter=' ',skipinitialspace=True)
  for r in rd:
   if r['chr']==str(ch) and start<=int(r['bp'])<=end:out.append(r)
 return out
def eur_samples():
 with open(PANEL) as f:return [r['sample'] for r in csv.DictReader(f,delimiter='\t') if r['super_pop']=='EUR']
def cache_vcf(ch,start,end,out):
 out.mkdir(parents=True,exist_ok=True); raw=out/'region.vcf'; gz=out/'region.vcf.gz'
 if not gz.exists():
  ca=certifi.where();os.environ['CURL_CA_BUNDLE']=ca;os.environ['SSL_CERT_FILE']=ca;u=URL.format(c=ch);src=pysam.VariantFile(u);dst=pysam.VariantFile(str(gz),'wz',header=src.header)
  for rec in src.fetch(str(ch),start-1,end):dst.write(rec)
  dst.close();src.close();pysam.tabix_index(str(gz),preset='vcf',force=True)
 return gz
def parse_vcf(p,keep):
 V={}; samples=[]
 with gzip.open(p,'rt') as f:
  for line in f:
   if line.startswith('##'):continue
   if line.startswith('#CHROM'):
    samples=line.rstrip().split('\t')[9:]; idx=[i for i,s in enumerate(samples) if s in keep]; kept=[samples[i] for i in idx];continue
   a=line.rstrip().split('\t'); ref=a[3].upper(); alts=a[4].upper().split(',')
   if len(alts)!=1:continue
   alt=alts[0]; ds=[]
   for i in idx:
    gt=a[9+i].split(':',1)[0]
    if '.' in gt:ds.append(np.nan)
    else:
     z=gt.replace('|','/').split('/'); ds.append(float(sum(int(q)==1 for q in z)))
   V[(int(a[1]),frozenset((ref,alt)))]=(ref,alt,np.array(ds),a[2])
 return V,kept
def prov(asset,label):
 x=SRC[asset];return {'label':label,'url':x['stable_source_url'],'version_or_build':x['genome_build'],'license_or_access':'public_aggregate_summary_statistics_original_terms_apply','retrieval_date':'2026-09-01','sha256':x['sha256'],'byte_count':x['byte_size'],'structural_check':x['integrity_status']}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--regions',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
 regs=list(csv.DictReader(open(a.regions),delimiter='\t')); r=regs[0]; rid=r['region_id']; ch=int(r['chromosome']); st=int(r['region_start']); en=int(r['region_end'])
 lo19,lo38=LiftOver(str(FWD)),LiftOver(str(REV)); ms=map_point(lo19,lo38,ch,st);me=map_point(lo19,lo38,ch,en)
 if not ms or not me or ms[0]!=me[0] or me[1]<ms[1]:raise SystemExit('E_WINDOW_LIFTOVER')
 lift={'region_id':rid,'source_build':'GRCh37','source':[ch,st,en],'target_build':'GRCh38','target':[ms[0],ms[1],me[1]],'status':'PASS_EXACT_ROUNDTRIP'};json.dump(lift,open(out/'liftover.json','w'),indent=2)
 ers={t:tabix_rows(p,f'{ms[0]}:{ms[1]}-{me[1]}') for t,p in EQ.items()}
 sets=[]
 for t,rows in ers.items():
  c={}
  for q in rows:
   if complete(q):c[q['gene_id']]=c.get(q['gene_id'],0)+1
  sets.append({g for g,n in c.items() if n>=20})
 genes=sorted(set.intersection(*sets));
 if not genes:raise SystemExit('E_NO_ELIGIBLE_GENE')
 gene=genes[0];json.dump({'gene_id':gene,'rule':'lexicographically_first_shared_gene_min20','candidate_count':len(genes)},open(out/'gene_selection.json','w'),indent=2)
 cad=scan_cad(ch,st,en); hf=tabix_rows(HF,f'{ch}:{st}-{en}'); GW={'CAD':cad,'HF':hf}
 vcf=cache_vcf(ch,st,en,out/'vcf_cache'); V,samples=parse_vcf(vcf,set(eur_samples()));
 if len(samples)!=503:raise SystemExit('E_EUR_SAMPLE_COUNT')
 attempts=[]
 for t,allrows in ers.items():
  eqrows=[q for q in allrows if q['gene_id']==gene and complete(q)]; em={}
  for q in eqrows:
   mp=map_point(lo38,lo19,q['chromosome'],int(q['position']))
   if mp and int(mp[0])==ch:em[(mp[1],frozenset((q['ref'].upper(),q['alt'].upper())))]=q
  for outcome,gr in GW.items():
   gm={}
   for q in gr:
    pos=int(q['bp'] if outcome=='CAD' else q['pos_b37']); a1=(q['a1'] if outcome=='CAD' else q['A1']).upper();a2=(q['a2'] if outcome=='CAD' else q['A2']).upper();gm[(pos,frozenset((a1,a2)))]=q
   keys=sorted(set(em)&set(gm)&set(V),key=lambda z:(z[0],sorted(z[1])))
   uid=f'{rid}__{t}__{outcome}'; u=out/'prepared'/uid;u.mkdir(parents=True,exist_ok=True); rows=[]; G=[]; R=[]; miss=[]; vids=[]
   for k in keys:
    e,g=em[k],gm[k];ref,alt,dos,rs=V[k]; mf=float(np.mean(np.isnan(dos)));miss.append({'target_variant':f'{ch}:{k[0]}:{ref}:{alt}','missing_fraction':mf})
    if mf>0.02:continue
    if np.all(np.isnan(dos)) or np.nanstd(dos)==0:continue
    R.append(dos.copy());x=dos.copy();x[np.isnan(x)]=np.nanmean(x);G.append(x);vid=f'{ch}:{k[0]}:{ref}:{alt}';vids.append(vid)
    eb=float(e['beta'])*(1 if e['alt'].upper()==alt else -1); a1=(g['a1'] if outcome=='CAD' else g['A1']).upper();gb=float(g['beta'] if outcome=='CAD' else g['A1_beta'])*(1 if a1==alt else -1)
    gn=int(float(g['N'] if outcome=='CAD' else g['N_total'])); cases=34541 if outcome=='CAD' else 139533
    af=float(np.nanmean(dos)/2);ldmaf=min(af,1-af);eqmaf=min(float(e['maf']),1-float(e['maf']));gf=float(g['af'] if outcome=='CAD' else g['A1_freq']);gmaf=min(gf,1-gf);cases=int(g['N_case']) if outcome=='HF' else 34541
    rows.append({'target_variant':vid,'position':k[0],'eqtl_beta':eb,'eqtl_se':e['se'],'gwas_beta':gb,'gwas_se':g['se'],'maf':ldmaf,'eqtl_maf':eqmaf,'gwas_maf':gmaf,'eqtl_n':int(e['an'])//2,'gwas_n':gn,'gwas_cases':cases,'alignment_status':'aligned' if e['alt'].upper()==alt and a1==alt else 'flipped','vcf_status':'present','gene_id':gene,'tissue_id':t,'outcome_id':outcome})
   if len(G)>=2:
    X=np.column_stack(G); Raw=np.column_stack(R); ld=np.corrcoef(X,rowvar=False); cc=Raw[~np.isnan(Raw).any(axis=1)];ldcc=np.corrcoef(cc,rowvar=False) if len(cc)>=2 else np.full((len(G),len(G)),np.nan)
   else:ld=np.eye(len(G));ldcc=np.eye(len(G));X=np.column_stack(G) if G else np.empty((503,0));cc=X
   with open(u/'summary.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['target_variant']);w.writeheader();w.writerows(rows)
   np.savez_compressed(u/'ld.npz',ld=ld,variant_ids=np.array(vids));np.savez_compressed(u/'ld_complete_case.npz',ld=ldcc,variant_ids=np.array(vids))
   preq=len(em);prg=len(gm);pld=len(V);json.dump({'n_eqtl_variants_pre':preq,'n_gwas_variants_pre':prg,'n_ld_variants_pre':pld,'n_ordered_overlap':len(vids),'ordered_variant_sha256':hbytes(vids)},open(u/'overlap.json','w'),indent=2)
   pv=[prov(t,t),prov('GCST005195' if outcome=='CAD' else 'HERMES2_EUR_outer',outcome),{'label':'1000G_region_cache','url':URL.format(c=ch),'version_or_build':'GRCh37_phase3_EUR503','license_or_access':'public_reference_panel','retrieval_date':'2026-09-01','sha256':hfile(vcf),'byte_count':vcf.stat().st_size,'structural_check':'regional_derivative_bgzip_tabix_EUR503_full_source_hash_not_claimed'},prov('UCSC_hg38ToHg19_chain','liftover_reverse'),prov('UCSC_hg19ToHg38_chain','liftover_forward')]
   if outcome=='HF':pv.append({'label':'HERMES2_Pheno1_member','url':SRC['HERMES2_EUR_outer']['stable_source_url']+'#FORMAT-METAL_Pheno1_EUR.tsv.gz','version_or_build':'GRCh37_b37_hg19','license_or_access':'public_aggregate_summary_statistics_HERMES_conditions_apply','retrieval_date':'2026-09-01','sha256':hfile(HF),'byte_count':HF.stat().st_size,'structural_check':'gzip_tabix_embedded_md5_PASS'})
   json.dump(pv,open(u/'provenance.json','w'),indent=2);json.dump({'policy':'mean_dosage_after_per_variant_missingness_le_0.02','variants':miss,'primary_samples':len(samples),'complete_case_samples':int(len(cc)),'complete_case_file':'ld_complete_case.npz'},open(u/'ld_missingness_trace.json','w'),indent=2)
   attempts.append({'unit_id':uid,'region_id':rid,'anchor_outcome':r['anchor_outcome'],'tissue_id':t,'outcome_id':outcome,'gene_id':gene,'pre_gate_status':'PENDING','pre_gate_code':'PENDING','n_overlap':len(vids)})
 with open(out/'attempts.tsv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(attempts[0]),delimiter='\t');w.writeheader();w.writerows(attempts)
 json.dump({'status':'PASS','region':r,'gene':gene,'units':len(attempts)},open(out/'run_summary.json','w'),indent=2)
if __name__=='__main__':main()

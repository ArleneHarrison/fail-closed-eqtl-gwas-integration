"""Persistent gzip-stream fallback for a read-only 1000G regional audit.

Used only when a local tabix/htslib executable is unavailable.  It retains the
same per-record decision fields as the tabix audit and records that no index
was used; it does not calculate LD.
"""
from __future__ import annotations
import argparse, csv, gzip
from pathlib import Path
from audit_1000g_eur_reference import classify_record

p=argparse.ArgumentParser()
p.add_argument('--vcf', required=True); p.add_argument('--panel', required=True)
p.add_argument('--region', required=True); p.add_argument('--variant-audit-out', required=True)
p.add_argument('--summary-out', required=True); p.add_argument('--selected-samples-out', required=True)
p.add_argument('--max-missing', type=float, default=.05); a=p.parse_args()
chrom, span=a.region.split(':'); start,end=map(int,span.split('-'))
panel=list(csv.DictReader(open(a.panel, encoding='utf-8'), delimiter='\t'))
eur=[r['sample'] for r in panel if r['super_pop']=='EUR']
if len(eur)!=503: raise ValueError(f'expected 503 EUR samples, found {len(eur)}')
rows=[]; samples=None; selected=None
with gzip.open(a.vcf, 'rt', encoding='utf-8') as h:
    for line in h:
        if line.startswith('#CHROM'):
            samples=line.rstrip('\n').split('\t')[9:]; ix={x:i for i,x in enumerate(samples)}
            absent=[x for x in eur if x not in ix]
            if absent: raise ValueError(f'{len(absent)} EUR panel samples absent from VCF')
            selected=[ix[x] for x in eur]; continue
        if line.startswith('#'): continue
        if samples is None: raise ValueError('VCF has no #CHROM header')
        f=line.rstrip('\n').split('\t')
        if f[0].removeprefix('chr') != chrom.removeprefix('chr'): continue
        pos=int(f[1])
        if pos < start: continue
        if pos > end: break
        rows.append(classify_record(f, selected, a.max_missing))
for target in (a.variant_audit_out,a.summary_out,a.selected_samples_out): Path(target).parent.mkdir(parents=True,exist_ok=True)
fields=['chrom','position','ref','alt','filter','eur_n','eur_missing_n','eur_missing_fraction','is_biallelic','include_for_ld','exclusion_reason']
with open(a.variant_audit_out,'w',newline='',encoding='utf-8') as h:
    w=csv.DictWriter(h,fieldnames=fields); w.writeheader(); w.writerows(rows)
with open(a.selected_samples_out,'w',newline='',encoding='utf-8') as h:
    w=csv.writer(h); w.writerow(['sample','super_pop']); w.writerows((x,'EUR') for x in eur)
s={'region':a.region,'eur_samples':len(eur),'records_returned':len(rows),'biallelic_records':sum(r['is_biallelic'] for r in rows),'multiallelic_records':sum(not r['is_biallelic'] for r in rows),'records_eligible_for_ld':sum(r['include_for_ld'] for r in rows),'records_excluded':sum(not r['include_for_ld'] for r in rows),'max_missing_fraction':a.max_missing,'access_method':'gzip_stream_no_index','index_not_used_reason':'persistent tabix/htslib executable unavailable','missing_genotype_handling':'records above threshold retained in audit and excluded; no imputation','status':'audit_complete_no_ld_calculated'}
with open(a.summary_out,'w',newline='',encoding='utf-8') as h:
    w=csv.DictWriter(h,fieldnames=s); w.writeheader(); w.writerow(s)

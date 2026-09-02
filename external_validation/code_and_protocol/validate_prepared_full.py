#!/usr/bin/env python3
"""Independent read-only integrity audit of the 24 x 4 x 2 prepared registry."""
import argparse, csv, hashlib, json
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np

TISSUES=("QTD000131","QTD000136","QTD000251","QTD000256"); OUTCOMES=("CAD","HF")
REQUIRED={"summary.csv","ld.npz","ld_complete_case.npz","overlap.json","provenance.json","ld_missingness_trace.json","preparation_audit.json"}
def horder(xs):return hashlib.sha256(("\n".join(xs)+"\n").encode()).hexdigest()
def readtsv(p):
 with p.open(encoding="utf-8") as f:return list(csv.DictReader(f,delimiter="\t"))
def main():
 p=argparse.ArgumentParser();p.add_argument("--prepared-root",type=Path,required=True);p.add_argument("--regions",type=Path,required=True);p.add_argument("--metadata-lock",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args();r=a.prepared_root
 summary=json.loads((r/"run_summary.json").read_text()); regions=readtsv(a.regions);genes=readtsv(r/"genes.tsv");attempts=readtsv(r/"attempts.tsv");lock=json.loads(a.metadata_lock.read_text())
 assert summary["status"]=="PASS" and summary["mode"]=="full" and summary["regions"]==24 and summary["genes_selected"]==24 and summary["attempts"]==192 and summary["prepared_units"]==192
 assert len(regions)==len(genes)==24 and len(attempts)==192 and all(x["gene_selection_status"]=="SELECTED" for x in genes)
 expected={f"{q['region_id']}__{t}__{o}" for q in regions for t in TISSUES for o in OUTCOMES}; observed={x["unit_id"] for x in attempts}; dirs={q.name for q in (r/"prepared").iterdir() if q.is_dir()}; assert observed==dirs==expected
 source=readtsv(r/"source_verification.tsv");assert source and all(x["status"].startswith("PASS") for x in source)
 dup=readtsv(r/"duplicate_key_audit.tsv"); quarantined=sum(int(x["discordant_duplicates"]) for x in dup if x["source"]=="1000G");assert quarantined==summary["reference_ambiguous_keys_quarantined"] and summary["discordant_duplicate_keys_remaining_in_analysis_view"]==0
 failures=[]; zero=0; total_rows=0; ambiguous_hits=0; primary_inodes={}; cc_inodes={}; outcome_rows=Counter()
 for i,u in enumerate(sorted(expected),1):
  d=r/"prepared"/u; missing=sorted(REQUIRED-{q.name for q in d.iterdir() if q.is_file()})
  if missing:failures.append(f"{u}:missing:{','.join(missing)}");continue
  with (d/"summary.csv").open(encoding="utf-8") as f:rows=list(csv.DictReader(f))
  ov=json.loads((d/"overlap.json").read_text());trace=json.loads((d/"ld_missingness_trace.json").read_text());prov=json.loads((d/"provenance.json").read_text())
  ids=[x["target_variant"] for x in rows]; total_rows+=len(rows);outcome=u.rsplit("__",1)[1];outcome_rows[outcome]+=len(rows)
  if len(rows)==0:zero+=1;assert ov["preparation_code"] in {"E_NO_ORDERED_OVERLAP","E_NO_ELIGIBLE_GENE"}
  assert len(rows)==int(ov["n_ordered_overlap"]) and ov["ordered_variant_sha256"]==horder(ids) and trace["ordered_variant_sha256"]==horder(ids)
  ambiguous_hits+=int(ov["n_matched_discordant_reference_keys"]);assert len(trace["reference_ambiguous_keys_excluded"])==int(ov["n_matched_discordant_reference_keys"])
  assert all(all(x.get(k) not in (None,"") for k in ("url","version_or_build","license_or_access","retrieval_date","sha256","byte_count","structural_check")) for x in prov)
  locked=float(lock["outcomes"][outcome]["case_fraction"])
  for row in rows:
   assert row["outcome_id"]==outcome and abs(float(row["gwas_case_fraction_locked_study"])-locked)<1e-15
   assert 0<float(row["eqtl_maf"])<=.5 and 0<float(row["gwas_maf"])<=.5 and 0<float(row["ld_maf"])<=.5
   assert int(float(row["gwas_cases"]))<=int(float(row["gwas_n"]))
  for name,seen,allow_nan in (("ld.npz",primary_inodes,False),("ld_complete_case.npz",cc_inodes,True)):
   path=d/name;inode=(path.stat().st_dev,path.stat().st_ino)
   if inode not in seen:
    seen[inode]=(path,ids,allow_nan)
  if i%24==0:print(f"validated metadata {i}/192",flush=True)
 def check_npz(task):
  path,ids,allow_nan=task;z=np.load(path,allow_pickle=False);m=z["ld"];v=[str(x) for x in z["variant_ids"]];assert v==ids and m.shape==(len(ids),len(ids))
  if not allow_nan:assert np.isfinite(m).all() and (len(ids)==0 or np.max(np.abs(m-m.T))<1e-8) and (len(ids)==0 or np.max(np.abs(np.diag(m)-1))<1e-8)
  return str(path)
 tasks=list(primary_inodes.values())+list(cc_inodes.values())
 with ThreadPoolExecutor(max_workers=8) as pool:
  for i,_ in enumerate(pool.map(check_npz,tasks),1):
   if i%48==0:print(f"validated LD archives {i}/{len(tasks)}",flush=True)
 status="GREEN" if not failures else "BLOCK"
 report={"status":status,"regions":len(regions),"genes":len(genes),"attempts":len(attempts),"units":len(dirs),"total_summary_rows":total_rows,"summary_rows_by_outcome":dict(outcome_rows),"zero_overlap_units":zero,"reference_ambiguous_keys_quarantined":quarantined,"ambiguous_key_unit_intersections":ambiguous_hits,"unique_primary_LD_files":len(primary_inodes),"unique_complete_case_LD_files":len(cc_inodes),"failures":failures}
 a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n");print(json.dumps(report,sort_keys=True))
 if failures:raise SystemExit(2)
if __name__=="__main__":main()

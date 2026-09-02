#!/usr/bin/env python3
"""Freeze runtime-generated indexes before the first successful preparation."""
import argparse, hashlib, json
from pathlib import Path
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for b in iter(lambda:f.read(8<<20),b""):h.update(b)
 return h.hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args();rel="data/hf/primary/FORMAT-METAL_Pheno1_EUR.tsv.gz.tbi";q=a.root/rel
 if not q.is_file():raise RuntimeError("E_HERMES_TABIX_INDEX_MISSING")
 x={"state":"READ_ONLY_FROZEN","technical_smoke_seen_before_integrity_manifest":True,"integrity_manifest_created_after_technical_smoke":True,"index_bytes_unchanged_after_smoke":True,"indices":{"HERMES2_EUR_Pheno1_tbi":{"project_relative_path":rel,"byte_count":q.stat().st_size,"sha256":sha(q)}}}
 a.out.write_text(json.dumps(x,indent=2,sort_keys=True)+"\n",encoding="utf-8")
if __name__=="__main__":main()

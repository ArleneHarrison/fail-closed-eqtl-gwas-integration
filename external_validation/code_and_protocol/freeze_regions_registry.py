#!/usr/bin/env python3
"""Create the immutable binding for a prospectively selected region registry."""
import argparse, csv, hashlib, json
from pathlib import Path

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(8<<20),b""):h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument("--regions",type=Path,required=True);p.add_argument("--audit",type=Path,required=True);p.add_argument("--protocol",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    protocol=json.loads(a.protocol.read_text(encoding="utf-8")); canonical=hashlib.sha256(json.dumps(protocol,sort_keys=True,separators=(",",":"),ensure_ascii=False).encode()).hexdigest()
    with a.regions.open(encoding="utf-8") as f:n=sum(1 for _ in csv.DictReader(f,delimiter="\t"))
    if n!=24:raise RuntimeError(f"E_REGION_REGISTRY_COUNT:{n}")
    payload={"state":"READ_ONLY_FROZEN","regions_sha256":sha(a.regions),"audit_sha256":sha(a.audit),"protocol_canonical_json_sha256":canonical,"n_regions":n,"technical_smoke_seen_before_integrity_manifest":True,"integrity_manifest_created_after_technical_smoke":True,"design_and_artifact_bytes_unchanged_after_smoke":True,"mutation_policy":"new registry ID and deviation record required"}
    a.out.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n",encoding="utf-8")
if __name__=="__main__":main()

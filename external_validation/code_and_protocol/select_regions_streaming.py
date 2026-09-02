#!/usr/bin/env python3
"""Streaming, deterministic implementation of the frozen 12+12 locus rule."""
import argparse, csv, gzip, hashlib, json
from collections import defaultdict
from pathlib import Path


def scan(path, outcome):
    hits=[]
    with gzip.open(path,"rt",encoding="utf-8") as f:
        reader=csv.DictReader(f,delimiter=" " if outcome=="CAD" else "\t",skipinitialspace=outcome=="CAD")
        for row in reader:
            try:
                chrom=str(row["chr"]).removeprefix("chr").lstrip("0") or "0"
                pos=int(row["bp"] if outcome=="CAD" else row["pos_b37"])
                p=float(row["pval"])
                vid=row["uniqid"] if outcome=="CAD" else row["#key"]
            except (KeyError,TypeError,ValueError):
                continue
            if 0<p<=5e-8 and chrom in {str(i) for i in range(1,23)}:
                hits.append((p,int(chrom),pos,str(vid)))
    hits.sort()
    chosen=[]; bychr=defaultdict(list)
    for p,ch,pos,vid in hits:
        if any(abs(pos-old)<1_000_000 for old in bychr[ch]): continue
        bychr[ch].append(pos); rank=len(chosen)+1
        chosen.append({"region_id":f"{outcome}_L{rank:02d}","anchor_outcome":outcome,"lead_variant":vid,"chromosome":str(ch),"lead_position":pos,"lead_p_value":p,"region_start":max(1,pos-500_000),"region_end":pos+500_000,"selection_rank":rank,"selection_basis":"GWAS_chromosome_position_and_p_value_only"})
        if rank==12: break
    if len(chosen)!=12: raise RuntimeError(f"E_LEAD_LOCUS_SHORTFALL:{outcome}:{len(chosen)}")
    return chosen


def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(8<<20),b""):h.update(b)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument("--cad",type=Path,required=True);p.add_argument("--hf",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    rows=scan(a.cad,"CAD")+scan(a.hf,"HF");a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
    audit={"n_regions":len(rows),"selection_basis":"GWAS_chromosome_position_and_p_value_only","CAD_sha256":sha(a.cad),"HF_sha256":sha(a.hf),"downstream_results_consulted":False}
    a.out.with_suffix(".audit.json").write_text(json.dumps(audit,indent=2)+"\n")
if __name__=="__main__":main()

#!/usr/bin/env bash
set -u

run_root=SERVER_ACCOUNT_ROOT/hcsmr_cvd_reanalysis_20260823/multiregion_transfer_20260828
data_root=SERVER_ACCOUNT_ROOT/hcsmr_cvd_reanalysis_20260823
log_file="$run_root/logs/benchmark.log"

cd "$run_root" || exit 90
export PYTHONPATH="$run_root/scripts"
python3 -u "$run_root/scripts/run_multiregion_transfer_benchmark.py" \
  --qtd "$data_root/raw/QTD000216.all.tsv.gz" \
  --qtd-index "$data_root/raw/QTD000216.all.tsv.gz.tbi" \
  --cad "$data_root/GCST90132314_buildGRCh37.tsv" \
  --vcf "$data_root/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" \
  --panel "$data_root/raw/ld_reference/integrated_call_samples_v3.20130502.ALL.panel" \
  --chain "$data_root/raw/liftover/hg38ToHg19.over.chain.gz" \
  --provenance-json "$run_root/config/input_provenance.json" \
  --regions-json "$run_root/config/locked_regions.json" \
  --out-dir "$run_root/outputs" \
  --tabix /usr/bin/tabix \
  --max-missing 0.05 \
  >"$log_file" 2>&1
exit_code=$?
printf 'BENCHMARK_EXIT_CODE=%s\n' "$exit_code" >>"$log_file"
exit "$exit_code"

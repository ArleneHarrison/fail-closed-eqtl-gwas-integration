#!/usr/bin/env Rscript
# Review amendment: preserve all SuSiE outcomes, evaluate ABF independently.
suppressPackageStartupMessages({library(coloc); library(susieR)})
args <- commandArgs(TRUE)
stopifnot(length(args) == 3)
code_dir <- args[[1]]; ld_dir <- args[[2]]; out_dir <- args[[3]]
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
for (expr in parse(file.path(code_dir, 'run_multisignal_coloc_robustness.R'))) {
  if (is.call(expr) && identical(expr[[1]], as.name('<-')) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name('function'))) eval(expr)
}
eur <- read_ld('EUR'); rg <- nearest_psd(eur$matrix); p <- nrow(rg)
architectures <- c('one_shared','two_shared','shared_plus_specific','distinct_only')
ns <- c(2000L,5000L,20000L)
abf <- list(); k <- 1L
for (replicate in 1:50) for(n in ns) for(architecture in architectures) {
  set.seed(20260920L + replicate * 1000L + match(n,ns)*100L + match(architecture,architectures))
  truth <- causal_configuration(replicate,architecture,p)
  z1 <- simulate_z(rg,n,truth$t1,rep(.085/sqrt(length(truth$t1)),length(truth$t1)))
  z2 <- simulate_z(rg,n,truth$t2,rep(.075/sqrt(length(truth$t2)),length(truth$t2)))
  v <- regional_pph4_abf(make_dataset(z1,rg,n,eur$variants),make_dataset(z2,rg,n,eur$variants))
  abf[[k]] <- data.frame(replicate=replicate,sample_size=n,architecture=architecture,pp_h4_abf=v)
  k <- k+1L
}
abf <- do.call(rbind,abf)
old_dir <- file.path(dirname(code_dir),'multisignal_coloc')
old <- read.delim(file.path(old_dir,'multisignal_replicates.tsv'))
keys <- function(x) paste(x$replicate,x$sample_size,x$architecture,sep=':')
idx <- match(keys(old),keys(abf)); stopifnot(!anyNA(idx))
fresh <- abf$pp_h4_abf[idx]; seen <- is.finite(old$pp_h4_abf)
max_difference <- max(abs(fresh[seen]-old$pp_h4_abf[seen]))
stopifnot(max_difference < 1e-10)
old$pp_h4_abf <- fresh; old$abf_completed_independently <- TRUE
write.table(old,file.path(out_dir,'multisignal_replicates.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
summary <- read.delim(file.path(old_dir,'multisignal_summary.tsv'))
for(i in seq_len(nrow(summary))) {
  x <- old[old$architecture==summary$architecture[i] & old$sample_size==summary$sample_size[i] & old$ld_policy==summary$ld_policy[i],]
  summary$failed[i] <- sum(!x$completed)
  summary$completed_without_finite_pph4[i] <- sum(x$completed & !is.finite(x$max_pp_h4_susie))
  summary$finite_pph4[i] <- sum(is.finite(x$max_pp_h4_susie))
  summary$observed_high_pph4_per_attempt[i] <- sum(x$max_pp_h4_susie>=.8,na.rm=TRUE)/nrow(x)
  nhigh <- sum(x$pp_h4_abf>=.8); ci <- wilson_interval(nhigh,nrow(x))
  summary$median_pp_h4_abf[i] <- median(x$pp_h4_abf)
  summary$h4_ge_0_8_abf[i] <- nhigh
  summary$h4_ge_0_8_abf_denominator[i] <- nrow(x)
  summary$h4_ge_0_8_abf_rate[i] <- nhigh/nrow(x)
  summary$h4_ge_0_8_abf_ci_low[i] <- ci[1]; summary$h4_ge_0_8_abf_ci_high[i] <- ci[2]
}
write.table(summary,file.path(out_dir,'multisignal_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(abf,file.path(out_dir,'independent_abf_600.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(c(paste('ABF overlap maximum absolute difference:',format(max_difference,digits=17)),
 'Preserved all original SuSiE outcomes, including failures and absent comparable pairs.',
 'ABF recomputed on 600 seeded datasets independently of SuSiE; mapped to 1800 paired conditions.',
 capture.output(sessionInfo())),file.path(out_dir,'REVISION_AUDIT.txt'))
print(c(attempted=nrow(old),failed=sum(!old$completed),finite_pph4=sum(is.finite(old$max_pp_h4_susie)),abf_finite=sum(is.finite(old$pp_h4_abf))))

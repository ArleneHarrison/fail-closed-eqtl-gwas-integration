#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(susieR))
args <- commandArgs(TRUE); stopifnot(length(args)==1)
out <- args[[1]]
stopifnot(as.character(packageVersion('susieR'))=='0.14.2')
paths <- list.files(file.path(out,'ld_inputs'),pattern='[.]tsv$',full.names=TRUE)
matrices <- lapply(paths,function(path){
  r <- as.matrix(read.delim(path,header=FALSE)); attr(r,'eigen') <- eigen(r,symmetric=TRUE); r
})
names(matrices) <- sub('[.]tsv$','',basename(paths))
results <- list()
for (split in c('validation','challenge')) {
  metadata <- read.delim(file.path(out,paste0(split,'_metadata.tsv')))
  z1 <- as.matrix(read.delim(file.path(out,paste0(split,'_z1.tsv')),header=FALSE))
  z2 <- as.matrix(read.delim(file.path(out,paste0(split,'_z2.tsv')),header=FALSE))
  keep <- if(split=='validation') which(tolower(as.character(metadata$is_valid))=='true') else seq_len(nrow(metadata))
  rows <- lapply(keep,function(i){
    score <- tryCatch({
      r <- matrices[[metadata$matrix[i]]]
      a <- suppressWarnings(kriging_rss(z1[i,],r,n=metadata$sample_size[i]))
      b <- suppressWarnings(kriging_rss(z2[i,],r,n=metadata$sample_size[i]))
      list(score=max(abs(c(a$conditional_dist$z_std_diff,b$conditional_dist$z_std_diff))),error='')
    },error=function(e)list(score=NA_real_,error=conditionMessage(e)))
    data.frame(row_index=i,split=split,fault=metadata$fault[i],score=score$score,error=score$error)
  })
  results[[split]] <- do.call(rbind,rows)
  write.table(results[[split]],file.path(out,paste0(split,'_kriging_scores.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
  cat(split,':',nrow(results[[split]]),'evaluated\n')
}
writeLines(c('Binary score = maximum |z_std_diff| over both traits and variants.',
  'Used kriging_rss default null-MLE shrinkage; threshold selected on validation-valid inputs only.',
  capture.output(sessionInfo())),file.path(out,'kriging_environment.txt'))

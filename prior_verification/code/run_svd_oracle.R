args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: Rscript run_svd_oracle.R manifest.tsv output.tsv")
}

manifest <- read.delim(args[[1]], stringsAsFactors = FALSE, check.names = FALSE)
rows <- lapply(seq_len(nrow(manifest)), function(i) {
  matrix_values <- as.matrix(read.delim(
    manifest$matrix_path[[i]], header = FALSE, sep = "\t", check.names = FALSE
  ))
  singular_values <- svd(matrix_values, nu = 0, nv = 0)$d
  tolerance <- manifest$rank_tolerance[[i]]
  data.frame(
    case_id = manifest$case_id[[i]],
    r_svd_rank = sum(singular_values > tolerance),
    r_svd_condition_number_2 = max(singular_values) / min(singular_values),
    stringsAsFactors = FALSE
  )
})
output <- do.call(rbind, rows)
write.table(output, args[[2]], sep = "\t", row.names = FALSE, quote = FALSE)

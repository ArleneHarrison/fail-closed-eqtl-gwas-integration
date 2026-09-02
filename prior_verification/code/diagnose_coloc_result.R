#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(susieR); library(coloc)})
args <- commandArgs(trailingOnly = TRUE)
read_ld <- function(path) {
  x <- read.delim(path, check.names = FALSE)
  m <- as.matrix(x[-1]); storage.mode(m) <- "double"
  m <- (m + t(m)) / 2
  eig <- eigen(m, symmetric = TRUE)
  y <- eig$vectors %*% (pmax(eig$values, 1e-8) * t(eig$vectors))
  d <- sqrt(diag(y)); y <- y / tcrossprod(d); y <- (y + t(y)) / 2; diag(y) <- 1
  y
}
r <- read_ld(args[[1]])
set.seed(20260832)
causal <- 37L
draw <- function(x) {e <- eigen(x, symmetric=TRUE); as.vector(e$vectors %*% (sqrt(pmax(e$values,0))*rnorm(nrow(x))))}
z1 <- sqrt(5000)*0.12*r[,causal]+draw(r)
z2 <- sqrt(5000)*0.10*r[,causal]+draw(r)
f1 <- susie_rss(z1, r, n=5000, L=3, estimate_residual_variance=FALSE)
f2 <- susie_rss(z2, r, n=5000, L=3, estimate_residual_variance=FALSE)
result <- coloc.susie(f1, f2)
str(result, max.level = 2)
print(result)

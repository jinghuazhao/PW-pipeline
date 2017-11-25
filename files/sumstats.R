# 25-11-2017 MRC-Epid JHZ

options(digits = 3, scipen=20)

sumstats <- Sys.getenv("sumstats")
if (!file.exists("sumstats.rda")) {
  library(Rmpfr)
  d <- read.table(sumstats,as.is=TRUE)
  colnames(d) <- c("SNP", "A1", "A2", "AF1", "b", "se", "p", "N", "chr", "pos")
  d <- within(d, {
    z_score <- b/se
    P <- format(2*pnorm(mpfr(abs(z_score),100),lower.tail=FALSE))
  })
  save(d,file="sumstats.rda")
}
#To be too specific might be confusing:
#!/bin/env Rscript --vanilla -q
#args <- commandArgs(trailingOnly=TRUE)
#sumstats <- args[1]


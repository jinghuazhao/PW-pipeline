# 28-11-2017 MRC-Epid JHZ

options(digits = 3, scipen=20)

sumstats <- Sys.getenv("sumstats")
mp <- Sys.getenv("mp")
if (!file.exists("sumstats.rda")) {
  d <- read.table(sumstats,as.is=TRUE)
  colnames(d) <- c("SNP", "A1", "A2", "AF1", "b", "se", "p", "N", "chr", "pos")
  d <- within(d, {
    z_score <- b/se
    if(mp=="1")
    {
      library(Rmpfr)
      P <- format(2*pnorm(mpfr(abs(z_score),100),lower.tail=FALSE))
      logP <- as.numeric(-log10(mpfr(P,100)))
    } else {
      P <- p
      logP <- log10(P)
    }
  })
  save(d,file="sumstats.rda")
}
#To be too specific might be confusing:
#!/bin/env Rscript --vanilla -q
#args <- commandArgs(trailingOnly=TRUE)
#sumstats <- args[1]

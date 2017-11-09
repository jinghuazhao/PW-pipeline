# 12-11-2017 MRC-Epid JHZ

options(digits = 3, scipen=20)

if(file.exists("sumstats.rda")) load("sumstats.rda")
else {
  library(Rmpfr)
  sumstats <- Sys.getnev("sumstats")
  d <- read.table(sumstats,as.is=TRUE)
  colnames(d) <- c("SNP", "A1", "A2", "AF1", "b", "se", "P", "N", "chr", "pos")
  d <- within(d, {
    z_score <- b/se
    p <- format(2*pnorm(mpfr(abs(z_score),100),lower.tail=FALSE))
  })
  save("sumstats.rda")
}

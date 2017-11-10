# 4-7-2017 MRC-Epid JHZ

options(digits = 3, scipen=20)
library(Rmpfr)
d <- read.table("doc/mm1kg.dat",as.is=TRUE)
colnames(d) <- c("chr", "pos", "Freq1", "Effect", "StdErr", "p", "NOBS", "SNP")
d <- within(d, {
  Marker <- sprintf("%s:%d",chr,pos)
  z_score <- Effect/StdErr
  P <- format(2*pnorm(mpfr(abs(z_score),100),lower.tail=FALSE))
})[c("SNP","P","Marker")]
snp.is.dot <- with(d,SNP==".")
d[snp.is.dot,'SNP'] <- d[snp.is.dot,'Marker']
names(d) <- c("SNP","P","Marker")
write.table(d[c("SNP","P")],file='vegas2v2',quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")

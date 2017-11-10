# 4-7-2017 MRC-Epid JHZ

options(digits = 3, scipen=20)
library(Rmpfr)
d <- read.table("doc/mm1kg.dat",as.is=TRUE)
colnames(d) <- c("Chr", "Pos", "Freq1", "Effect", "StdErr", "p", "TOTALSAMPLESIZE", "SNP")
d <- within(d, {
  Marker <- sprintf("%s:%d",Chr,Pos)
  z_score <- Effect/StdErr
  P <- format(2*pnorm(mpfr(abs(z_score),100),lower.tail=FALSE))
  logP <- as.numeric(-log10(mpfr(P,100)))
})[c("SNP","Chr","Pos","P","logP","Marker")]
gwas_threshold <- as.numeric(-log10(mpfr(5e-8,100)))
z <- gzfile("depict.txt.gz")
write.table(subset(d,logP>=gwas_threshold),file=z,quote=FALSE,row.name=FALSE,sep="\t")

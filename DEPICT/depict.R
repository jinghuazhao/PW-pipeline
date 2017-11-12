# 12-11-2017 MRC-Epid JHZ

options(digits = 3, scipen=20)
library(Rmpfr)
f <- Sys.getnev("f")
d <- read.table(f,as.is=TRUE)
colnames(d) <- c("SNP", "A1", "A2", "AF1", "b", "se", "p", "N", "Chr", "Pos")
d <- within(d, {
  z_score <- b/se
  P <- format(2*pnorm(mpfr(abs(z_score),100),lower.tail=FALSE))
  logP <- as.numeric(-log10(mpfr(P,100)))
  Marker <- sprintf("%s:%d",Chr,Pos)
})[c("SNP","Chr","Pos","P","logP","Marker")]
gwas_threshold <- as.numeric(-log10(mpfr(5e-8,100)))
z <- gzfile("DEPICT/depict.txt.gz")
write.table(subset(d,logP>=gwas_threshold),file=z,quote=FALSE,row.name=FALSE,sep="\t")

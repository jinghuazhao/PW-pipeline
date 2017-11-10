# 4-7-2017 MRC-Epid JHZ

library(Rmpfr)
d <- read.table("doc/mm1kg.dat",as.is=TRUE)
colnames(d) <- c("CHR", "BP", "Freq1", "Effect", "StdErr", "p", "NOBS", "SNP")
d <- within(d, {
  z_score <- Effect/StdErr
  P <- format(2*pnorm(mpfr(abs(z_score),100),lower.tail=FALSE))
})[c("SNP","CHR","BP","P","NOBS")]
d <- subset(d, SNP!=".")
write.table(d[c("SNP","CHR","BP")],file="MAGMA/magma.snploc",quote=FALSE,row.name=FALSE,col.names=FALSE,sep="\t")
write.table(d[c("SNP","P","NOBS")],file="MAGMA/magma.pval",quote=FALSE,row.name=FALSE,sep="\t")

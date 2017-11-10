# 5-7-2017 MRC-Epid JHZ

library(Rmpfr)
d <- read.table("doc/mm1kg.dat",as.is=TRUE)
colnames(d) <- c("chr", "pos", "Freq1", "Effect", "StdErr", "P", "TOTALSAMPLESIZE", "rsid")
d <- within(d, {
  z_score <- Effect/StdErr
  p <- format(2*pnorm(mpfr(abs(z_score),100),lower.tail=FALSE))
})
write.table(d[,c("chr","pos","z_score","p")],file="magenta",
            quote=FALSE,col.names=FALSE,row.name=FALSE,sep="\t")

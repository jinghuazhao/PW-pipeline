# 12-11-2017 MRC-Epid JHZ

library(Rmpfr)
f <- Sys.getnev("f")
d <- read.table(f,as.is=TRUE)
colnames(d) <- c("SNP", "A1", "A2", "AF1", "b", "se", "P", "N", "chr", "pos")
d <- within(d, {
  z_score <- b/se
  p <- format(2*pnorm(mpfr(abs(z_score),100),lower.tail=FALSE))
})
write.table(d[,c("chr","pos","z_score","p")],file="MAGENTA/magenta",
            quote=FALSE,col.names=FALSE,row.name=FALSE,sep="\t")

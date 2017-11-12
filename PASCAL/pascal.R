# 12-11-2017 MRC-Epid JHZ

sumstats <- Sys.getnev("sumstats")
source(sumstats)
colnames(d) <- c("SNP", "A1", "A2", "AF1", "b", "se", "p", "N", "chr", "pos")
d <- within(d, {
  z_score <- b/se
  P <- format(2*pnorm(mpfr(abs(z_score),100),lower.tail=FALSE))
  Marker <- sprintf("%s:%d",chr,pos)
})[c("SNP","P","Marker")]
snp.is.dot <- with(d,SNP==".")
d[snp.is.dot,'SNP'] <- d[snp.is.dot,'Marker']
names(d) <- c("SNP","P","Marker")
write.table(d[c("SNP","P")],file='PASCAL/vegas2v2',quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")

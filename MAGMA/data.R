# 16-11-2017 MRC-Epid JHZ

sumstats <- Sys.getenv("sumstats_rda")
load(sumstats)
d <- within(d, {
   CHR <- chr
   BP <- pos
   NOBS <- as.integer(N)
})
d <- subset(d, SNP!=".")
write.table(d[c("SNP","CHR","BP")],file="magma.snploc",quote=FALSE,row.name=FALSE,col.names=FALSE,sep="\t")
write.table(d[c("SNP","P","NOBS")],file="magma.pval",quote=FALSE,row.name=FALSE,sep="\t")

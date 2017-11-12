# 12-11-2017 MRC-Epid JHZ

sumstats <- Sys.getnev("sumstats")
source(sumstats)
colnames(d) <- c("SNP", "A1", "A2", "AF1", "b", "se", "P", "NOBS", "CHR", "BP")
d <- d[c("SNP","CHR","BP","P","NOBS")]
d <- subset(d, SNP!=".")
write.table(d[c("SNP","CHR","BP")],file="MAGMA/magma.snploc",quote=FALSE,row.name=FALSE,col.names=FALSE,sep="\t")
write.table(d[c("SNP","P","NOBS")],file="MAGMA/magma.pval",quote=FALSE,row.name=FALSE,sep="\t")

# 15-11-2017 MRC-Epid JHZ

sumstats <- Sys.getenv("sumstats_rda")
load(sumstats)
colnames(d) <- c("SNP", "A1", "A2", "AF1", "b", "se", "P", "NOBS", "CHR", "BP")
d <- d[c("SNP","CHR","BP","P","NOBS")]
d <- within(d, {NOBS <- as.integer(NOBS)})
d <- subset(d, SNP!=".")
write.table(d[c("SNP","CHR","BP")],file="magma.snploc",quote=FALSE,row.name=FALSE,col.names=FALSE,sep="\t")
write.table(d[c("SNP","P","NOBS")],file="magma.pval",quote=FALSE,row.name=FALSE,sep="\t")

# 24-11-2017 MRC-Epid JHZ

sumstats <- Sys.getenv("sumstats_rda")
load(sumstats)
d <- within(d, {
   CHR <- chr
   BP <- pos
   NOBS <- as.integer(N)
})
d <- subset(d, SNP!=".")

db <- Sys.getenv("_db")
write.table(d[c("SNP","CHR","BP")],file=paste0(db,".snploc"),quote=FALSE,row.name=FALSE,col.names=FALSE,sep="\t")
write.table(d[c("SNP","P","NOBS")],file=paste0(db,".pval"),quote=FALSE,row.name=FALSE,sep="\t")

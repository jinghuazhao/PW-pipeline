# 16-11-2017 MRC-Epid JHZ

sumstats <- Sys.getenv("sumstats_rda")
load(sumstats)
d <- within(d, {Marker <- sprintf("%s:%d",chr,pos)})[c("SNP","P","Marker")]
snp.is.dot <- with(d,SNP==".")
d[snp.is.dot,'SNP'] <- d[snp.is.dot,'Marker']
write.table(d[c("SNP","P")],file='vegas2v2',quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")

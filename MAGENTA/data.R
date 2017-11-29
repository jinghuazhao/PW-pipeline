# 16-11-2017 MRC-Epid JHZ

sumstats <- Sys.getenv("sumstats_rda")
load(sumstats)
d <- d[c("chr","pos","z_score","P")]
write.table(d,file="magenta",quote=FALSE,col.names=FALSE,row.name=FALSE,sep="\t")

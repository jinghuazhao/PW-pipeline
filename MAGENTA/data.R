# 15-11-2017 MRC-Epid JHZ

sumstats <- Sys.getenv("sumstats")
source(sumstats)
write.table(d[,c("chr","pos","z_score","p")],file="magenta",
            quote=FALSE,col.names=FALSE,row.name=FALSE,sep="\t")

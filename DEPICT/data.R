# 26-11-2017 MRC-Epid JHZ

sumstats <- Sys.getenv("sumstats_rda")
load(sumstats)
d <- within(d, {
  Chr <- chr
  Pos <- pos
  Marker <- sprintf("%s:%d",Chr,Pos)
})[c("SNP","Chr","Pos","P","logP","Marker")]
gwas_threshold <- as.numeric(-log10(mpfr(5e-8,100)))
z <- gzfile("depict.txt.gz")
write.table(subset(d,logP>=gwas_threshold),file=z,quote=FALSE,row.name=FALSE,sep="\t")

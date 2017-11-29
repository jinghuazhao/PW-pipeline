# 25-11-2017 MRC-Epid JHZ

sumstats <- Sys.getenv("sumstats_rda")
load(sumstats)
mp <- Sys.getenv("mp")
p_threshold <- as.numeric(Sys.getenv("p_threshold"))
d <- within(d, {
  Chr <- chr
  Pos <- pos
  Marker <- sprintf("%s:%d",Chr,Pos)
})[c("SNP","Chr","Pos","P","logP","Marker")]
if(mp=="1") {
  library(Rmpfr)
  gwas_threshold <- as.numeric(-log10(mpfr(p_threshold,100)))
} else gwas_threshold <- -log10(p_threshold)
z <- gzfile("depict.txt.gz")
write.table(subset(d,logP>=gwas_threshold),file=z,quote=FALSE,row.name=FALSE,sep="\t")

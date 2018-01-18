# 18-1-2018 MRC-Epid JHZ

options(digits=3, scipen=20, width=200)
library(openxlsx)
db <- Sys.getenv("db")
xlsx <- paste0(db,".xlsx")
unlink(xlsx, recursive = FALSE, force = FALSE)
wb <- createWorkbook(xlsx)
library(plyr)
load(paste0("MAGENTA/",db,".rda"))
load(paste0("MAGMA/",db,".rda"))
load(paste0("PASCAL/",db,".rda"))
for (tbl in c("_genesetenrichment.txt", "_geneprioritization.txt", "_loci.txt", "_tissueenrichment.txt",".clumped", "_depict.tab"))
{
  file <- paste0(db,tbl)
  sep <- ifelse(tbl==".clumped","", "\t")
  assign(file,read.table(paste0("DEPICT/",file),as.is=TRUE,header=TRUE,sep=sep,quote=""))
}
set <- rename(set,c("FULL_NAME"="GS"))
ps <- rename(ps,c("Name"="GS"))
depict <- rename(get(paste0(db,"_genesetenrichment.txt")),c("Original.gene.set.ID"="GS"))
mmpd <- join_all(list(MAGENTA,set,ps,depict),"GS")
n <- nrow(mmpd)
mmpd <- within(mmpd,
{
   p_MAGENTA <- NOMINAL_GSEA_PVAL_95PERC_CUTOFF
   log10p_MAGENTA <- -log10(p_MAGENTA)
   p_MAGMA <- P
   log10p_MAGMA <- -log10(p_MAGMA)
   p_PASCAL <- chi2Pvalue
   log10p_PASCAL <- -log10(p_PASCAL)
   p_DEPICT <- Nominal.P.value
   log10p_DEPICT <- -log10(p_DEPICT)
   fdr_MAGENTA <- p.adjust(p_MAGENTA, "fdr", n)
   fdr_MAGMA <- p.adjust(p_MAGMA, "fdr", n)
   fdr_PASCAL <- p.adjust(p_PASCAL, "fdr", n)
   fdr_DEPICT <- p.adjust(p_DEPICT, "fdr", n)
   r_MAGENTA <- rank(p_MAGENTA)
   r_MAGMA <- rank(p_MAGMA)
   r_PASCAL <- rank(p_PASCAL)
   r_DEPICT <- rank(p_DEPICT)
})
addWorksheet(wb, db)
writeDataTable(wb, db, mmpd)
save(mmpd,file=paste0(db,".rda"))
kruskal.test(mmpd[c("p_MAGENTA","p_MAGMA","p_PASCAL","p_DEPICT")])
cor(mmpd[c("log10p_MAGENTA","log10p_MAGMA","log10p_PASCAL","log10p_DEPICT")],method="pearson",use="complete.obs")
png(paste0(db,".png"),res=300,height=11,width=8,units="in")
par(mfrow=c(2,3),oma=c(0,0,2,0))
with(mmpd, {
# MAGENTA-MAGMA
  plot(log10p_MAGENTA,log10p_MAGMA,xlab="MAGENTA",ylab="MAGMA",main="MAGENTA-MAGMA comparison")
  abline(0,1)
# MAGENTA-PASCAL
  plot(log10p_MAGENTA,log10p_PASCAL, xlab="MAGENTA",ylab="PASCAL",main="MAGENTA-PASCAL comparison")
  abline(0,1)
# MAGENTA-DEPICT
  plot(log10p_MAGENTA,log10p_DEPICT,xlab="MAGENTA",ylab="DEPICT",main="MAGENTA-DEPICT comparison") 
  abline(0,1)
# MAGMA-PASCAL
  plot(log10p_MAGMA,log10p_PASCAL, xlab="MAGMA",ylab="PASCAL",main="MAGMA-PASCAL comparison")
  abline(0,1)
# MAGMA-DEPICT
  plot(log10p_MAGMA,log10p_DEPICT, xlab="MAGMA",ylab="DEPICT",main="MAGMA-DEPICT comparison") 
  abline(0,1)   
# PASCAL-DEPICT
  plot(log10p_PASCAL,log10p_DEPICT, xlab="PASCAL",ylab="DEPICT",main="PASCAL-DEPICT comparison") 
  abline(0,1)   
# title("Comparison of -log10(p) between software",outer=TRUE)
})
dev.off()
addWorksheet(wb, "depict_discretized.plot")
insertImage(wb, "depict_discretized.plot", paste0(db,".png"),width=16,height=10)
ID <- get(paste0(db,"_genesetenrichment.txt"))[c("Original.gene.set.ID","Original.gene.set.description")]
MAGENTA <- subset(MAGENTA,FDR_95PERC_CUTOFF<0.05)
MAGENTA <- merge(ID,MAGENTA,by.x="Original.gene.set.ID",by.y="GS",all.y=TRUE)
ord <- with(MAGENTA,order(FDR_95PERC_CUTOFF))
addWorksheet(wb, "MAGENTA")
writeDataTable(wb, "MAGENTA", MAGENTA[ord,])
set <- subset(set,P<0.05)
set <- merge(ID,set,by.x="Original.gene.set.ID",by.y="GS",all.y=TRUE)
ord <- with(set,order(P))
addWorksheet(wb, "MAGMA")
writeDataTable(wb, "MAGMA", set[ord,])
ps <- subset(ps,chi2Pvalue<0.05)
ps <- merge(ID,ps,by.x="Original.gene.set.ID",by.y="GS",all.y=TRUE)
ord <- with(ps,order(chi2Pvalue))
addWorksheet(wb, "PASCAL")
writeDataTable(wb, "PASCAL", ps[ord,])
depict <- subset(depict,False.discovery.rate<0.05)
ord <- with(depict,order(False.discovery.rate))
addWorksheet(wb, "DEPICT")
writeDataTable(wb, "DEPICT", depict[ord,])
FDR_all <- subset(mmpd,fdr_MAGENTA<0.05&fdr_PASCAL<0.05&fdr_DEPICT<0.05)[c("GS",paste0("fdr_",c("MAGENTA","PASCAL","DEPICT")))]
addWorksheet(wb, "FDR_all")
writeDataTable(wb, "FDR_all", merge(FDR_all,ID,by.x="GS",by.y="Original.gene.set.ID"))

cat("See\nhttps://github.com/perslab/depict/wiki/DEPICT-result-files-format\n for header information\n")
saveWorkbook(wb, file=xlsx, overwrite=TRUE)
unlink(paste0(db,".png"))

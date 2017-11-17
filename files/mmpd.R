# 17-11-2017 MRC-Epid JHZ

options(digits=3, scipen=20, width=200)
library(openxlsx)
xlsx <- "mmpd.xlsx"
unlink(xlsx, recursive = FALSE, force = FALSE)
wb <- createWorkbook(xlsx)

# depict_discretized_cutoff3.2 results
t <- -log10(0.05/14463)
load("MAGENTA/MAGENTA.rda")
load("MAGMA/MAGMA.rda")
load("PASCAL/PASCAL.rda")
for (tbl in c("_genesetenrichment.txt", "_geneprioritization.txt", "_loci.txt", "_tissueenrichment.txt",".clumped", "_depict.tab"))
{
  file <- paste0("depict",tbl)
  sep <- ifelse(tbl==".clumped","", "\t")
  assign(file,read.table(paste0("DEPICT/",file),as.is=TRUE,header=TRUE,sep=sep,quote=""))
}
mm4 <- merge(MAGENTA,set,by.x=c("GS"),by.y=c("FULL_NAME"))
mmp4 <- merge(mm4,ps,by.x=c("GS"),by.y=c("Name"))
mmpd <- merge(mmp4,depict_genesetenrichment.txt,by.x=c("GS"),by.y=c("Original.gene.set.ID"))
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
with(mmpd,{
   cat("Number of pathways reaching Bonferroni threshold (",0.05/n,"):",
      "MAGENTA=", length(p_MAGENTA[p_MAGENTA<0.05/n]), "MAGMA=",length(p_MAGMA[p_MAGMA<0.05/n]),
      "PASCAL=", length(p_PASCAL[p_PASCAL<0.05/n]), "DEPICT=",length(p_DEPICT[p_DEPICT<0.05/n]), "\n")
  cat("FDR < 0.05:", "MAGENTA=", length(fdr_MAGENTA[fdr_MAGENTA<0.05]), "MAGMA=", length(fdr_MAGMA[fdr_MAGMA<0.05]),
     "PASCAL=",length(fdr_PASCAL[fdr_PASCAL<0.05]),"DEPICT=",length(fdr_DEPICT[fdr_DEPICT<0.05]), "\n")
})
addWorksheet(wb, "depict_discretized_cutoff3.2")
writeDataTable(wb, "depict_discretized_cutoff3.2", mmpd)
save(mmpd,file="mmpd.rda")

h <- c("GS","p_MAGENTA","p_MAGMA","p_PASCAL","p_DEPICT","r_MAGENTA","r_MAGMA","r_PASCAL","r_DEPICT")
mmpd <- mmpd[with(mmpd,order(p_MAGENTA)),c(h,setdiff(names(mmpd),h))]
dim(mmpd)
head(mmpd[h],20)
cor(mmpd[c("p_MAGENTA","p_MAGMA","p_PASCAL","p_DEPICT")],method="pearson",use="complete.obs")
kruskal.test(mmpd[c("p_MAGENTA","p_MAGMA","p_PASCAL","p_DEPICT")])
cor(mmpd[c("log10p_MAGENTA","log10p_MAGMA","log10p_PASCAL","log10p_DEPICT")],method="pearson",use="complete.obs")
png("mmpd.png",res=300,height=11,width=8,units="in")
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
insertImage(wb, "depict_discretized.plot", paste0("mmpd.png"),width=16,height=10)

cat("See\nhttps://github.com/perslab/depict/wiki/DEPICT-result-files-format\n for header information\n")
saveWorkbook(wb, file=xlsx, overwrite=TRUE)
unlink("mmpd.png")

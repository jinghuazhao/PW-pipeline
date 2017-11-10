# 10-11-2017 MRC-Epid JHZ

options(digits=3, scipen=20, width=200)
setwd("/genetics/data/gwas/4-7-17")
library(openxlsx)
xlsx <- "mmp.xlsx"
unlink(xlsx, recursive = FALSE, force = FALSE)
wb <- createWorkbook(xlsx)

# MAGENTA.db results
load("MAGENTA.MAGENTA.db/MAGENTA.rda")
load("MAGMA/MAGMA.rda")
load("PASCAL.MAGENTA_10.db/PASCAL.rda")

mm <- merge(MAGENTA,set2,by.x=c("DB","GS"),by.y=c("FULL_NAME","pathway"))
mm <- within(mm,
{
   p_MAGENTA <- NOMINAL_GSEA_PVAL_95PERC_CUTOFF
   p_MAGMA <- P
   log10p_MAGENTA <- -log10(p_MAGENTA)
   log10p_MAGMA <- -log10(p_MAGMA)
})
mmp <- merge(mm,ps_10,by.x=c("DB","GS"),by.y=c("Name","pathway"))
mmp <- within(mmp,
{
   p_PASCAL <- chi2Pvalue
   log10p_PASCAL <- -log10(p_PASCAL)
})
n1 <- nrow(mmp)
mmp <- within(mmp,
{
   rank_MAGENTA <- rank(p_MAGENTA)
   rank_MAGMA <- rank(p_MAGMA)
   rank_PASCAL <- rank(p_PASCAL)
   fdr_MAGENTA <- p.adjust(p_MAGENTA, "fdr", n1)
   fdr_MAGMA <- p.adjust(p_MAGMA, "fdr", n1)
   fdr_PASCAL <- p.adjust(p_PASCAL, "fdr", n1)
})
h <- c("DB","GS","p_MAGENTA","p_MAGMA","p_PASCAL","rank_MAGENTA","rank_MAGMA","rank_PASCAL")
mmp <- mmp[with(mmp,order(p_MAGENTA)),c(h,setdiff(names(mmp),h))]
dim(mmp)
head(mmp[h],20)
cor(mmp[c("p_MAGENTA","p_MAGMA","p_PASCAL")],method="pearson",use="complete.obs")
kruskal.test(mmp[c("p_MAGENTA","p_MAGMA","p_PASCAL")])
cor(mmp[c("log10p_MAGENTA","log10p_MAGMA","log10p_PASCAL")],method="pearson",use="complete.obs")

addWorksheet(wb, "MAGENTA.db")
writeDataTable(wb, "MAGENTA.db", mmp)

t <- -log10(0.05/2529)
png("mmp.png",res=300,height=11,width=8,units="in")
par(mfrow=c(2,2),cex=0.6,pch=20)
# MAGENTA-MAGMA
with(mmp,plot(log10p_MAGENTA,log10p_MAGMA,xlab="MAGENTA",ylab="MAGMA",main="MAGENTA-MAGMA comparison"))
abline(0,1)
# MAGENTA-PASCAL
with(mmp,plot(log10p_MAGENTA,log10p_PASCAL, xlab="MAGENTA",ylab="PASCAL",main="MAGENTA-PASCAL comparison"))
abline(0,1)
# MAGMA-PASCAL
with(mmp,plot(log10p_MAGMA,log10p_PASCAL, xlab="MAGMA",ylab="PASCAL",main="MAGMA-PASCAL comparison"))
abline(0,1)
dev.off()
addWorksheet(wb, "MAGENTA.plot")
insertImage(wb, "MAGENTA.plot", paste0("mmp.png"),width=15,height=18)

# c2.all.v6.0 results
t <- -log10(0.05/4731)
load("MAGENTA.c2.all.v6.0/MAGENTA_c2.rda")
load("PASCAL.c2.all.v6.0/PASCAL.rda")

mm2 <- merge(MAGENTA_c2,c2,by.x=c("DB"),by.y=c("FULL_NAME"))
mm2 <- within(mm2,
{  
   p_MAGENTA <- NOMINAL_GSEA_PVAL_95PERC_CUTOFF
   p_MAGMA <- P
   log10p_MAGENTA <- -log10(p_MAGENTA)
   log10p_MAGMA <- -log10(p_MAGMA)
})
mmp2 <- merge(mm2,ps,by.x=c("DB"),by.y=c("Name"))
mmp2 <- within(mmp2,
{
   p_PASCAL <- chi2Pvalue
   log10p_PASCAL <- -log10(p_PASCAL)
})
n2 <- nrow(mmp2)
mmp2 <- within(mmp2,
{
   rank_MAGENTA <- rank(p_MAGENTA)
   rank_MAGMA <- rank(p_MAGMA)
   rank_PASCAL <- rank(p_PASCAL)
   fdr_MAGENTA <- p.adjust(p_MAGENTA, "fdr", n2)
   fdr_MAGMA <- p.adjust(p_MAGMA, "fdr", n2)
   fdr_PASCAL <- p.adjust(p_PASCAL, "fdr", n2)
})
h <- c("DB","GS","p_MAGENTA","p_MAGMA","p_PASCAL","rank_MAGENTA","rank_MAGMA","rank_PASCAL")
mmp2 <- mmp2[with(mmp2,order(p_MAGENTA)),c(h,setdiff(names(mmp2),h))]
dim(mmp2)
head(mmp2[h],20)
cor(mmp2[c("p_MAGENTA","p_MAGMA","p_PASCAL")],method="pearson",use="complete.obs")
kruskal.test(mmp2[c("p_MAGENTA","p_MAGMA","p_PASCAL")])
cor(mmp2[c("log10p_MAGENTA","log10p_MAGMA","log10p_PASCAL")],method="pearson",use="complete.obs")
  
addWorksheet(wb, "c2.all")
writeDataTable(wb, "c2.all", mmp2)

png("mmp2.png",res=300,height=11,width=8,units="in")
par(mfrow=c(2,2))
# MAGENTA-MAGMA
with(mmp2,plot(log10p_MAGENTA,log10p_MAGMA,xlab="MAGENTA",ylab="MAGMA",main="MAGENTA-MAGMA comparison"))
abline(0,1)
# MAGENTA-PASCAL
with(mmp2,plot(log10p_MAGENTA,log10p_PASCAL, xlab="MAGENTA",ylab="PASCAL",main="MAGENTA-PASCAL comparison"))
abline(0,1)
# MAGMA-PASCAL
with(mmp2,plot(log10p_MAGMA,log10p_PASCAL, xlab="MAGMA",ylab="PASCAL",main="MAGMA-PASCAL comparison"))
abline(0,1)
dev.off()
addWorksheet(wb, "c2.all.plot")
insertImage(wb, "c2.all.plot", paste0("mmp2.png"),width=15,height=18)

# msigdb.v6.0 results
t <- -log10(0.05/17779)
load("MAGENTA.msigdb.v6.0/MAGENTA.rda")
load("PASCAL.msigdb.v6.0/PASCAL.rda")

mm3 <- merge(MAGENTA,set1,by.x=c("DB"),by.y=c("FULL_NAME"))
mm3 <- within(mm3,
{  
   p_MAGENTA <- NOMINAL_GSEA_PVAL_95PERC_CUTOFF
   p_MAGMA <- P
   log10p_MAGENTA <- -log10(p_MAGENTA)
   log10p_MAGMA <- -log10(p_MAGMA)
})
mmp3 <- merge(mm3,ps,by.x=c("DB"),by.y=c("Name"))
mmp3 <- within(mmp3,
{
   p_PASCAL <- chi2Pvalue
   log10p_PASCAL <- -log10(p_PASCAL)
})
n3 <- nrow(mmp3)
mmp3 <- within(mmp3,
{
   rank_MAGENTA <- rank(p_MAGENTA)
   rank_MAGMA <- rank(p_MAGMA)
   rank_PASCAL <- rank(p_PASCAL)
   fdr_MAGENTA <- p.adjust(p_MAGENTA, "fdr", n3)
   fdr_MAGMA <- p.adjust(p_MAGMA, "fdr", n3)
   fdr_PASCAL <- p.adjust(p_PASCAL, "fdr", n3)
})
h <- c("DB","GS","p_MAGENTA","p_MAGMA","p_PASCAL","rank_MAGENTA","rank_MAGMA","rank_PASCAL")
mmp3 <- mmp3[with(mmp3,order(p_MAGENTA)),c(h,setdiff(names(mmp3),h))]
dim(mmp3)
head(mmp3[h],20)
cor(mmp3[c("p_MAGENTA","p_MAGMA","p_PASCAL")],method="pearson",use="complete.obs")
kruskal.test(mmp3[c("p_MAGENTA","p_MAGMA","p_PASCAL")])
cor(mmp3[c("log10p_MAGENTA","log10p_MAGMA","log10p_PASCAL")],method="pearson",use="complete.obs")
  
addWorksheet(wb, "msigdb")
writeDataTable(wb, "msigdb", mmp3)

png("mmp3.png",res=300,height=11,width=8,units="in")
par(mfrow=c(2,2))
# MAGENTA-MAGMA
with(mmp3,plot(log10p_MAGENTA,log10p_MAGMA,xlab="MAGENTA",ylab="MAGMA",main="MAGENTA-MAGMA comparison"))
abline(0,1)
# MAGENTA-PASCAL
with(mmp3,plot(log10p_MAGENTA,log10p_PASCAL, xlab="MAGENTA",ylab="PASCAL",main="MAGENTA-PASCAL comparison"))
abline(0,1)
# MAGMA-PASCAL
with(mmp3,plot(log10p_MAGMA,log10p_PASCAL, xlab="MAGMA",ylab="PASCAL",main="MAGMA-PASCAL comparison"))
abline(0,1)
dev.off()
addWorksheet(wb, "msigdb.plot")
insertImage(wb, "msigdb.plot", paste0("mmp3.png"),width=15,height=18)

saveWorkbook(wb, file=xlsx, overwrite=TRUE)
unlink("mmp.png")
unlink("mmp2.png")
unlink("mmp3.png")

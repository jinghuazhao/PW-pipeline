# 24-11-2017 MRC-Epid JHZ

options(digits=3, scipen=20, width=200)
library(openxlsx)
db <- Sys.getenv("db")
xlsx <- paste0(db,".xlsx")
unlink(xlsx, recursive = FALSE, force = FALSE)
wb <- createWorkbook(xlsx)

# MAGENTA.db results
load(paste0("MAGENTA/",db,".rda"))
load(paste0("MAGMA/",db,".rda"))
load(paste0("PASCAL/",db,".rda"))

mm <- merge(MAGENTA,set,by.x=c("DB","GS"),by.y=c("FULL_NAME","pathway"))
mm <- within(mm,
{
   p_MAGENTA <- NOMINAL_GSEA_PVAL_95PERC_CUTOFF
   p_MAGMA <- P
   log10p_MAGENTA <- -log10(p_MAGENTA)
   log10p_MAGMA <- -log10(p_MAGMA)
})
mmp <- merge(mm,ps,by.x=c("DB","GS"),by.y=c("Name","pathway"))
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

png(paste0(db,".png"),res=300,height=11,width=8,units="in")
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

saveWorkbook(wb, file=xlsx, overwrite=TRUE)
unlink(paste0(db,".png"))

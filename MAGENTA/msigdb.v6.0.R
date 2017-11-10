# 20-7-2017 MRC-Epid JHZ

library(openxlsx)
options(digits=3, scipen=20)
setwd("/genetics/data/gwas/4-7-17")
xlsx <- "msigdb.v6.0.xlsx"
unlink(xlsx, recursive = FALSE, force = FALSE)
wb <- createWorkbook(xlsx)

# MAGENTA -- gold standard?
load("MAGENTA.msigdb.v6.0/MAGENTA.rda")
m <- MAGENTA
addWorksheet(wb, "MAGENTA")
writeDataTable(wb, "MAGENTA", MAGENTA)

# MAGMA
load("MAGMA/MAGMA.rda")
addWorksheet(wb, "MAGMA")
writeDataTable(wb, "MAGMA", set1)

# PASCAL
load("PASCAL.msigdb.v6.0/PASCAL.rda")
addWorksheet(wb, "PASCAL_genescores")
writeDataTable(wb, "PASCAL_genescores", gs)
addWorksheet(wb, "PASCAL_fusion_genescores")
writeDataTable(wb, "PASCAL_fusion_genescores", fg)
addWorksheet(wb, "PASCAL_PathwaySet")
writeDataTable(wb, "PASCAL_PathwaySet", ps)

for (tbl in c("_genesetenrichment.txt", "_geneprioritization.txt", "_loci.txt", "_tissueenrichment.txt",".clumped", "_depict.tab"))
{
  file <- paste0("depict",tbl)
  sep <- ifelse(tbl==".clumped","", "\t")
  assign(file,read.table(paste0("DEPICT/",file),as.is=TRUE,header=TRUE,sep=sep,quote=""))
  addWorksheet(wb, paste0("DEPICT",tbl))
  dat <- get(file)
  if(tbl=="_genesetenrichment.txt") dat <- within(dat,dat[order(Nominal.P.value),])
  writeDataTable(wb,paste0("DEPICT",tbl),dat)
}
for(s in c("cells","multiplot","system", "tissues"))
{
  i <- paste0("DEPICT_",s)
  addWorksheet(wb, i)
  insertImage(wb, i, paste0("DEPICT/tissue_plot_depict_genenetwork_",s,".png"), width=12, height=6)
}
cat("See\nhttps://github.com/perslab/depict/wiki/DEPICT-result-files-format\n for header information\n")
saveWorkbook(wb, file=xlsx, overwrite=TRUE)

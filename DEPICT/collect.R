# 23-11-2017 MRC-Epid JHZ

options(digits=3, scipen=20, width=200)
library(openxlsx)
xlsx <- "depict.xlsx"
unlink(xlsx, recursive = FALSE, force = FALSE)
wb <- createWorkbook(xlsx)
db <- Sys.getenv("db")

for (tbl in c("_genesetenrichment.txt", "_geneprioritization.txt", "_loci.txt", "_tissueenrichment.txt",".clumped", "_depict.tab"))
{
  file <- paste0(db,tbl)
  sep <- ifelse(tbl==".clumped","", "\t")
  assign(file,read.table(file,as.is=TRUE,header=TRUE,sep=sep,quote=""))
  addWorksheet(wb, paste0("DEPICT",tbl))
  dat <- get(file)
  if(tbl=="_genesetenrichment.txt") dat <- within(dat,dat[order(Nominal.P.value),])
  writeDataTable(wb,paste0("DEPICT",tbl),dat)
}
for(s in c("cells","multiplot","system", "tissues"))
{
  i <- paste0("DEPICT_",s)
  addWorksheet(wb, i)
  insertImage(wb, i, paste0("tissue_plot_depict_genenetwork_",s,".png"), width=12, height=6)
}
cat("See\nhttps://github.com/perslab/depict/wiki/DEPICT-result-files-format\n for header information\n")
saveWorkbook(wb, file=xlsx, overwrite=TRUE)



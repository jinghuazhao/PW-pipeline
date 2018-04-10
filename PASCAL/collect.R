# 10-4-2018 MRC-Epid JHZ

gsps <- function(f,db="MAGENTA",method="sum")
{
  gs <- read.table(paste0(f,".",method,".genescores.txt"),as.is=TRUE,header=TRUE,sep="\t")
  ord <- with(gs,order(pvalue))
  gs <- gs[ord,]
  ps <- read.table(paste0(f,".PathwaySet--",db,"--",method,".txt"),as.is=TRUE,header=TRUE,sep="\t")
  ord <- with(ps,order(chi2Pvalue))
  ps <- ps[ord,]
  save(gs,ps,file=paste0(db,".rda"))
}

db <- Sys.getenv("db")
gsps("vegas2v2",db=db)

options(digits=3, scipen=20, width=200)
library(openxlsx)
xlsx <- paste0(db,".xlsx")
unlink(xlsx, recursive = FALSE, force = FALSE)
wb <- createWorkbook(xlsx)
addWorksheet(wb, "gs")
writeDataTable(wb, "gs", gs)
addWorksheet(wb, "ps")
writeDataTable(wb, "ps", ps)
saveWorkbook(wb, file=xlsx, overwrite=TRUE)


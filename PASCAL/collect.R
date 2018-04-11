# 11-4-2018 MRC-Epid JHZ

gsps <- function(f,db="MAGENTA",method="sum")
{
  gs <- read.table(paste0(f,".",method,".genescores.txt"),as.is=TRUE,header=TRUE,sep="\t")
  ord <- with(gs,order(pvalue))
  gs <- gs[ord,]
  ps <- read.table(paste0(f,".PathwaySet--",db,"--",method,".txt"),as.is=TRUE,header=TRUE,sep="\t")
  ord <- with(ps,order(chi2Pvalue))
  ps <- ps[ord,]
  save(gs,ps,file=paste0(db,".rda"))
  n <- nrow(ps)
  ps <- within(ps,{fdr <- p.adjust(chi2Pvalue,"fdr",n)})
  write.table(ps[c("Name","chi2Pvalue","fdr")],file=paste0(db,".dat"),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
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

if (db=="depict_discretized_cutoff3.2")
{
   prefix <- "depict"
   for (tbl in c("_cluster_results.txt","_summary.txt","_network_table.txt","_nodeattributes.txt"))
   {
     file <- paste0(prefix,tbl)
     assign(file,read.table(file,as.is=TRUE,header=TRUE,sep="\t",quote=""))
     addWorksheet(wb, paste0("PASCAL",tbl))
     dat <- get(file)
     writeDataTable(wb,paste0("PASCAL",tbl),dat)
   }
   # addWorksheet(wb, "PASCAL_network_diagram")
   # insertImage(wb, "PASCAL_network_diagram", paste0(db,"_network_diagram.png"),width=12,height=6)
}

saveWorkbook(wb, file=xlsx, overwrite=TRUE)


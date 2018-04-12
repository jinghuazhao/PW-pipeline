# 12-4-2018 MRC-Epid JHZ

db <- Sys.getenv("db")
suffix <- Sys.getenv("suffix")

get.results <- function(db,date=suffix)
{
  file <- paste0("MAGENTA_pval_GeneSetEnrichAnalysis_",db,"_110kb_upstr_40kb_downstr",date,".results")
  body <- read.table(paste0(db,date,"/",file),sep="\t",as.is=TRUE,skip=1,quote="")
  print(dim(body))
  body
}
r <- get.results(db)

header <- read.table(paste0(db,".dat"),sep="\t",as.is=TRUE,header=TRUE)
names(r) <- names(header)[1:ncol(r)]
dim(r)
ord <- with(r,order(NOMINAL_GSEA_PVAL_95PERC_CUTOFF))
MAGENTA <- r[ord,]
write.table(MAGENTA,paste0(db,".dat"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
sp <- subset(MAGENTA,FDR_95PERC_CUTOFF<0.05,select="GS",drop=TRUE)
write.table(sp,paste0(db,".colnames"),col.names=FALSE,quote=FALSE,row.names=FALSE)
save(MAGENTA,file=paste0(db,".rda"))

options(digits=3, scipen=20, width=200)
library(openxlsx)
xlsx <- paste0(db,".xlsx")
unlink(xlsx, recursive = FALSE, force = FALSE)
wb <- createWorkbook(xlsx)
addWorksheet(wb, "MAGENTA")
writeDataTable(wb, "MAGENTA", MAGENTA)

if (db=="depict_discretized_cutoff3.2")
{
   for (tbl in c("_APCluster_info","_APCluster_cluster","_APCluster_iid"))
   {
     file <- paste0(db,tbl,".txt")
     assign(file,read.table(file,as.is=TRUE,header=TRUE,sep="\t",quote=""))
     addWorksheet(wb, paste0("MAGENTA",tbl))
     dat <- get(file)
     writeDataTable(wb,paste0("MAGENTA",tbl),dat)
   }
   prefix <- "depict"
   for (tbl in c("_cluster_results","_summary","_network_table","_nodeattributes"))
   {
     file <- paste0(prefix,tbl,".txt")
     assign(file,read.table(file,as.is=TRUE,header=TRUE,sep="\t",quote=""))
     addWorksheet(wb, paste0("MAGENTA",tbl))
     dat <- get(file)
     writeDataTable(wb,paste0("MAGENTA",tbl),dat)
   }
   # addWorksheet(wb, "MAGENTA_network_diagram")
   # insertImage(wb, "MAGENTA_network_diagram", paste0(db,"_network_diagram.png"),width=12,height=6)
}
saveWorkbook(wb, file=xlsx, overwrite=TRUE)

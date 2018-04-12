# 12-4-2018 MRC-Epid JHZ

options(width=250)
sets <- function(db)
{
  sets.out <- read.table(paste0(db,".sets.out"), as.is=TRUE, skip=3, header=TRUE)
  ordered <- with(sets.out, order(P))
  require(gap)
  with(sets.out,{
    summary(P)
    pdf(paste0(db,".sets.pdf"))
    qqunif(P)
    dev.off()
  })
  keep_var <- !(names(sets.out)%in%"SET")
  sets.out[ordered, keep_var]
}

db <- Sys.getenv("db")
magma <- sets(db)
ord <- with(magma,order(P))
set <- magma[ord,]
save(set,file=paste0(db,".rda"))

n <- nrow(set)
set <- within(set,{fdr <- p.adjust(P,"fdr",n)})
write.table(set[c("FULL_NAME","P","fdr")],file=paste0(db,".dat"),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

options(digits=3, scipen=20, width=200)
library(openxlsx)
xlsx <- paste0(db,".xlsx")
unlink(xlsx, recursive = FALSE, force = FALSE)
wb <- createWorkbook(xlsx)
addWorksheet(wb, "MAGMA")
writeDataTable(wb, "MAGMA", set)

if (db=="depict_discretized_cutoff3.2")
{
   for (tbl in c("_APCluster_info","_APCluster_cluster","_APCluster_iid"))
   {
     file <- paste0(prefix,tbl,".txt")
     assign(file,read.table(file,as.is=TRUE,header=TRUE,sep="\t",quote=""))
     addWorksheet(wb, paste0("MAGMA",tbl))
     dat <- get(file)
     writeDataTable(wb,paste0("MAGMA",tbl),dat)
   }
   prefix <- "depict"
   for (tbl in c("_cluster_results.txt","_summary.txt","_network_table.txt","_nodeattributes.txt"))
   {
     file <- paste0(prefix,tbl)
     assign(file,read.table(file,as.is=TRUE,header=TRUE,sep="\t",quote=""))
     addWorksheet(wb, paste0("MAGMA",tbl))
     dat <- get(file)
     writeDataTable(wb,paste0("MAGMA",tbl),dat)
   }
   # addWorksheet(wb, "MAGMA_network_diagram")
   # insertImage(wb, "MAGMA_network_diagram", paste0(db,"_network_diagram.png"),width=12,height=6)
}

saveWorkbook(wb, file=xlsx, overwrite=TRUE)


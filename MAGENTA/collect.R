# 24-11-2017 MRC-Epid JHZ

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

header <- read.table("header.dat",sep="\t",as.is=TRUE,header=TRUE)
names(r) <- names(header)[1:ncol(r)]
dim(r)
ord <- with(r,order(NOMINAL_GSEA_PVAL_95PERC_CUTOFF))
MAGENTA <- r[ord,]
write.table(MAGENTA,paste0(db,".dat"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
save(MAGENTA,file=paste0(db".rda"))

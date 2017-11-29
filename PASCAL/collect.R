# 24-11-2017 MRC-Epid JHZ

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


# 17-11-2017 MRC-Epid JHZ

gsps <- function(f,db="MAGENTA",method="sum")
{
  gs <- read.table(paste0(f,".",method,".genescores.txt"),as.is=TRUE,header=TRUE,sep="\t")
  ord <- with(gs,order(pvalue))
  gs <- gs[ord,]
  fg <- read.table(paste0(f,".",method,".fusion.genescores.txt"),as.is=TRUE,header=TRUE,sep="\t")
  ord <- with(fg,order(pvalue))
  fg <- fg[ord,]
  ps <- read.table(paste0(f,".PathwaySet--",db,"--",method,".txt"),as.is=TRUE,header=TRUE,sep="\t")
  ord <- with(ps,order(chi2Pvalue))
  ps <- ps[ord,]
  save(gs,fg,ps,file="PASCAL.rda")
}
gsps("vegas2v2")


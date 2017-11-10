# 14-7-2017 MRC-Epid JHZ

gsps <- function(f,db="MAGENTA_10",method="sum")
{
  gs <- read.table(paste0(f,".",method,".genescores.txt"),as.is=TRUE,header=TRUE,sep="\t")
  ord <- with(gs,order(pvalue))
  gs_10 <- gs[ord,]
  fg <- read.table(paste0(f,".",method,".fusion.genescores.txt"),as.is=TRUE,header=TRUE,sep="\t")
  ord <- with(fg,order(pvalue))
  fg_10 <- fg[ord,]
  ps <- read.table(paste0(f,".PathwaySet--",db,"--",method,".txt"),as.is=TRUE,header=TRUE,sep="\t")
  id <- read.table("MAGENTA_10.id",as.is=TRUE,col.names=c("FULL_NAME","pathway"),sep="\t",quote="")
  ps <- merge(ps,id,by.x="Name",by.y="FULL_NAME",all=TRUE)
  ps <- within(ps,{Name=substring(Name,7)})
  ord <- with(ps,order(chi2Pvalue))
  ps_10 <- ps[ord,]
  save(gs_10,fg_10,ps_10,file="PASCAL.rda")
}
gsps("vegas2v2")


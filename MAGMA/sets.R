# 3-10-2017 MRC-Epid JHZ

sets <- function(id)
{
  sets.out <- read.table(paste0(id,".sets.out"), as.is=TRUE, skip=3, header=TRUE)
  ordered <- with(sets.out, order(P))
  require(gap)
  with(sets.out,{
    summary(P)
    pdf(paste0(id,".sets.pdf"))
    qqunif(P)
    dev.off()
  })
  keep_var <- !(names(sets.out)%in%"SET")
  sets.out[ordered, keep_var]
}

options(width=250)
magma <- sets("magma")
ord1 <- with(magma,order(P))
set1 <- magma[ord1,]
c2 <- sets("c2")
ordc2 <- with(c2,order(P))
c2 <- c2[ordc2,]
## add FULL_NAME
M <- sets("MAGENTA_10")
id <- read.table("MAGENTA_10.id",as.is=TRUE,col.names=c("FULL_NAME","pathway"),sep="\t",quote="")
# library(readstata13)
# id <- read.dta13("MAGENTA.dta")
dim(id)
MAGENTA <- merge(M,id,by="FULL_NAME",all=TRUE)
MAGENTA <- within(MAGENTA,{FULL_NAME=substring(FULL_NAME,7)})
dim(MAGENTA)
ord2 <- with(MAGENTA,order(P))
set2 <- MAGENTA[ord2,]
depict <- sets("depict")
ord3 <- with(depict,order(P))
set3 <- depict[ord3,]
save(set1,c2,set2,set3,file="MAGMA.rda")

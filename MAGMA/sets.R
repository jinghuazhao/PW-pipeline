# 24-11-2017 MRC-Epid JHZ

options(width=250)
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

db <- Sys.getenv("db")
magma <- sets(db)
ord <- with(magma,order(P))
set <- magma[ord,]
save(set,file=paste0(db,".rda"))

# 24-11-2017 MRC-Epid JHZ

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

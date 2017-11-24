#24-11-2017 MRC-Epid JHZ

db <- Sys.getenv("db")
load(paste0(db,".rda"))

GO <- function(prefix)
{
  n <- nrow(mmpd)
  SOFTWARE <- c("MAGENTA","MAGMA","PASCAL","DEPICT")
  for (software in SOFTWARE)
  {
    v <- paste0(prefix,"_",software)
    t <- ifelse(prefix=="p",0.05/n,0.05)
    s <- subset(mmpd, !is.na(mmpd[v]) & mmpd[v]<t)
    o <- with(s, order(s[v]))
    s <- s[o,]
    z <- s[c("Original.gene.set.description",paste0(prefix,"_",SOFTWARE))]
    write.table(z,file=paste0(software,".txt"),quote=FALSE,row.names=FALSE,sep="\t")
    cat(software, ": n (", prefix,") = ",nrow(z), sep="")
    if (prefix=="fdr" & software %in% c("MAGENTA","DEPICT"))
       {
         v <- ifelse(software=="MAGENTA", "FDR_95PERC_CUTOFF", "False.discovery.rate")
         s <- subset(mmpd,!is.na(mmpd[v]) & mmpd[v]<0.05)
         o <- with(s, order(s[v]))
         s <- s[o,]
         cat(" vs", nrow(s), "[provided FDR < 0.05]",sep=" ")
       }
    cat("\n")
  }
}

GO("fdr")

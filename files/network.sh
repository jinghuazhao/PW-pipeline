# 9-12-2017 MRC-Epid JHZ

export prefix=/genetics/bin/DEPICT/data/reconstituted_genesets
export BP=${prefix}/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.txt.gz
export columns=${prefix}/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary.columns.txt
export rows=${prefix}/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary.rows.txt

# DEPICT result with FDR<0.05 for 1418 pathways
export _db=depict
awk 'NR>1{print $1}' ${_db}_genesetenrichment.txt | \
head -1418 > ${_db}.colnames
zgrep -n -T -x -f ${_db}.colnames ${columns} | \
cut -f1 > ${_db}.colid
fn=$(cat ${_db}.colid)
echo $fn > ${_db}.cat
fields=$(sed 's/ /,/g' ${_db}.cat)
gunzip -c ${BP} | \
cut -f1 > ${_db}.rownames
gunzip -c ${BP} | \
cut -f1 --complement | \
cut -f$fields | \
paste ${_db}.rownames - > ${_db}.network
gunzip -c id_descrip.txt.gz  | \
cut -f2 > descrip
fd=$(cat descrip)
echo $fd | \
cut -d' ' -f$fields --output-delimiter=$'\t' > ${_db}.descrip

## v slow!
gunzip -c ${BP} | awk -vFS="\t" -vOFS="\t" -f xpose.awk | grep -x -f ${_db}.colnames | awk -f xpose.awk > ${_db}.network

R --no-save <<END
  require(apcluster)
  require(cluster)
  require(factoextra)
  require(graph)
  require(igraph)
  require(network)
  require(reshape)
  require(RCytoscape)
  require(NbClust)
  set.seed(31415625)
  db <- Sys.getenv("_db")
  nw <- read.table(paste0(db,".network"),as.is=TRUE,header=TRUE,quote="")
  Raw <- nw[,-1]
  tRaw <- t(Raw)
  corRaw <- cor(Raw)
  colnames(corRaw) <- rownames(corRaw) <- names(Raw)
  r <- melt(corRaw)
  e <- cbind(r[1],r[3],r[2])
  write.table(e,file="network.sif",col.names=FALSE,row.names=FALSE,quote=FALSE)
  m <- (abs(corRaw)>0.7)+0
  diag(m) <- 0
  g <- graph.adjacency(m)
  write_graph(g,"network.el","edgelist")
  z <- gzfile("id_descrip.txt.gz")
  id_descrip <- read.table(z,sep="\t",as.is=TRUE,header=TRUE,quote="")
  descrip <- id_descrip[names(Raw)%in%with(id_descrip,Original.gene.set.ID),"Original.gene.set.description"]
  names(Raw) <- descrip
  apres <- apcluster(corSimMat,tRaw,details=TRUE)
  show(apres)
  pdf("network.pdf")
  plot(apres,tRaw)
  heatmap(apres)
  distance <- as.dist(1-abs(corRaw))
  fviz_dist(distance,gradient=list(low="#00AFBB",mid="white",high="#FC4E07"))
  g <- network(m, directed=FALSE)
  plot(pam(distance,9))
  plot(g)
  cl <- kmeans(distance,5,nstart=20)
  fviz_cluster(cl,data=distance)
  gmat <- new("graphAM", adjMat=m, edgemode='directed')
  glist <- as(gmat, 'graphNEL')
  plot(glist)
  dev.off()
  pdf("nb.pdf")
  gap_stat <- clusGap(tRaw, FUN=kmeans, nstart=25, K.max=10, B=5)
  fviz_gap_stat(gap_stat)
  fviz_nbclust(tRaw, kmeans, method="gap_stat", nboot=5)
  nb <- NbClust(tRaw,distance="euclidean",min.nc=2,max.nc=10,method="kmeans",index="all")
  dev.off()
# uninformative
# require(spatstat)
# plot(im(corRaw[nrow(corRaw):1,]),main="Correlation Matrix Map")
# too slow
# require(pvclust)
# cl.bootstrap <- pvclust(Raw,nboot=1000,method.dist="correlation")
# plot(cl.bootstrap)
# pvrect(cl.bootstrap)
  close(z)
END

# http://research.stowers.org/mcm/efg/R/Visualization/cor-cluster/index.htm
# http://www.sthda.com/english/wiki/print.php?id=239

# 6-12-2017 MRC-Epid JHZ

export prefix=/genetics/bin/DEPICT/data/reconstituted_genesets
export BP=${prefix}/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.txt.gz
export columns=${prefix}/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary.columns.txt
export rows=${prefix}/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary.rows.txt

gunzip -c ${BP} | \
awk -vFS="\t" -vOFS="\t" -f xpose.awk | \
grep -f ${_db}.colnames | \
awk -f xpose.awk > ${_db}.network
zgrep -n -T -f ${_db}.colnames ${columns} | \
cut -f1 > ${_db}.colid

R --no-save <<END
  require(cluster)
  require(network)
  require(factoextra)
  require(NbClust)
  require(graph)
  require(RCytoscape)
  set.seed(31415625)
  db <- Sys.getenv("_db")
  nw <- read.table(paste0(db,".network"),as.is=TRUE,header=TRUE)
  Raw <- nw[,-1]
  corRaw <- cor(Raw)
  colnames(corRaw) <- rownames(corRaw) <- names(Raw)
  distance <- as.dist(1-abs(corRaw))
  m <- (abs(corRaw)>0.7)+0
  diag(m) <- 0
  g <- network(m, directed=FALSE)
  pdf("network.pdf")
  plot(pam(distance,9))
  plot(g)
  fviz_dist(distance,gradient=list(low="#00AFBB",mid="white",high="#FC4E07"))
  cl <- kmeans(distance,9,nstart=20)
  fviz_cluster(cl,data=distance)
  tRaw <- t(Raw)
  gap_stat <- clusGap(tRaw, FUN=kmeans, nstart=25, K.max=10, B=5)
  fviz_gap_stat(gap_stat)
  fviz_nbclust(tRaw, kmeans, method="gap_stat", nboot=5)
  nb <- NbClust(tRaw,distance="euclidean",min.nc=2,max.nc=10,method="kmeans",index="all")
  graphmat <- new("graphAM", adjMat=m, edgemode='directed')
  graphlist <- as(graphmat, 'graphNEL')
  dev.off()
# uninformative
# require(spatstat)
# plot(im(corRaw[nrow(corRaw):1,]),main="Correlation Matrix Map")
# too slow
# require(pvclust)
# cl.bootstrap <- pvclust(Raw,nboot=1000,method.dist="correlation")
# plot(cl.bootstrap)
# pvrect(cl.bootstrap)
END

# http://research.stowers.org/mcm/efg/R/Visualization/cor-cluster/index.htm
# http://www.sthda.com/english/wiki/print.php?id=239

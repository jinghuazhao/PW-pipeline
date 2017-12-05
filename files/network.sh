# 5-12-2017 MRC-Epid JHZ

export prefix=/genetics/bin/DEPICT/data/reconstituted_genesets
export BP=${prefix}/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.txt.gz
export columns=${prefix}/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary.columns.txt
export rows=${prefix}/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary.rows.txt

gunzip -c ${BP} | \
awk -vFS="\t" -vOFS="\t" -f xpose.awk | \
grep -f depict_discretized_cutoff3.2.colnames | \
awk -f xpose.awk > depict_discretized_cutoff3.2.network

R --no-save <<END
  nw <- read.table("depict_discretized_cutoff3.2.network",as.is=TRUE,header=TRUE)
  Raw <- nw[,-1]
  corRaw <- cor(Raw)
  dissimilarity <- 1-abs(cor(corRaw))
  distance <- as.dist(dissimilarity)
  require(factoextra)
  fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
  tRaw <- t(Raw)
  set.seed(31415625)
  fviz_nbclust(tRaw, kmeans, method = "gap_stat")
  gap_stat <- clusGap(tRaw, FUN = kmeans, nstart = 25, K.max = 10, B = 50)
  fviz_gap_stat(gap_stat)
  cl <- kmeans(distance, 4, nstart = 20)
  fviz_cluster(cl,data=dissimilarity)
# missing basics?
  require(cluster)
  plot(pam(distance,4))
# not very informative
# require(spatstat)
# plot(im(corRaw[nrow(corRaw):1,]), main="Correlation Matrix Map")
# too slow p-value?
# require(pvclust)
# cl.bootstrap <- pvclust(Raw, nboot=1000,method.dist="correlation")
# plot(cl.bootstrap)
# pvrect(cl.bootstrap)
  require(NbClust)
nb <- NbClust(tRaw, distance = "euclidean", min.nc = 2, max.nc = 10, method = "complete", index ="all")
END
zgrep -n -T -f depict_discretized_cutoff3.2.colnames ${columns} | cut -f1 > depict_discretized_cutoff3.2.colid

# http://research.stowers.org/mcm/efg/R/Visualization/cor-cluster/index.htm
# http://www.sthda.com/english/wiki/print.php?id=239

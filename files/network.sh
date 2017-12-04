# 4-12-2017 MRC-Epid JHZ

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
  require(spatstat)
  plot(im(corRaw[nrow(corRaw):1,]), main="Correlation Matrix Map")
  dissimilarity <- 1-abs(cor(corRaw))
  distance <- as.dist(dissimilarity)
  require(cluster)
  plot(pam(distance,5))
# require(pvclust)
# cluster.bootstrap <- pvclust(Raw, nboot=1000,method.dist="correlation")
# plot(cluster.bootstrap)
# pvrect(cluster.bootstrap)
  cl <- kmeans(dissimilarity, 5, nstart = 25)
END
zgrep -n -T -f depict_discretized_cutoff3.2.colnames ${columns} | cut -f1 > depict_discretized_cutoff3.2.colid

# http://research.stowers.org/mcm/efg/R/Visualization/cor-cluster/index.htm

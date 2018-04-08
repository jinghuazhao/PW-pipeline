# 8-4-2018 MRC-Epid JHZ

db <- Sys.getenv("db")
software <- Sys.getenv("software")
PW_location <- Sys.getenv("PW_location")
options(width=200,stringsAsFactors=FALSE)
set.seed(31415625)
misc_runs <- FALSE

nw <- read.table(paste0(db,".network"),as.is=TRUE,header=TRUE,quote="")
Raw <- nw[,-1]
rownames(Raw) <- nw[,1]
corRaw <- cor(Raw)
distance <- as.dist(1-abs(corRaw))
colnames(corRaw) <- rownames(corRaw) <- names(Raw)
# gephi csv
cat(";",file=paste0(software,".csv"))
write.table(signif(corRaw,4),file=paste0(software,".csv"),append=TRUE,col.names=TRUE,row.names=TRUE,quote=FALSE,sep=";")
require(reshape)
r <- melt(corRaw)
l <- as.numeric(cut(m$value,breaks=c(-1,0.4,0.7,1),right=FALSE,include.lowest=TRUE))
e <- cbind(r[1],"interact",r[2],r[3],l)
names(e) <- c("Source","interact","Target","Pearson_correlation", "Pearson_correlation_discrete")
# Cytoscape sif
write.table(e,file=paste0(software,".sif"),col.names=FALSE,row.names=FALSE,quote=FALSE)
z <- gzfile(paste0(PW_location,"/files/id_descrip.txt.gz"))
id_descrip <- within(read.table(z,sep="\t",as.is=TRUE,header=TRUE,quote=""), 
{
   Original.gene.set.ID <- gsub(":",".",Original.gene.set.ID)
   Original.gene.set.description <- gsub(" ",".",Original.gene.set.description)
})
namesRaw <- data.frame(Original.gene.set.ID=names(Raw),iid=1:ncol(Raw))
descrip <- merge(namesRaw,id_descrip,by="Original.gene.set.ID")
descrip <- descrip[with(descrip,order(iid)),]
names(Raw) <- with(descrip, Original.gene.set.description)
tRaw <- t(Raw)
pdf(paste0(software,".pdf"))
require(apcluster)
apres <- apcluster(corSimMat,tRaw,details=TRUE)
heatmap(apres,Rowv=FALSE,Colv=FALSE,cexRow=0.25,cexCol=0.25)
aggres <- aggExCluster(x=apres)
plot(aggres, cex=0.3, horiz=TRUE, nodePar=list(pch=NA, lab.cex=0.25))
# visualisation
# for (k in 20:2) plot(aggres, tRaw[,features], k=k, main=paste(k, "clusters"))
cutres <- cutree(aggres,k=13)
# apresK <- apclusterK(corSimMat, tRaw, K=13, details=TRUE)
cluster_info <- function(z, features=1:15, showClusters=TRUE, output=TRUE, tag="APCluster")
{
   if(showClusters) show(z)
 # plot(z,tRaw[,features])
   exemplars <- z@exemplars
   idx <- z@idx
   clusters <- z@clusters
   sizes <- unlist(lapply(clusters,length),use.names=FALSE)
   m <- lapply(clusters,"[")
   d <- data.frame(labels(clusters),exemplars,sizes,I(m))
   names(d) <- c("label","exemplar","size","nodes")
   M <- matrix(NA,nrow=nrow(d),ncol=max(sizes),dimnames=list(NULL,paste0("node",1:max(sizes))))
   M[cbind(rep(sequence(nrow(d)),sizes),sequence(sizes))] <- unlist(d[["nodes"]],use.names=FALSE)
   names(clusters) <- labels(clusters)
   i <- data.frame(idx,stack(clusters))
   names(i) <- c("label","member","cluster")
   if (output)
   {
    # note that nodes is dropped below and _fl_ in the name is undesirable
      require(splitstackshape)
      dM <- listCol_w(d,"nodes")
      names(dM)[4:(3+max(sizes))] <- paste0("node",1:max(sizes))
      write.table(dM,file=paste0(db,"_",tag,"_cluster.txt"),quote=FALSE,row.names=FALSE)
      write.table(i,file=paste0(db,"_",tag,"_info.txt"),quote=FALSE,row.names=FALSE)
   }
   list(info=i,cluster=d,M=M)
}
apres_info <- cluster_info(apres)
cutres_info <- cluster_info(cutres,tag="cutree")
if (misc_runs)
{
   require(cluster)
   plot(pam(distance,5))
   cl <- kmeans(tRaw,5,nstart=20)
   plot(tRaw,col=cl$cluster)
   points(cl$centers,col=1:5,pch=21)
   require(factoextra)
   fviz_cluster(cl,data=distance)
   clid <- data.frame(cluster=cl$cluster,id=rownames(tRaw))
   rownames(clid) <- NULL
   clid[with(clid,order(cluster)),]
   require(igraph)
   m <- abs(corRaw)
   diag(m) <- 0
   g <- graph.adjacency(m)
   plot(g)
   write_graph(g,paste0(software,".el"),"edgelist")
   require(network)
   n <- network(m, directed=FALSE)
   plot(n)
   require(GOstats)
   require(Rgraphviz)
 # gene.set.id (without duplicates) is used for corrGraph and it may also be slow to draw the diagram.
   gData <- new("ExpressionSet", exprs=t(nw[,-1]))
   corrGraph = compCorrGraph(gData, tau=0.7)
   edgemode(corrGraph) <- "undirected"
   plot(corrGraph)
   require(graph)
   gmat <- new("graphAM", adjMat=m, edgemode='undirected')
   glist <- as(gmat, 'graphNEL')
   plot(glist)
   gap_stat <- clusGap(tRaw, FUN=kmeans, nstart=25, K.max=10, B=5)
   fviz_dist(distance,gradient=list(low="#00AFBB",mid="white",high="#FC4E07"))
   fviz_gap_stat(gap_stat)
   fviz_nbclust(tRaw, kmeans, method="gap_stat", nboot=5)
   require(NbClust)
   nb <- NbClust(tRaw,distance="euclidean",min.nc=2,max.nc=10,method="kmeans",index="all")
   require(RCytoscape)
   require(spatstat)
   plot(im(corRaw[nrow(corRaw):1,]),main="Correlation Matrix Map")
   require(pvclust)
   cl.bootstrap <- pvclust(Raw,nboot=1000,method.dist="correlation")
   plot(cl.bootstrap)
   pvrect(cl.bootstrap)
}
dev.off()

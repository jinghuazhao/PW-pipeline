# 5-4-2018 MRC-Epid JHZ

db <- Sys.getenv("db")
software <- Sys.getenv("software")
PW_location <- Sys.getenv("PW_location")
options(width=200,digits=2)
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
write.table(format(corRaw,digits=getOption("digits")),file=paste0(software,".csv"),append=TRUE,col.names=TRUE,row.names=TRUE,quote=FALSE,sep=";")
require(reshape)
r <- format(melt(corRaw),digits=getOption("digits"))
library(splitstackshape)
l <- listCol_w(r, "value")[, lapply(.SD, as.numeric), by = .(X1, X2)]
l <- as.numeric(cut(with(l,value_fl_1),breaks=c(0,0.4,0.7,1),right=FALSE,include.lowest=TRUE))
e <- cbind(r[1],"interact",r[2],r[3],l)
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
show(apres)
heatmap(apres,Rowv=FALSE,Colv=FALSE,cexRow=0.4,cexCol=0.25)
aggres <- aggExCluster(x=apres)
plot(aggres, cex=0.3, horiz=TRUE, nodePar=list(pch=NA, lab.cex=0.4))
cutres <- cutree(aggres,k=13)
show(cutres)
apresK <- apclusterK(corSimMat, tRaw, K=13, verbose=TRUE)
show(apresK)
# features <- 1:15
# plot(cutres,tRaw[,features])
# par(mfrow=c(2,2))
# for (k in 20:2) plot(aggres, tRaw[,features], k=k, main=paste(k, "clusters"))
# plot(apresK, tRaw[,features])
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

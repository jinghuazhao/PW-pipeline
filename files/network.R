# 2-3-2018 MRC-Epid JHZ

db <- Sys.getenv("db")
software <- Sys.getenv("software")
PW_location <- Sys.getenv("PW_location")
options(width=200)
set.seed(31415625)

nw <- read.table(paste0(db,".network"),as.is=TRUE,header=TRUE,quote="")
Raw <- nw[,-1]
corRaw <- cor(Raw)
distance <- as.dist(1-abs(corRaw))
colnames(corRaw) <- rownames(corRaw) <- names(Raw)
require(reshape)
r <- melt(corRaw)
e <- cbind(r[1],r[3],r[2])
write.table(e,file=paste0(software,".sif"),col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(subset(e, value>=0.7),file=paste0(software,"-1.sif"),col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(subset(e, value>=0.4 & value<0.7),file=paste0(software,"-2.sif"),col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(subset(e, value<0.4),file=paste0(software,"-3.sif"),col.names=FALSE,row.names=FALSE,quote=FALSE)
m <- (abs(corRaw)>0.7)+0
diag(m) <- 0
tRaw <- t(Raw)
pdf(paste0(software,".pdf"))
require(apcluster)
apres <- apcluster(corSimMat,tRaw,details=TRUE)
show(apres)
heatmap(apres)
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
g <- graph.adjacency(m)
plot(g)
write_graph(g,paste0(software,".el"),"edgelist")
require(network)
n <- network(m, directed=FALSE)
plot(n)
## change as appropriate
if(0) { 
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

z <- gzfile(paste0(PW_location,"/files/id_descrip.txt.gz"))
id_descrip <- within(read.table(z,sep="\t",as.is=TRUE,header=TRUE,quote=""), Original.gene.set.ID <- gsub(":",".",Original.gene.set.ID))
namesRaw <- data.frame(Original.gene.set.ID=names(Raw),iid=1:ncol(Raw))
descrip <- merge(namesRaw,id_descrip,by="Original.gene.set.ID")
descrip <- descrip[with(descrip,order(iid)),]
names(Raw) <- gsub(" ",".",descrip["Original.gene.set.description"])


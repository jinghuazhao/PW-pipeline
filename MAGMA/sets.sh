# 22-6-2017 MRC-Epid JHZ

f=magma.setgenes
for s in $(seq 3); do
    grep _SET${s}_ $f.out | grep -w SET | awk '{print $NF}' > $f-$s.names
    grep -f $f-$s.names /genetics/src/MSigDB/msigdb_v6.0_GMTs/msigdb.v6.0.entrez.gmt > $f-$s.entrez
    grep -f $f-$s.names /genetics/src/MSigDB/msigdb_v6.0_GMTs/msigdb.v6.0.symbols.gmt > $f-$s.symbol
R --vanilla $f-$s <<END
    options(width=200)
    args <- commandArgs(trailingOnly=FALSE)
    entrez <- scan(paste0(args[3],".entrez"),"")
    symbol <- scan(paste0(args[3],".symbol"),"")
    a <- as.data.frame(cbind(entrez,symbol))
    names(a) <- c("ENTREZID","SYMBOL")
    s <-read.table(args[3],header=T)
    m <- merge(s,a,by.x="GENE",by.y="ENTREZID",all.x=TRUE)
    n <- scan(paste0(args[3],".names"),"")
    write.table(m,file=paste0(args[3],".txt"),quote=FALSE,row.names=FALSE,sep="\t")
    csv <- with(m,data.frame(source=n,target=SYMBOL,interaction="GG",value=-log10(P)))
    write.table(csv,file=paste0(args[3],".csv"),quote=FALSE,col.names=FALSE,row.names=FALSE,sep=",")
END
echo -e "source,target,interaction,value" > $f.csv
cat $f-*.csv >> $f.csv
done
## it doesn't seem to match entrezID
f=MAGENTA.setgenes
for s in $(seq 14); do
    grep _SET${s}_ $f.out | grep -w SET | awk '{print $NF}' > $f-$s.tmp
    grep -f $f-$s.tmp MAGENTA.id | awk -vFS="\t" '{print $NF}' > $f-$s.names
    grep -f $f-$s.names /genetics/src/MSigDB/msigdb_v6.0_GMTs/msigdb.v6.0.entrez.gmt > $f-$s.entrez
    grep -f $f-$s.names /genetics/src/MSigDB/msigdb_v6.0_GMTs/msigdb.v6.0.symbols.gmt > $f-$s.symbol
R --vanilla $f-$s <<END
    options(width=200)
    args <- commandArgs(trailingOnly=FALSE)
    entrez <- scan(paste0(args[3],".entrez"),"")
    symbol <- scan(paste0(args[3],".symbol"),"")
    a <- as.data.frame(cbind(entrez,symbol))
    names(a) <- c("ENTREZID","SYMBOL")
    s <-read.table(args[3],header=T)
    m <- merge(s,a,by.x="GENE",by.y="ENTREZID",all.x=TRUE)
    n <- scan(paste0(args[3],".names"),"")
    write.table(m,file=paste0(args[3],".txt"),quote=FALSE,row.names=FALSE,sep="\t")
    csv <- with(m,data.frame(source=n,target=SYMBOL,interaction="GG",value=-log10(P)))
    write.table(csv,file=paste0(args[3],".csv"),quote=FALSE,col.names=FALSE,row.names=FALSE,sep=",")
END
echo -e "source,target,interaction,value" > $f.csv
cat $f-*.csv >> $f.csv
done
exit

## the following looks appealing but has leakage in gene symbols

f=magma.setgenes
for s in $(seq 3); do
grep _SET${s}_ $f.out > $f-$s
R --vanilla $f-$s <<END
    options(width=100)
    require(EnsDb.Hsapiens.v86)
    args <- commandArgs(trailingOnly=FALSE)
    print(args[3])
    s <-read.table(args[3],header=T)
    k <- with(s,paste0(GENE))
#   more accurate but leading to redundant gene symbols
#   select(EnsDb.Hsapiens.v86, keys = k, columns = c("ENTREZID", "SYMBOL", "GENEID"), keytype = "ENTREZID")
    a <- select(EnsDb.Hsapiens.v86, keys = k, columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID")
    m <- merge(s,a,by.x="GENE",by.y="ENTREZID",all.x=TRUE)
    na <- with(m,is.na(SYMBOL))
    m[na,"SYMBOL"] <- m[na,"GENE"]
    write.table(m,file=paste0(args[3],".txt"),quote=FALSE,row.names=FALSE,sep="\t")
END
done
for s in $(seq 14)
do
  f=MAGENTA.setgenes
  grep _SET${s}_ $f.out > $f-$s
done

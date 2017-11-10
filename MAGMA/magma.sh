# 30-9-2017 MRC-Epid JHZ

cd /genetics/data/gwas/4-7-17/MAGMA

MAGMA=/genetics/bin/MAGMA
MSigDB=/genetics/src/MSigDB/msigdb_v6.0_GMTs/

# Annotation
magma --annotate window=50,50 --snp-loc magma.snploc --gene-loc $MAGMA/NCBI37.3.gene.loc --out magma

# Gene analysis - SNP p-values
magma --bfile $MAGMA/g1000_eur --pval magma.pval ncol=NOBS --gene-annot magma.genes.annot --out magma

# Gene-set analysis, msigdb.v6.0.entrez.gmt
magma --gene-results magma.genes.raw --set-annot $MSigDB/msigdb.v6.0.entrez.gmt self-contained --out magma
# Gene-set analysis, c2.all.v6.0.entrez.gmt
magma --gene-results magma.genes.raw --set-annot $MSigDB/c2.all.v6.0.entrez.gmt self-contained --out c2

# MAGENTA database
awk -vFS="\t" '{$1=sprintf("%00005d-%s",NR,$1);$2=""};1' ../MAGENTA.db | awk '{$2=$2};1'> MAGENTA.db
cut -f1,2 ../MAGENTA.db | awk -vFS="\t" -vOFS="\t" '{$1=sprintf("%00005d-%s",NR,$1);gsub(/\t/," ",$2)};1' > MAGENTA.id
magma --gene-results magma.genes.raw --set-annot MAGENTA.db self-contained --out MAGENTA

# MAGENTA database and those with at least 10 genes
awk -vFS="\t" '{$1=sprintf("%00005d-%s",NR,$1);$2=""};1' ../MAGENTA.db | awk '{$2=$2};1'|awk 'NF>10' > MAGENTA_10.db
awk -vFS="\t" -vOFS="\t" '{$1=sprintf("%00005d-%s",NR,$1);gsub(/\t/," ",$2)};1' ../MAGENTA.db|awk -vFS="\t" -vOFS="\t" 'NF>11{print $1,$2}' > MAGENTA_10.id
magma --gene-results magma.genes.raw --set-annot MAGENTA_10.db self-contained --out MAGENTA_10

# depict.db
magma --gene-results magma.genes.raw --set-annot depict_discretized_cutoff3.2.gmt self-contained --out depict


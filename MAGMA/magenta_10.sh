# MAGENTA database and those with at least 10 genes
awk -vFS="\t" '{$1=sprintf("%00005d-%s",NR,$1);$2=""};1' ../MAGENTA.db | awk '{$2=$2};1'|awk 'NF>10' > MAGENTA_10.db
awk -vFS="\t" -vOFS="\t" '{$1=sprintf("%00005d-%s",NR,$1);gsub(/\t/," ",$2)};1' ../MAGENTA.db|awk -vFS="\t" -vOFS="\t" 'NF>11{print $1,$2}' > MAGENTA_10.id
magma --gene-results magma.genes.raw --set-annot MAGENTA_10.db self-contained --out MAGENTA_10

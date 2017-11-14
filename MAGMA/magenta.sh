# MAGENTA database

awk -vFS="\t" '{$1=sprintf("%00005d-%s",NR,$1);$2=""};1' ../MAGENTA.db | awk '{$2=$2};1'> magenta.db
cut -f1,2 ../magenta.db | awk -vFS="\t" -vOFS="\t" '{$1=sprintf("%00005d-%s",NR,$1);gsub(/\t/," ",$2)};1' > magenta.id
magma --gene-results magma.genes.raw --set-annot magenta.db self-contained --out magenta

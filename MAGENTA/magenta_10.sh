# 12-7-2017

awk -vOFS="\t" '{gsub(/ /, "\t", $0);print $1, $0}' ../MAGMA/MAGENTA_10.db > MAGENTA_10.db

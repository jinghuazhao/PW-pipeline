# 6-7-2017 MRC-Epid JHZ

awk '(NR==1){gsub(/\#/,"",$0);print}' MAGENTA.db_10000perm_Jul05_17/\
MAGENTA_pval_GeneSetEnrichAnalysis_MAGENTA.db_110kb_upstr_40kb_downstr_10000perm_Jul05_17.results > header.dat

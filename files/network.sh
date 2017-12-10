# 10-12-2017 MRC-Epid JHZ

export prefix=/genetics/bin/DEPICT/data/reconstituted_genesets
export BP=${prefix}/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.txt.gz
export columns=${prefix}/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary.columns.txt
export rows=${prefix}/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary.rows.txt

# DEPICT result with FDR<0.05 for 1418 pathways
export _db=depict
export N=1418

## fast cut!
awk 'NR>1{print $1}' ${_db}_genesetenrichment.txt | \
head -$N > ${_db}.colnames
zgrep -n -T -x -f ${_db}.colnames ${columns} | \
cut -f1 > ${_db}.colid
fn=$(cat ${_db}.colid)
echo $fn > ${_db}.cat
fields=$(sed 's/ /,/g' ${_db}.cat)
gunzip -c ${BP} | \
cut -f1 > ${_db}.rownames
gunzip -c ${BP} | \
cut -f1 --complement | \
cut -f$fields | \
paste ${_db}.rownames - > ${_db}.network

## vv slow!
gunzip -c ${BP} | awk -vFS="\t" -vOFS="\t" -f xpose.awk | grep -x -f ${_db}.colnames | awk -f xpose.awk > ${_db}.network

R --no-save -q < network.R > network.log

# http://research.stowers.org/mcm/efg/R/Visualization/cor-cluster/index.htm
# http://www.sthda.com/english/wiki/print.php?id=239

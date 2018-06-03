# 3-6-2018 MRC-Epid JHZ

export prefix=/genetics/bin/DEPICT/data/reconstituted_genesets
export BP=${prefix}/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.txt.gz
export columns=${prefix}/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary.columns.txt
export rows=${prefix}/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary.rows.txt

export software=$1
## fast cut!
case $1 in
magenta)
  echo MAGENTA network analysis
  export N=$(awk '(NR==1||$11<0.05)' ${db}.dat|awk 'END{print NR}')
  awk 'NR>1{print $2}' ${db}.dat | \
  head -$N > ${_db}.colnames
  ;;
magma)
  echo MAGMA network analysis
  export N=$(awk 'NR>7&&$6<0.05' ${db}.sets.out|awk 'END{print NR}')
  awk 'NR>7{print $6}' ${db}.sets.out | \
  head -$N > ${_db}.colnames
  ;;
pascal)
  echo PASCAL network analysis
  export N=$(awk 'NR==1||$2<0.05' vegas2v2.PathwaySet--${db}--sum.txt|awk 'END{print NR}')
  awk 'NR>1{print $1}' vegas2v2.PathwaySet--${db}--sum.txt | \
  head -$N > ${_db}.colnames
  ;;
depict)
  echo depict network analysis
  export N=$(awk -vFS="\t" '($4<0.05||$4=="<0.01"||$4=="<0.05")' ${db}_genesetenrichment.txt|awk 'END{print NR}')
  awk 'NR>1{print $1}' ${db}_genesetenrichment.txt | \
  head -$N > ${_db}.colnames
  ;;
*)
  echo not implemented
esac

zgrep -n -T -x -f ${_db}.colnames ${columns} | \
sed 's/://g' | \
cut -f1 > ${_db}.colid
fn=$(cat ${_db}.colid)
echo $fn > ${_db}.cat
fields=$(sed 's/ /,/g' ${_db}.cat)
gunzip -c ${BP} | \
cut -f1 > ${_db}.rownames
gunzip -c ${BP} | \
cut -f1 --complement | \
cut -f$fields | \
paste ${_db}.rownames - > ${db}.network

## vv slow!
# gunzip -c ${BP} | awk -vFS="\t" -vOFS="\t" -f xpose.awk | grep -x -f ${db}.colnames | awk -f xpose.awk > ${db}.network

R --no-save -q < $PW_location/files/network.R > ${db}_network.log

# http://research.stowers.org/mcm/efg/R/Visualization/cor-cluster/index.htm
# http://www.sthda.com/english/wiki/print.php?id=239

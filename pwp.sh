#!/bin/bash
# 12-11-2017 MRC-Epid JHZ

## settings

export PW_location=/genetics/bin/PW-pipeline
export use_UCSC=0

### software

export magenta=1
export magma=1
export pascal=1
export depict=1

### databases

export c2=/genetics/src/MSigDB/msigdb_v6.0_GMTs/c2.all.v6.0.entrez.gmt
export msigdb=/genetics/src/MSigDB/msigdb_v6.0_GMTs/msigdb.v6.0.entrez.gmt
export depict2=/genetics/bin/PASCAL/resources/genesets/depict_discretized_cutoff3.2.gmt
MAGMA=/genetics/bin/MAGMA
MSigDB=/genetics/src/MSigDB/msigdb_v6.0_GMTs/
export msigdb_c2=0
export msigdb=0
export magenta_db=1
export depict_db=1
export depict_db2=1

## indidivual analyses according to request

export sumstats=${PW_location}/files/sumstats.R

if [ $magenta -eq 1 ]; then
   echo "MAGENTA"
   if [ ! -d "MAGENTA" ]; then
      mkdir MAGENTA
   fi
   R -q --no-save < ${PW_location}/MAGENTA/data.R
   if [ $magenta_db -eq 1 ]; then
      export db=magenta
      qsub ${PW_location}/MAGENTA/magenta.sh
   elif [ $msigdb_c2 -eq 1 ]; then
      export db=c2
      awk '{$2=$1; $1="c2"; print}' $c2 > MAGENTA/c2.db
      qsub ${PW_location}/MAGENTA/c2.sh
   elif [ $msigdb -eq 1 ]; then
      export db=msigdb
      awk '{$2=$1; $1="msigdb"; print}' $msigdb > MAGENTA/msigdb.db
      qsub ${PW_location}/MAGENTA/msigdb.sh
   else
      export db=depict2
      awk '{FS=OFS="\t";$2=$1;$1="depict";print}' $depict2 > MAGENTA/depict.db
      qsub ${PW_location}/MAGENTA/depict2.sh
   fi
   export suffix=MAGENTA.db_10000perm_Jul05_17
   awk '(NR==1){gsub(/\#/,"",$0);print}' $suffix/MAGENTA_pval_GeneSetEnrichAnalysis_${db}.db_110kb_upstr_40kb_downstr_${suffix}.results > MAGENTA/header.dat
   R -q --no-save < ${PW_location}/MAGENTA/collect.R > MAGENTA/collect.log
fi

if [ $magma -eq 1 ]; then
   echo "MAGMA"
   if [ ! -d "MAGMA" ]; then
      mkdir MAGMA
   fi
   R -q --no-save < ${PW_location}/MAGMA/data.R
   # Annotation
   magma --annotate window=50,50 --snp-loc magma.snploc --gene-loc $MAGMA/NCBI37.3.gene.loc --out MAGMA/magma
   # Gene analysis - SNP p-values
   magma --bfile $MAGMA/g1000_eur --pval magma.pval ncol=NOBS --gene-annot magma.genes.annot --out MGAMA/magma
   if [ $magenta_db -eq 1 ]; then
      qsub ${PW_location}/MAGMA/magenta.sh
   elif [ $msigdb_c2 -eq 1 ]; then
      qsub ${PW_location}/MAGMA/c2.sh
   elif [ $msigbdb -eq 1 ]; then
      qsub ${PW_location}/MAGMA/msigdb.sh
   else
      qsub ${PW_location}/MAGMA/depict.sh
   fi
   R -q --no-save < ${PW_location}/MAGMA/collect.R > MAGMA/collect.log
fi

if [ $pascal -eq 1 ]; then
   echo "PASCAL"
   if [ ! -d "PASCAL" ]; then
      mkdir PASCAL
   fi
   R -q --no-save < ${PW_location}/PASCAL/data.R
   if [ $magenta_db -eq 1 ]; then
      qsub ${PW_location}/PASCAL/magenta.sh
   elif [ $msigdb_c2 -eq 1 ]; then
      qsub ${PW_location}/PASCAL/c2.sh
   elif [ $msigbdb -eq 1 ]; then
      qsub ${PW_location}/PASCAL/msigdb.sh
   else
      qsub ${PW_location}/PASCAL/depict.sh
   fi
   R -q --no-save < $PW_location}/PASCAL/collect.R > PASCAL/collect.log
fi

if [ $depict -eq 1 ]; then
   echo "DEPICT"
   if [ ! -d "DEPICT" ]; then
      mkdir DEPICT
   fi
   R -q --no-save < ${PW_location}/DEPICT/data.R
   if [ $depict_db -eq 1 ]; then
      qsub ${PW_location}/DEPICT/depict.sh
   elif [ $depict_db2 -eq 1 ]; then
      qsub ${PW_location}/DEPICT/depict2.sh
   fi
   R -q --no-save < ${PW_location}/DEPICT/collect.R > DEPICT/collect.log
fi

## collection into Excel

R -q --no-save < ${PW_location}/files/mmp.R > mmp.log

if [ $depict_db2 -eq 1]; then
    R -q --no-save < ${PW_location}/files/mmpd.R > mmpd.log
fi

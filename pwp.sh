#!/bin/bash
# 14-11-2017 MRC-Epid JHZ

## SETTINGS

export MAGENTA=/genetics/bin/MAGENTA_software_package_vs2_July2011
export MAGMA=/genetics/bin/MAGMA
export MSigDB=/genetics/src/MSigDB/msigdb_v6.0_GMTs/
export PW_location=/genetics/bin/PW-pipeline
export use_UCSC=0

### software

export magenta=1
export magma=1
export pascal=1
export depict=1

### databases

export c2=$MSigDB/msigdb_v6.0_GMTs/c2.all.v6.0.entrez.gmt
export msigdb=$MSigDB/msigdb_v6.0_GMTs/msigdb.v6.0.entrez.gmt
export depict2=/genetics/bin/PASCAL/resources/genesets/depict_discretized_cutoff3.2.gmt
export msigdb_c2=0
export msigdb=0
export magenta_db=1
export depict_db=1
export depict_db2=1

## indidivual analyses according to request

export sumstats=${PW_location}/files/sumstats.R

if [ $magenta_db -eq 1 ]; then
   if [ ! -d "MAGENTA" ]; then
      mkdir MAGENTA
   fi
   cd MAGENTA
   for f in $(ls $MAGENTA/*_db); do ln -sf $f; done
   cat GO_terms_BioProc_MolFunc_db Ingenuity_pathways_db KEGG_pathways_db PANTHER_BioProc_db PANTHER_MolFunc_db PANTHER_pathways_db | \
   awk '{$1="MAGENTA"};1' FS="\t" OFS="\t" > magenta.db
   cd -
fi
if [ $magenta -eq 1 ]; then
   echo "MAGENTA"
   if [ ! -d "MAGENTA" ]; then
      mkdir MAGENTA
   fi
   cd MAGENTA
   R -q --no-save < ${PW_location}/MAGENTA/data.R
   for f in $(ls $MAGENTA); do ln -sf $f; done
   for f in ($ls $PW_location/MAGENTA); do ln -sf $f; done
   if [ $magenta_db -eq 1 ]; then
      export db=magenta
      qsub -V -sync y ${PW_location}/MAGENTA/magenta.sh
   elif [ $msigdb_c2 -eq 1 ]; then
      export db=c2
      awk '{$2=$1; $1="c2"; print}' $c2 > c2.db
      awk '{$2=$1; $1="c2"};1' FS="\t" OFS="\t" $c2 > c2.db
      qsub -V -sync y ${PW_location}/MAGENTA/c2.sh
   elif [ $msigdb -eq 1 ]; then
      export db=msigdb
      awk '{$2=$1; $1="msigdb"};1' FS="\t" OFS="\t" $msigdb > msigdb.db
      qsub -V -sync y ${PW_location}/MAGENTA/msigdb.sh
   else
      export db=depict2
      awk '{FS=OFS="\t";$2=$1;$1="depict";print}' $depict2 > MAGENTA/depict.db
      qsub -V -sync y ${PW_location}/MAGENTA/depict2.sh
   fi
   export suffix=MAGENTA.db_10000perm_Jul05_17
   awk '(NR==1){gsub(/\#/,"",$0);print}' $suffix/MAGENTA_pval_GeneSetEnrichAnalysis_${db}.db_110kb_upstr_40kb_downstr_${suffix}.results > MAGENTA/header.dat
   R -q --no-save < ${PW_location}/MAGENTA/collect.R > collect.log
   cd -
fi

if [ $magma -eq 1 ]; then
   echo "MAGMA"
   if [ ! -d "MAGMA" ]; then
      mkdir MAGMA
   fi
   cd MAGMA
   R -q --no-save < ${PW_location}/MAGMA/data.R
   # Annotation
   magma --annotate window=50,50 --snp-loc magma.snploc --gene-loc $MAGMA/NCBI37.3.gene.loc --out magma
   # Gene analysis - SNP p-values
   magma --bfile $MAGMA/g1000_eur --pval magma.pval ncol=NOBS --gene-annot magma.genes.annot --out magma
   if [ $magenta_db -eq 1 ]; then
      awk -vFS="\t" '{$1=$2;$2=""};1' $PW_location/MAGENTTA/magenta.db | awk '{$2=$2};1'> magenta.db
      cut -f1 magenta.db | awk -vFS="\t" -vOFS="\t" '{print $1)}' > magenta.id
      qsub -V -sync y ${PW_location}/MAGMA/magenta.sh
   elif [ $msigdb_c2 -eq 1 ]; then
      qsub -V -sync y ${PW_location}/MAGMA/c2.sh
   elif [ $msigbdb -eq 1 ]; then
      qsub -V -sync y ${PW_location}/MAGMA/msigdb.sh
   else
      qsub -V -sync y ${PW_location}/MAGMA/depict.sh
   fi
   R -q --no-save < ${PW_location}/MAGMA/collect.R > collect.log
   cd -
fi

if [ $pascal -eq 1 ]; then
   echo "PASCAL"
   if [ ! -d "PASCAL" ]; then
      mkdir PASCAL
   fi
   cd PASCAL
   R -q --no-save < ${PW_location}/PASCAL/data.R
   if [ $magenta_db -eq 1 ]; then
      qsub -V -sync y ${PW_location}/PASCAL/magenta.sh
   elif [ $msigdb_c2 -eq 1 ]; then
      qsub -V -sync y ${PW_location}/PASCAL/c2.sh
   elif [ $msigbdb -eq 1 ]; then
      qsub -V -sync y ${PW_location}/PASCAL/msigdb.sh
   else
      qsub -V -sync y ${PW_location}/PASCAL/depict.sh
   fi
   R -q --no-save < $PW_location}/PASCAL/collect.R > collect.log
   cd -
fi

if [ $depict -eq 1 ]; then
   echo "DEPICT"
   if [ ! -d "DEPICT" ]; then
      mkdir DEPICT
   fi
   cd DEPICT
   R -q --no-save < ${PW_location}/DEPICT/data.R
   if [ $depict_db -eq 1 ]; then
      qsub -V -sync y ${PW_location}/DEPICT/depict.sh
   elif [ $depict_db2 -eq 1 ]; then
      qsub -V -sync y ${PW_location}/DEPICT/depict2.sh
   fi
   cd -
fi

## collection into Excel workbook(s)

if [ $depict_db2 -eq 1]; then
    R -q --no-save < ${PW_location}/files/mmpd.R > mmpd.log
else
    R -q --no-save < ${PW_location}/files/mmp.R > mmp.log
fi

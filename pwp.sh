#!/bin/bash
# 15-11-2017 MRC-Epid JHZ

## SETTINGS

export MAGENTA=/genetics/bin/MAGENTA_software_package_vs2_July2011
export MAGMA=/genetics/bin/MAGMA
export MSigDB=/genetics/src/MSigDB/msigdb_v6.0_GMTs/
export PASCAL=/genetics/bin/PASCAL
export DEPICT=/genetics/bin
export PLINK_EXECUTABLE=/genetics/bin/plink-1.9
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
export depict2=$PASCAL/resources/genesets/depict_discretized_cutoff3.2.gmt
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
   export MAGENTA_db=${PWD}/magenta.db
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
      export db=magenta.db
   elif [ $msigdb_c2 -eq 1 ]; then
      export db=c2.db
      awk '{$2=$1; $1="c2"; print}' $c2 > c2.db
      awk '{$2=$1; $1="c2"};1' FS="\t" OFS="\t" $c2 > c2.db
   elif [ $msigdb -eq 1 ]; then
      export db=msigdb.db
      awk '{$2=$1; $1="msigdb"};1' FS="\t" OFS="\t" $msigdb > msigdb.db
   else
      export db=depict.db
      awk '{FS=OFS="\t";$2=$1;$1="depict";print}' $depict2 > depict.db
   fi
   sed -i 's|magenta.db|'"$db"'|g' magenta.m
   qsub -V -sync y ${PW_location}/MAGENTA/magenta.sh
   export suffix=MAGENTA.db_10000perm_$(date +'%b%d_%y')
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
      awk -vFS="\t" '{$1=$2;$2=""};1' $MAGENTTA_db | awk '{$2=$2};1'> magenta.db
      export MAGMA_DB=$PWD/magenta.db
   elif [ $msigdb_c2 -eq 1 ]; then
      export MAGMA_DB=$c2
   elif [ $msigbdb -eq 1 ]; then
      export MAGMA_DB=$msigdb
   else
      export MAGMA_DB=$depict2
   fi
   qsub -V -sync y ${PW_location}/MAGMA/magma.sh
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
   cp $PW_location/PASCAL/* .
   sed -i 's|OUTPUTDIRECTORY|'"$PWD"'|g' settings.txt
   sed -i 's|PASCAL_location|'"$PASCAL"'|g' settings.txt
   if [ $magenta_db -eq 1 ]; then
      sed -i 's|GENESETFILE|'"$MAGENTA_db"'|g' settings.txt
   elif [ $msigdb_c2 -eq 1 ]; then
      sed -i 's|GENESETFILE|'"$c2"'|g' settings.txt
   elif [ $msigdb -eq 1 ]; then
      sed -i 's|GENESETFILE|'"$MSigDB"'|g' settings.txt
   else
      sed -i 's|GENESETFILE|'"$depict2"'|g' settings.txt
   fi
   qsub -V -sync y ${PW_location}/PASCAL/pascal.sh
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
   cp $PW_location/DEPICT/* .
   sed -i 's|ANALYSIS_PATH|'"$PWD"'|g' depict.cfg
   sed -i 's|PLINK_EXECUTABLE|'"$PLINK_EXECUTABLE"'|g' depict.cfg
   if [ $depict_db -eq 1 ]; then
      sed -i 's|RECONSTITUTED_GENESETS_FILE|data/reconstituted_genesets/reconstituted_genesets_150901.binary|g' depict.cfg
   elif [ $depict_db2 -eq 1 ]; then
      sed -i 's|RECONSTITUTED_GENESETS_FILE|data/reconstituted_genesets/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary|g' depict.cfg
   fi
   qsub -V -sync y ${PW_location}/DEPICT/depict.sh
   cd -
fi

## collection into Excel workbook(s)

if [ $depict_db -eq 1]; then
    R -q --no-save < ${PW_location}/files/mmp.R > mmp.log
elif [ $depict_db2 -eq 1]; then
    R -q --no-save < ${PW_location}/files/mmpd.R > mmpd.log
fi

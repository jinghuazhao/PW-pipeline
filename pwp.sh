#!/bin/bash
# 15-11-2017 MRC-Epid JHZ

## SETTINGS

export MAGENTA=/genetics/bin/MAGENTA_software_package_vs2_July2011
export MAGMA=/genetics/bin/MAGMA
export MSigDB=/genetics/src/MSigDB/msigdb_v6.0_GMTs/
export PASCAL=/genetics/bin/PASCAL
export DEPICT=/genetics/bin/DEPICT/src/python
export PLINK_EXECUTABLE=/genetics/bin/plink-1.9
export PW_location=/genetics/bin/PW-pipeline
export use_UCSC=0

### software flags

export magenta=1
export magma=1
export pascal=1
export depict=1

### database flags (magenta, c2, msigdb, depict) and names

export _db=c2
export c2_db=$MSigDB/msigdb_v6.0_GMTs/c2.all.v6.0.entrez.gmt
export msigdb_db=$MSigDB/msigdb_v6.0_GMTs/msigdb.v6.0.entrez.gmt
export depict_db=$PASCAL/resources/genesets/depict_discretized_cutoff3.2.gmt

## indidivual analyses according to request

R -q --no-save --slave --file=${PW_location}/files/sumstats.R --args ${PWD}/$1

if [ $_db == "magenta" ]; then
   if [ ! -d "MAGENTA" ]; then
      mkdir MAGENTA
   fi
   cd MAGENTA
   for f in $(ls $MAGENTA/*_db); do ln -sf $f; done
   cat GO_terms_BioProc_MolFunc_db Ingenuity_pathways_db KEGG_pathways_db PANTHER_BioProc_db PANTHER_MolFunc_db PANTHER_pathways_db | \
   awk '{$1="magenta"};1' FS="\t" OFS="\t" > magenta.db
   export magenta_db=${PWD}/magenta.db
   cd -
fi
if [ $magenta -eq 1 ]; then
   echo "MAGENTA"
   if [ ! -d "MAGENTA" ]; then
      mkdir MAGENTA
   fi
   cd MAGENTA
   for f in $(find $MAGENTA -name "*"); do ln -sf $f; done
   for f in $(find $PW_location/MAGENTA -name "*"); do ln -sf $f; done
   R -q --no-save < data.R > data.log
   if [ $_db == "magenta" ]; then
      export db=magenta.db
   elif [ $_db == "c2" ]; then
      export db=c2.db
      awk '{$2=$1; $1="c2"};1' FS="\t" OFS="\t" ${c2_db} > c2.db
   elif [ $_db == "msigdb" ]; then
      export db=msigdb.db
      awk '{$2=$1; $1="msigdb"};1' FS="\t" OFS="\t" ${msigdb_db} > msigdb.db
   else
      export db=depict.db
      awk '{FS=OFS="\t";$2=$1;$1="depict";print}' ${depict_db} > depict.db
   fi
   sed -i 's|magenta.db|'"$db"'|g' magenta.m
   qsub -cwd -V -sync y ${PW_location}/MAGENTA/magenta.sh
   export suffix=magenta.db_10000perm_$(date +'%b%d_%y')
   awk '(NR==1){gsub(/\#/,"",$0);print}' $suffix/MAGENTA_pval_GeneSetEnrichAnalysis_${db}.db_110kb_upstr_40kb_downstr_${suffix}.results > header.dat
   R -q --no-save < collect.R > collect.log
   cd -
fi

if [ $magma -eq 1 ]; then
   echo "MAGMA"
   if [ ! -d "MAGMA" ]; then
      mkdir MAGMA
   fi
   cd MAGMA
   R -q --no-save < ${PW_location}/MAGMA/data.R > data.log
   # Annotation
   magma --annotate window=50,50 --snp-loc magma.snploc --gene-loc $MAGMA/NCBI37.3.gene.loc --out magma
   # Gene analysis - SNP p-values
   magma --bfile $MAGMA/g1000_eur --pval magma.pval ncol=NOBS --gene-annot magma.genes.annot --out magma
   if [ $_db == "magenta" ]; then
      awk -vFS="\t" '{$1=$2;$2=""};1' ${magenta_db} | awk '{$2=$2};1'> magenta.db
      export db=magenta.db
   elif [ $_db == "c2" ]; then
      export db=${c2_db}
   elif [ $_db == "msigdb" ]; then
      export db=${msigdb_db}
   else
      export db=${depict_db}
   fi
   qsub -cwd -V -sync y ${PW_location}/MAGMA/magma.sh
   R -q --no-save < collect.R > collect.log
   cd -
fi

if [ $pascal -eq 1 ]; then
   echo "PASCAL"
   if [ ! -d "PASCAL" ]; then
      mkdir PASCAL
   fi
   cd PASCAL
   cp $PW_location/PASCAL/* .
   R -q --no-save < data.R > data.log
   sed -i 's|OUTPUTDIRECTORY|'"$PWD"'|g; s|PASCAL_location|'"$PASCAL"'|g' settings.txt
   if [ $_db == "magenta" ]; then
      sed -i 's|GENESETFILE|'"${magenta_db}"'|g' settings.txt
   elif [ $_db == "c2" ]; then
      sed -i 's|GENESETFILE|'"${c2_db}"'|g' settings.txt
   elif [ $_db == "msigdb" ]; then
      sed -i 's|GENESETFILE|'"${msigdb_db}"'|g' settings.txt
   else
      sed -i 's|GENESETFILE|'"${depict_db}"'|g' settings.txt
   fi
   qsub -cwd -V -N PASCAL_${_db} -sync y ${PW_location}/PASCAL/pascal.sh
   R -q --no-save < collect.R > collect.log
   cd -
fi

if [ $depict -eq 1 ]; then
   echo "DEPICT"
   if [ ! -d "DEPICT" ]; then
      mkdir DEPICT
   fi
   cd DEPICT
   R -q --no-save < ${PW_location}/DEPICT/data.R > data.log
   cp $PW_location/DEPICT/* .
   sed -i 's|ANALYSIS_PATH|'"$PWD"'|g; s|PLINK_EXECUTABLE|'"$PLINK_EXECUTABLE"'|g' depict.cfg
   if [ $_db == "depict" ]; then
      sed -i 's|RECONSTITUTED_GENESETS_FILE|data/reconstituted_genesets/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary|g' depict.cfg
   else
      sed -i 's|RECONSTITUTED_GENESETS_FILE|data/reconstituted_genesets/reconstituted_genesets_150901.binary|g' depict.cfg
   fi
   qsub -cwd -V -sync y ${PW_location}/DEPICT/depict.sh
   cd -
fi

## collection into Excel workbook(s)

if [ $magenta -eq 1 ] && [ $magma -eq 1 ] && [ $pascal -eq 1 ] && [ $depict -eq 1 ] && [ $_db == "depict" ]; then
    R -q --no-save < ${PW_location}/files/mmpd.R > mmpd.log
elif [ $magenta -eq 1 ] && [ $magma -eq 1 ] && [ $pascal -eq 1 ] && [ $depict -eq 1 ]; then
    R -q --no-save < ${PW_location}/files/mmp.R > mmp.log
fi

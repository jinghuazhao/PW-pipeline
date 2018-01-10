#!/bin/bash
# 10-1-2018 MRC-Epid JHZ

## SETTINGS

# multiple precision flag; setting to 1 if needed
export mp=0
export R_LIBS=/genetics/bin/R:/usr/local/lib64/R/library
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/genetics/data/software/lib
export PATH=/genetics/bin/anaconda2/bin:$PATH
export PYTHONPATH=/genetics/bin/anaconda2/lib/python2.7/site-packages

export MAGENTA=/genetics/bin/MAGENTA_software_package_vs2_July2011
export MAGMA=/genetics/bin/MAGMA
export PASCAL=/genetics/bin/PASCAL
export DEPICT=/genetics/bin/DEPICT/src/python
export PLINK_EXECUTABLE=/genetics/bin/plink-1.9
export PW_location=/genetics/bin/PW-pipeline

export MSigDB=/genetics/src/MSigDB/msigdb_v6.0_GMTs
export c2_db=$MSigDB/c2.all.v6.0.entrez.gmt
export msigdb_db=$MSigDB/msigdb.v6.0.entrez.gmt
export depict_discretized=$PASCAL/resources/genesets/depict_discretized_cutoff3.2.gmt

### software flags

export magenta=1
export magma=1
export pascal=1
export depict=1

### min_gs_size for MAGENTA in line with other software

export min_gs_size=5
export max_gs_size=2000

### P value threshold for DEPICT

export p_threshold=0.00000005

### database flag (magenta, c2, msigdb, depict_discretized)

export _db=depict_discretized

if [ $# -lt 1 ] || [ "$args" == "-h" ]; then
    echo "Usage: pwp.sh <input>"
    echo "where <input> is in sumstats format:"
    echo "SNP A1 A2 freqA1 beta se P N chr pos"
    echo "where SNP is RSid, A1 is effect allele"
    echo "and the outputs will be in MAGENTA/MAGMA/PASCAL/DEPICT directory"
    exit
fi

## indidivual analyses according to request

sed 's/\t/ /g' $1 > $1.sumstats
export sumstats=$1.sumstats
R -q --no-save <${PW_location}/files/sumstats.R
export sumstats_rda=${PWD}/sumstats.rda

if [ $_db == "magenta" ]; then
   if [ ! -d "MAGENTA" ]; then
      mkdir MAGENTA
   fi
   cd MAGENTA
   for f in $(ls $MAGENTA/*_db); do ln -sf $f; done
   cat GO_terms_BioProc_MolFunc_db Ingenuity_pathways_db KEGG_pathways_db PANTHER_BioProc_db PANTHER_MolFunc_db PANTHER_pathways_db | \
   awk '{$1="magenta";gsub(/ /,"_",$2);$2=NR ":" $2};1' FS="\t" OFS="\t" > magenta.db
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
   R -q --no-save < data.R > ${_db}.data.log
   if [ $_db == "magenta" ]; then
      export db=magenta.db
   elif [ $_db == "c2" ]; then
      export db=c2.db
      awk '{$2=$1; $1="c2"};1' FS="\t" OFS="\t" ${c2_db} > $db
   elif [ $_db == "msigdb" ]; then
      export db=msigdb.db
      awk '{$2=$1; $1="msigdb"};1' FS="\t" OFS="\t" ${msigdb_db} > $db
   else
      export db=$(basename $depict_discretized .gmt)
      awk '{$2=$1;$1="depict";print}' FS="\t" OFS="\t" ${depict_discretized} > $db
   fi
   sed -i 's|GWAS_SNP_SCORE_FILE_NAME|magenta|g' Run_MAGENTA_vs2_July_2011.m
   sed -i 's|GENESET_DB_FILE_NAME|'"$db"'|g' Run_MAGENTA_vs2_July_2011.m
   sed -i 's|MIN_GS_SIZE|'"$min_gs_size"'|g' Run_MAGENTA_vs2_July_2011.m
   sed -i 's|MAX_GS_SIZE|'"$max_gs_size"'|g' Run_MAGENTA_vs2_July_2011.m
   export suffix=_10000perm_$(date +'%b%d_%y')
   qsub -cwd -N MAGENTA_${db} -V -sync y ${PW_location}/MAGENTA/magenta.sh
   awk '(NR==1){gsub(/\#/,"",$0);print}' ${db}${suffix}/MAGENTA_pval_GeneSetEnrichAnalysis_${db}_110kb_upstr_40kb_downstr${suffix}.results > ${db}.dat
#  sed -i 's/[[:digiti:]]\+\://g' ${db}${suffix}/MAGENTA_pval_GeneSetEnrichAnalysis_${db}_110kb_upstr_40kb_downstr${suffix}.results
   R -q --no-save < collect.R > ${_db}.collect.log
   $PW_location/files/network.sh magenta
   cd -
fi

if [ $magma -eq 1 ]; then
   echo "MAGMA"
   if [ ! -d "MAGMA" ]; then
      mkdir MAGMA
   fi
   cd MAGMA
   R -q --no-save < ${PW_location}/MAGMA/data.R > ${_db}.data.log
   # Annotation
   magma --annotate window=50,50 --snp-loc ${_db}.snploc --gene-loc $MAGMA/NCBI37.3.gene.loc --out ${_db}
   # Gene analysis - SNP p-values
   magma --bfile $MAGMA/g1000_eur --pval ${_db}.pval ncol=NOBS --gene-annot ${_db}.genes.annot --out ${_db}
   if [ $_db == "magenta" ]; then
      awk '{$1=$2;$2=""};1' FS="\t" ${magenta_db} | awk '{$2=$2};1'> magenta.db
      export db=magenta.db
   elif [ $_db == "c2" ]; then
      export db=${c2_db}
   elif [ $_db == "msigdb" ]; then
      export db=${msigdb_db}
   else
      export db=$(basename $depict_discretized .gmt)
      ln -sf ${depict_discretized} $db
   fi
   # Gene-set analysis
   qsub -cwd -N MAGMA_${_db} -V -sync y ${PW_location}/MAGMA/magma.sh
   export db=$(basename $db)
   R -q --no-save < ${PW_location}/MAGMA/sets.R > ${_db}.sets.log
   $PW_location/files/network.sh magma
   cd -
fi

if [ $pascal -eq 1 ]; then
   echo "PASCAL"
   if [ ! -d "PASCAL" ]; then
      mkdir PASCAL
   fi
   cd PASCAL
   cp $PW_location/PASCAL/* .
   R -q --no-save < data.R > ${_db}.data.log
   sed -i 's|OUTPUTDIRECTORY|'"$PWD"'|g; s|PASCAL_location|'"$PASCAL"'|g' settings.txt
   if [ $_db == "magenta" ]; then
      awk '{$1=$2;$2="."};1' FS="\t" ${magenta_db} > magenta.db
      export db=${PWD}/magenta.db
      sed -i 's|GENESETFILE|'"${db}"'|g' settings.txt
   elif [ $_db == "c2" ]; then
      export db=$(basename $c2_db .gmt)
      sed -i 's|GENESETFILE|'"${c2_db}"'|g' settings.txt
   elif [ $_db == "msigdb" ]; then
      export db=$(basename $msigdb_db .gmt)
      sed -i 's|GENESETFILE|'"${msigdb_db}"'|g' settings.txt
   else
      export db=$(basename $depict_discretized .gmt)
      sed -i 's|GENESETFILE|'"${depict_discretized}"'|g' settings.txt
   fi
   qsub -cwd -V -N PASCAL_${_db} -sync y ${PW_location}/PASCAL/pascal.sh
   R -q --no-save < collect.R > ${_db}.collect.log
   $PW_location/files/network.sh pascal
   cd -
fi

if [ $depict -eq 1 ]; then
   echo "DEPICT"
   if [ ! -d "DEPICT" ]; then
      mkdir DEPICT
   fi
   cd DEPICT
   cp $PW_location/DEPICT/* .
   R -q --no-save < data.R > ${_db}.data.log
   sed -i 's|ANALYSIS_PATH|'"$PWD"'|g; s|PLINK_EXECUTABLE|'"$PLINK_EXECUTABLE"'|g' depict.cfg
   if [ $_db == "depict_discretized" ]; then
      export db=$(basename $depict_discretized .gmt)
      sed -i 's|RECONSTITUTED_GENESETS_FILE|data/reconstituted_genesets/reconstituted_genesets_150901.binary|g' depict.cfg
      sed -i 's|LABEL_FOR_OUTPUT_FILES|depict_discretized_cutoff3.2|g' depict.cfg
   else
      export db=depict
      sed -i 's|RECONSTITUTED_GENESETS_FILE|data/reconstituted_genesets/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary|g' depict.cfg
      sed -i 's|LABEL_FOR_OUTPUT_FILES|depict|g' depict.cfg
   fi
   qsub -cwd -N DEPICT -V -sync y ${PW_location}/DEPICT/depict.sh
   bash tissue_plot.sh $db
   R -q --no-save < ${PW_location}/DEPICT/collect.R > ${_db}.collect.log
   $PW_location/files/network.sh depict
   cd -
fi

## collection into Excel workbook(s)

if [ $magenta -eq 1 ] && [ $magma -eq 1 ] && [ $pascal -eq 1 ] && [ $depict -eq 1 ] && [ $_db == "depict_discretized" ]; then
    R -q --no-save < ${PW_location}/files/mmpd.R > ${_db}.mmpd.log
    R -q --no-save < ${PW_location}/files/summary.R > ${_db}.summary.log
elif [ $magenta -eq 1 ] && [ $magma -eq 1 ] && [ $pascal -eq 1 ]; then
    R -q --no-save < ${PW_location}/files/mmp.R > ${_db}.mmp.log
fi

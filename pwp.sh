#!/bin/bash
# 17-9-2018 MRC-Epid JHZ

## SETTINGS

source pwp.ini

## ANALYSIS

# functions

function fdr_cutoff()
{
   sort -k1,1 | \
   awk '{
       FS=OFS="\t";
       fdr=$3;
       if(fdr>=0.2) { $3=">=0.2" } else { 
         $3="<0.2";
         if(fdr<0.05) $3="<0.05";
         if(fdr<0.01) $3="<0.01";
       }
       print
    }' > ${db}.txt; \ 
    echo -e "Original gene set ID\tOriginal gene set description\tNominal P value\tFalse discovery rate" > ${db}_genesetenrichment.txt; \
    gunzip -c $PW_location/files/id_descrip.txt.gz | awk 'NR>1' | sort -k1,1 | join -t $'\t' -j1 - ${db}.txt >> ${db}_genesetenrichment.txt
}

function network_plot()
{
     cp $PW_location/files/network_plot*.cfg $PW_location/files/network_plot_CytoscapeStyle_v1.xml .
   # minor changes to network_plot.py are necessary.
     export file_genesetenrichment=${db}_genesetenrichment.txt
     sed -i 's|FILE_GENESETENRICHMENT|'"$file_genesetenrichment"'|g' network_plot.cfg
     sed -i 's|OUTPUT_LABEL|'"$db"'|g' network_plot.cfg
     sed -i 's|CUTOFF_TYPE|'"$cutoff_type"'|g' network_plot.cfg
     sed -i 's|PVALUE_CUTOFF|'"$pvalue_cutoff"'|g' network_plot.cfg
     sed -i 's|CYTOSCAPE_LOC|'"$CYTOSCAPE"'|g' network_plot.cfg
     $PW_location/files/network_plot.py network_plot.cfg
   # The old version requires addtional changes as follows,
   # sed 's/flag_interactive_cytoscape_session/interactive_cytoscape_session/g' network_plot.cfg > network_plot_2015.cfg
   # sed -i 's|output_label: ./'"$db"'|output_label: network_plot_2015/'"$db"'|g' network_plot_2015.cfg
   # $PW_location/files/network_plot_2015.py network_plot_2015.cfg
   # pdftopng -r 300 ${db}_network_diagram.pdf ${db}_network_diagram
   # mv ${db}_network_diagram-000001.png ${db}_network_diagram.png
}

if [ $collection_only -eq 0 ]; then
   if [ $# -lt 1 ] || [ "$args" == "-h" ]; then
       echo "Usage: pwp.sh <input>"
       echo "where <input> is in sumstats format:"
       echo "SNP A1 A2 freqA1 beta se P N chr pos"
       echo "where SNP is RSid, A1 is effect allele"
       echo "and the outputs will be in MAGENTA/MAGMA/PASCAL/DEPICT directory"
       exit
   fi

 # individual analyses according to request

   sed 's/\t/ /g' $1 > $1.sumstats
   export sumstats=$1.sumstats
   R -q --no-save <${PW_location}/files/sumstats.R
   export sumstats_rda=${PWD}/sumstats.rda
fi

if [ $_db == "magenta" ]; then
   if [ ! -d "MAGENTA" ]; then
      mkdir MAGENTA
   fi
   cd MAGENTA
   if [ $collection_only -eq 0 ]; then
      for f in $(ls $MAGENTA/*_db); do ln -sf $f; done
      cat GO_terms_BioProc_MolFunc_db Ingenuity_pathways_db KEGG_pathways_db PANTHER_BioProc_db PANTHER_MolFunc_db PANTHER_pathways_db | \
      awk '{$1="magenta";gsub(/ /,"_",$2);$2=NR ":" $2};1' FS="\t" OFS="\t" > magenta.db
   fi
   export magenta_db=${PWD}/magenta.db
   cd -
fi

# Meta-Analysis Gene-set Enrichment of variaNT Associations (MAGENTA)

if [ $magenta -eq 1 ]; then
   echo "MAGENTA"
   if [ ! -d "MAGENTA" ]; then
      mkdir MAGENTA
   fi
   cd MAGENTA
   if [ $collection_only -eq 0 ]; then
      for f in $(find $MAGENTA -name "*"); do ln -sf $f; done
      for f in $(find $PW_location/MAGENTA -name "*"); do ln -sf $f; done
      R -q --no-save < data.R > ${_db}.data.log
   fi
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
   if [ $collection_only -eq 0 ]; then
      sed -i 's|GWAS_SNP_SCORE_FILE_NAME|magenta|g' Run_MAGENTA_vs2_July_2011.m
      sed -i 's|GENESET_DB_FILE_NAME|'"$db"'|g' Run_MAGENTA_vs2_July_2011.m
      sed -i 's|MIN_GS_SIZE|'"$min_gs_size"'|g' Run_MAGENTA_vs2_July_2011.m
      sed -i 's|MAX_GS_SIZE|'"$max_gs_size"'|g' Run_MAGENTA_vs2_July_2011.m
      export suffix=_10000perm_$(date +'%b%d_%y')
      if [ $use_seq -eq 1 ]; then 
         qsub -cwd -N MAGENTA_${db} -V -sync y ${PW_location}/MAGENTA/magenta.sh
      else
         ${PW_location}/MAGENTA/magenta.sh
      fi
      awk '(NR==1){gsub(/\#/,"",$0);print}' ${db}${suffix}/MAGENTA_pval_GeneSetEnrichAnalysis_${db}_110kb_upstr_40kb_downstr${suffix}.results > ${db}.dat
      #  sed -i 's/[[:digiti:]]\+\://g' ${db}${suffix}/MAGENTA_pval_GeneSetEnrichAnalysis_${db}_110kb_upstr_40kb_downstr${suffix}.results
      if [ $_db == "depict_discretized" ]; then
         $PW_location/files/network.sh magenta $DEPICT
         cut -f2,10,11 ${db}.dat | awk 'NR>1' | fdr_cutoff
         network_plot
      fi
      R -q --no-save < collect.R > ${_db}_collect.log
   fi
   cd -
fi

# Multi-marker Analysis of GenoMic Annotation (MAGMA)

if [ $magma -eq 1 ]; then
   echo "MAGMA"
   if [ ! -d "MAGMA" ]; then
      mkdir MAGMA
   fi
   cd MAGMA
   if [ $collection_only -eq 0 ]; then
      R -q --no-save < ${PW_location}/MAGMA/data.R > ${_db}.data.log
      # Annotation
      magma --annotate window=50,50 --snp-loc ${_db}.snploc --gene-loc $MAGMA/NCBI37.3.gene.loc --out ${_db}
      # Gene analysis - SNP p-values
      magma --bfile $MAGMA/g1000_eur --pval ${_db}.pval ncol=NOBS --gene-annot ${_db}.genes.annot --out ${_db}
   fi
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
   if [ $collection_only -eq 0 ]; then
      if [ $use_sge -eq 1 ]; then
         qsub -cwd -N MAGMA_${_db} -V -sync y ${PW_location}/MAGMA/magma.sh
      else
         ${PW_location}/MAGMA/magma.sh
      fi
      export db=$(basename $db)
      if [ $_db == "depict_discretized" ]; then
         $PW_location/files/network.sh magma $DEPICT
         fdr_cutoff
         network_plot
      fi
      R -q --no-save < ${PW_location}/MAGMA/sets.R > ${_db}.sets.log
   fi
   cd -
fi

# PAthway SCoring ALgorithm (PASCAL)

if [ $pascal -eq 1 ]; then
   echo "PASCAL"
   if [ ! -d "PASCAL" ]; then
      mkdir PASCAL
   fi
   cd PASCAL
   if [ $collection_only -eq 0 ]; then
      cp $PW_location/PASCAL/* .
      R -q --no-save < data.R > ${_db}.data.log
      sed -i 's|OUTPUTDIRECTORY|'"$PWD"'|g; s|PASCAL_location|'"$PASCAL"'|g' settings.txt
   fi
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
   if [ $collection_only -eq 0 ]; then
      if [ $use_sge -eq 1 ]; then
         qsub -cwd -V -N PASCAL_${_db} -sync y ${PW_location}/PASCAL/pascal.sh
      else
         ${PW_location}/PASCAL/pascal.sh
      fi
      if [ $_db == "depict_discretized" ]; then
         $PW_location/files/network.sh pascal $DEPICT
         fdr_cutoff
         network_plot
      fi
      R -q --no-save < collect.R > ${_db}_collect.log
   fi
   cd -
fi

# Data-Driven Expression Prioritized Integration for Complex Traits (DEPICT)

if [ $depict -eq 1 ]; then
   echo "DEPICT"
   if [ ! -d "DEPICT" ]; then
      mkdir DEPICT
   fi
   cd DEPICT
   if [ $collection_only -eq 0 ]; then
      cp $PW_location/DEPICT/* .
      R -q --no-save < data.R > ${_db}.data.log
      sed -i 's|ANALYSIS_PATH|'"$PWD"'|g; s|PLINK_EXECUTABLE|'"$PLINK_EXECUTABLE"'|g' depict.cfg
   fi
   if [ $_db == "depict" ]; then
      export db=depict
      sed -i 's|RECONSTITUTED_GENESETS_FILE|data/reconstituted_genesets/reconstituted_genesets_150901.binary|g' depict.cfg
      sed -i 's|LABEL_FOR_OUTPUT_FILES|depict|g' depict.cfg
   elif [ $_db == "depict_discretized" ]; then
      export db=$(basename $depict_discretized .gmt)
      sed -i 's|RECONSTITUTED_GENESETS_FILE|data/reconstituted_genesets/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.binary|g' depict.cfg
      sed -i 's|LABEL_FOR_OUTPUT_FILES|depict_discretized_cutoff3.2|g' depict.cfg
   fi
   if [ $collection_only -eq 0 ]; then
      sed -i 's|NUMBER_OF_THREADS|'"$number_of_threads"'|g' depict.sh
      sed -i 's|NUMBER_OF_THREADS|'"$number_of_threads"'|g' depict.cfg
      sed -i 's|ASSOCIATION_PVALUE_CUTOFF|'"$p_threshold"'|g' depict.cfg
      sed -i 's|NR_REPITITIONS|'"$nr_repititions"'|g' depict.cfg
      if [ $use_sge -eq 1 ]; then
         qsub -cwd -N DEPICT -V -sync y ./depict.sh
      else
         ./depict.sh
      fi
      bash tissue_plot.sh $db
      if [ _db == "depict" ] || [ $_db == "depict_discretized" ]; then
         $PW_location/files/network.sh depict $DEPICT
      fi
      network_plot
      R -q --no-save < collect.R > ${_db}_collect.log
   fi
   cd -
fi

# collection into Excel workbook(s)

if [ $magenta -eq 1 ] && [ $magma -eq 1 ] && [ $pascal -eq 1 ] && [ $depict -eq 1 ] && [ $_db == "depict_discretized" ]; then
    R -q --no-save < ${PW_location}/files/mmpd.R > ${_db}.mmpd.log
    R -q --no-save < ${PW_location}/files/summary.R > ${_db}.summary.log
elif [ $magenta -eq 1 ] && [ $magma -eq 1 ] && [ $pascal -eq 1 ]; then
    R -q --no-save < ${PW_location}/files/mmp.R > ${_db}.mmp.log
fi

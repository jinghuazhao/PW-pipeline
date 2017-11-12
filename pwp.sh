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
   R -q --no-save < ${PW_location}/MAGENTA/magenta.R
   if [ $magenta_db -eq 1 ]; then
      qsub ${PW_location}/MAGENTA/magenta.sh
   elif [ $msigdb_c2 -eq 1 ]; then
      qsub ${PW_location}/MAGENTA/c2.sh
   elif [ $msigbdb -eq 1 ]; then
      qsub ${PW_location}/MAGENTA/msigdb.sh
   else
      qsub ${PW_location}/MAGENTA/depict2.sh
   fi
   R -q --no-save < ${PW_location}/MAGENTA/collect.R > MAGENTA/collect.log
fi

if [ $magma -eq 1 ]; then
   echo "MAGMA"
   if [ ! -d "MAGMA" ]; then
      mkdir MAGMA
   fi
   R -q --no-save < ${PW_location}/MAGMA/magma.R
   if [ $magenta_db -eq 1 ]; then
      qsub ${PW_location}/MAGMA/magenta.sh
   elif [ $msigdb_c2 -eq 1 ]; then
      qsub ${PW_location}/MAGMA/c2.sh
   elif [ $msigbdb -eq 1 ]; then
      qsub ${PW_location}/MAGMA/msigdb.sh
   else
      qsub ${PW_location}/MAGMA/depict2.sh
   fi
   R -q --no-save < ${PW_location}/MAGMA/collect.R > MAGMA/collect.log
fi

if [ $pascal -eq 1 ]; then
   echo "PASCAL"
   if [ ! -d "PASCAL" ]; then
      mkdir PASCAL
   fi
   R -q --no-save < ${PW_location}/PASCAL/pascal.R
   if [ $magenta_db -eq 1 ]; then
      qsub ${PW_location}/PASCAL/magenta.sh
   elif [ $msigdb_c2 -eq 1 ]; then
      qsub ${PW_location}/PASCAL/c2.sh
   elif [ $msigbdb -eq 1 ]; then
      qsub ${PW_location}/PASCAL/msigdb.sh
   else
      qsub ${PW_location}/PASCAL/depict2.sh
   fi
   R -q --no-save < $PW_location}/PASCAL/collect.R > PASCAL/collect.log
fi

if [ $depict -eq 1 ]; then
   echo "DEPICT"
   if [ ! -d "DEPICT" ]; then
      mkdir DEPICT
   fi
   R -q --no-save < ${PW_location}/DEPICT/depict.R
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

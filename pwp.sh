# 10-11-2017 MRC-Epid JHZ

## settings

export PW_location=/genetics/bin/PW-pipeline

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

if [ $magenta -eq 1 ]; then
   echo "MAGENTA"
   if [ $magenta_db -eq 1 ]; then
      qsub ${PW_location}/MAGENTA/matlab.sh
   elif [ $msigdb_c2 -eq 1 ]; then
      qsub ${PW_location}/MAGENTA/c2.sh
   elif [ $msigbdb -eq 1 ]; then
      qsub ${PW_location}/MAGENTA/msigdb.sh
   else
      qsub ${PW_location}/MAGENTA/depict.sh
   fi
fi

if [ $magma -eq 1 ]; then
   echo "MAGMA"
   if [ $magenta_db -eq 1 ]; then
      qsub ${PW_location}/MAGENTA/matlab.sh
   elif [ $msigdb_c2 -eq 1 ]; then
      qsub ${PW_location}/MAGENTA/c2.sh
   elif [ $msigbdb -eq 1 ]; then
      qsub ${PW_location}/MAGENTA/msigdb.sh
   else
      qsub ${PW_location}/MAGENTA/depict.sh
   fi
fi

if [ $pascal -eq 1 ]; then
   echo "PASCAL"
   if [ $magenta_db -eq 1 ]; then
      qsub ${PW_location}/MAGENTA/matlab.sh
   elif [ $msigdb_c2 -eq 1 ]; then
      qsub ${PW_location}/MAGENTA/c2.sh
   elif [ $msigbdb -eq 1 ]; then
      qsub ${PW_location}/MAGENTA/msigdb.sh
   else
      qsub ${PW_location}/MAGENTA/depict.sh
   fi
fi

if [ $depict -eq 1 ]; then
   echo "DEPICT"
   if [ $depict_db -eq 1 ]; then
      qsub ${PW_location}/MAGENTA/matlab.sh
   elif [ $depict_db2 -eq 1 ]; then
      qsub ${PW_location}/MAGENTA/c2.sh
   fi
fi

## collection into Excel

R -q --no-save < ${PW_location}/files/mmp.R > mmp.log

if [ $depict -eq 1]; then
    R -q --no-save < ${PW_location}/files/mmpd.R > mmpd.log
fi


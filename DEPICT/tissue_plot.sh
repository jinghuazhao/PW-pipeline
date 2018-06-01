# 1-6-2018 MRC-Epid JHZ

# R CMD BATCH --no-save --no-restore '--args ${1}_tissueenrichment.txt' tissue_plot.R ${1}.out

Rscript tissue_plot.R ${1}_tissueenrichment.txt $1

# pdf --> png conversion via xpdf to facilitate generation of Excel workbook

for p in cells multiplot system tissues
do
   r=${1}_genenetwork_${p}
   if [ -f ${r}.png ]; then
      rm ${r}.png
   fi
   pdftopng -r 300 ${r}.pdf ${r}
   if [ -f ${r}-000001.png ]; then
      mv ${r}-000001.png ${r}.png
   fi
done

# conversion via ImageMagick:

#export prefix=${1}_genenetwork_
#for type in cells multiplot system tissues
#do
#   echo Converting ${prefix}{type} ...
#   convert -units PixelsPerInch ${prefix}${type}.pdf -density 300 ${prefix}${type}.png
#done

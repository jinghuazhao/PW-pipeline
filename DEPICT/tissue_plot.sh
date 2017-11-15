# 15-11-2017 MRC-Epid JHZ

# R CMD BATCH --no-save --no-restore '--args depict_tissueenrichment.txt' tissue_plot.R depict.out

Rscript tissue_plot.R depict_tissueenrichment.txt depict

# pdf --> png conversion via xpdf to facilitate generation of Excel workbook

for p in cells multiplot system tissues
do
   r=tissue_plot_depict_genenetwork_${p}
   if [ -f ${r}.png ]; then
      rm ${r}.png
   fi
   pdftopng ${r}.pdf -r 300 ${r}
   if [ -f ${r}-000001.png ]; then
      mv ${r}-000001.png ${r}.png
   fi
done

exit

# conversion via ImageMagick:

export prefix=tissue_plot_depict_genenetwork_
for type in cells multiplot system tissues
do
   echo Converting ${prefix}{type} ...
   convert -units PixelsPerInch ${prefix}${type}.pdf -density 300 ${prefix}${type}.png
done

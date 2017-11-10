# 5-7-2017 MRC-Epid JHZ

# R CMD BATCH --no-save --no-restore '--args depict_tissueenrichment.txt' tissue_plot.R depict.out

Rscript tissue_plot.R depict_tissueenrichment.txt depict

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

#convert tissue_plot_depict_genenetwork_cells.pdf tissue_plot_depict_genenetwork_cells.png
#convert tissue_plot_depict_genenetwork_multiplot.pdf tissue_plot_depict_genenetwork_multiplot.png
#convert tissue_plot_depict_genenetwork_system.pdf tissue_plot_depict_genenetwork_system.png
#convert tissue_plot_depict_genenetwork_tissues.pdf tissue_plot_depict_genenetwork_tissues.png


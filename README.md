# PW-pipeline

PathWay pipeline using GWAS summary statistics, named analogously after FM-pipepline I have implemented.

## INTRODUCTION

Pathway analysis becomes an important element in GWAS. Broadly, it involves SNP annotation, such as Variant Effect Predictor 
(VEP), gene analysis such as VEGAS2, and gene set analysis. Visualisation of a particular region has been facilitated with 
LocusZoom, while network(s) from pathway analysis via [gephi](https://gephi.org/) or [Cytoscape](http://www.cytoscape.org/), 
which uses genes and a collection of edges, directed or undirected, to build a network. Aspects to consider include part or all databases, 
individual level genotype data vs GWAS summary statistics, computing speed, with and without tissue enrichment.

![diagram from CytoScape/GeneMANIA](files/obesity.png)

## INSTALLATION

This pipeline involves several software for pathway analysis using GWAS summary statistics, as shown below,

**Full name** | **Abbreviation** | **Reference**
----------------------------------------------------------|--------------|----------
Meta-Analysis Gene-set Enrichment of variaNT Associations | MAGENTA | Segre, et al. (2010)
Multi-marker Analysis of GenoMic Annotation | MAGMA | de Leeuw, et al. (2015)
PAthway SCoring ALgorithm | PASCAL | Lamparter, et al. (2016)
Data-Driven Expression Prioritized Integration for Complex Traits | DEPICT | Pers, et al.(2015)

The full functionality of the pipeline requires availability of individual software for pathway analysis, whose requirements should be
fulfiled, e.g., [Matlab](https://www.mathworks.com/products/matlab.html) for MAGENTA. For PASCAL, minor changes
need to be made with Pascal in that by default DIR is where it is called so it needs to be changed into PASCAL installation directory 
and that instead of the relative path jars/pascalDeployed.jar an absolute prefix should be added.
The current version of pipeline also uses DEPICT from the GitHub but with data in the release version from the Broad,
[https://data.broadinstitute.org/mpg/depict/](https://data.broadinstitute.org/mpg/depict/depict_140721.tar.bz2). It is useful to install
[XpdfReader](https://www.xpdfreader.com/) or [ImageMagick](https://www.imagemagick.org/) to produce Excel workbook. By 
default [Sun grid engine](https://en.wikipedia.org/wiki/Oracle_Grid_Engine) is used but this can be any other mechanism such 
as [GNU parallel](https://www.gnu.org/software/parallel/) [note with its --env to pass environment variables]. As usual, 
[R](https://www.r-project.org/) is required.

On systems supporting modules, they can be loaded before hand but it is possible that appropriate module is loaded seamlessly, e.g.,
```bash
echo -e "function module (){eval \`/usr/bin/modulecmd bash $*\`}" > matlab
echo module load matlab/r2017b >> matlab
echo matlab \$\@ >> matlab
chmod +x matlab
```
NB "$@" enables matlab to be called as usual. Alternatively, the module coammands can be part of pwp.ini which is sourced with pwp.sh.

The pipeline itself can be installed from GitHub in the usual way.
```
git clone https://github.com/jinghuazhao/PW-pipeline
```

## USAGE

The pipeline requires users to specify both software and database to be used, which is now through [pwp.ini](pwp.ini) in the working directory.

The syntax is then
```bash
bash pwp.sh <input file>
```
### Input

The input will be GWAS summary statistics described at https://github.com/jinghuazhao/SUMSTATS **in that order without the header**,

### Output

The output will be available from individual directories named after the software you choose, and optionally in case all software are used the output can also be
an Excel workbook containing combined results.

For DEPICT databases, it is possible to call [netowrk_plot.py](files/network_plot.py) to generate network diagram and perform cluster analysis.

## EXAMPLES

The `BMI` and `ST4` from https://github.com/jinghuazhao/SUMSTATS can be called as follows,
```bash
gunzip -c bmi.tsv.gz > BMI
pwp.sh BMI
```
and
```
pwp.sh ST4 &
```

## ADDITIONAL TOPICS

The [wiki page](https://github.com/jinghuazhao/PW-pipeline/wiki) contains the following information,

* [Databases](https://github.com/jinghuazhao/PW-pipeline/wiki/Databases)
* [Features](https://github.com/jinghuazhao/PW-pipeline/wiki/Features)
* [Tissue and network plots](https://github.com/jinghuazhao/PW-pipeline/wiki/Tissue-and-network-plots)
* [Result collection](https://github.com/jinghuazhao/PW-pipeline/wiki/Result-collection)

You can make changes to the configuration files for each software in their own direcories. See also [software-notes](https://github.com/jinghuazhao/software-notes) on how to set up VEGAS2, as in [vegas2v2.sh](vegas2v2.sh).

## RELATED LINKS

* [BioGRID](https://thebiogrid.org/): an interaction repository with data compiled through comprehensive curation efforts.
* [Osprey](https://osprey.thebiogrid.org/): Network Visualization System.
* [GeneMANIA](http://genemania.org/): Imports interaction networks from public databases from a list of genes with their annotations and putative functions.
* [rGREAT](https://github.com/jokergoo): Client for GREAT Analysis
* [VisANT](http://visant.bu.edu/): Visual analyses of metabolic networks in cells and ecosystems.

## ACKNOWLEDGEMENTS

The work drives from comparison of software performances using our own GWAS data. The practicality of a common DEPICT database to all software 
here was due to PASCAL developer(s). At the end of our implementation it came to our attention that similar effort has been made, e.g., 
[DEPICT-pipeline](https://github.com/RebeccaFine/DEPICT-pipeline) and other adaptations.

## SOFTWARE AND REFERENCES

**[DEPICT](https://data.broadinstitute.org/mpg/depict/)** (**[GitHub](https://github.com/perslab/depict)**)

Pers TH et al.(2015) Biological interpretation of genome-wide association studies using predicted gene functions. *Nat Commun* 6:5890. doi: 10.1038/ncomms6890.

**[MAGENTA](https://software.broadinstitute.org/mpg/magenta/)**

Segre AV, et al (2010). Common Inherited Variation in Mitochondrial Genes Is Not Enriched for Associations with Type 2 Diabetes or Related Glycemic Traits. *PLoS 
Genet* 12;6(8). pii: e1001058. doi: 10.1371/journal.pgen.1001058

**[MAGMA](http://ctg.cncr.nl/software/magma)**

de Leeuw C, et al. (2015) MAGMA: Generalized Gene-Set Analysis of GWAS Data. *PLoS Comput Biol* 11(4): e1004219. doi:  10.1371/journal.pcbi.1004219

**[PASCAL](http://www2.unil.ch/cbg/images/3/3d/PASCAL.zip)** (**[GitHub](https://github.com/dlampart/Pascal)**)

Lamparter D, et al. (2016) Fast and Rigorous Computation of Gene and Pathway Scores from SNP-Based Summary Statistics. *PLoS Comput Biol* 12(1):e1004714. 
doi: 10.1371/journal.pcbi.1004714

**[VEGAS2](https://vegas2.qimrberghofer.edu.au/)** (Versatile Gene-based Association Study)

Liu JZ, et al. (2010). A versatile gene-based test for genome-wide association studies. *Am J Hum Genet* 87:139–145.

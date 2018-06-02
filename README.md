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

The full functionality of the pipeline requires availability of individual software for pathway analysis, which should fulfil 
their requirements, e.g., [Matlab](https://www.mathworks.com/products/matlab.html) for MAGENTA, PLINK. It is useful to 
install [xpdf](https://www.xpdfreader.com/) or [ImageMagick](https://www.imagemagick.org/) to produce Excel workbook. By 
default [Sun grid engine](https://en.wikipedia.org/wiki/Oracle_Grid_Engine) is used but this can be any other mechanism such 
as [GNU parallel](https://www.gnu.org/software/parallel/) [note with its --env to pass environment variables]. As usual, 
[R](https://www.r-project.org/) is required.

The pipeline itself can be installed from GitHub in the usual way.
```
git clone https://github.com/jinghuazhao/PW-pipeline
```
Features and databases are described at the repository's [wiki page](https://github.com/jinghuazhao/PW-pipeline/wiki).

## USAGE

The pipeline requires users to specify both software and database to be used. It is possible 
that a given database can be used for several software when appropriate.

The syntax is
```
bash pwp.sh <input file>
```
### Input

The input will be GWAS summary statistics described at [SUMSTATS](https://github.com/jinghuazhao/SUMSTATS) **in that order without the header**,

### Output

The output will be available from individual directories named after the software you choose, and optionally in case all software are used the output can also be
an Excel workbook containing combined results.

In case of DEPICT database, it is possible to call [netowrk_plot.py](files/network_plot.py) to generate network diagram and perform cluster analysis.

## EXAMPLES

The `bmi.txt` and `ST4` from SUMSTATS can be called as follows,
```
pwp.sh bmi.txt
```
and
```
pwp.sh ST4 &
```

## ACKNOWLEDGEMENTS

The work drives from comparison of software performances using our own GWAS data. The practicality of a common DEPICT database to all software 
here was due to PASCAL developer(s). At the end of our implementation it came to our attention that similar effort has been made, e.g., 
[DEPICT-pipeline](https://github.com/RebeccaFine/DEPICT-pipeline) and other adaptations.

## RELATED LINKS

* [BioGRID](https://thebiogrid.org/): an interaction repository with data compiled through comprehensive curation efforts.
* [Osprey](https://osprey.thebiogrid.org/): Network Visualization System.
* [GeneMANIA](http://genemania.org/): Imports interaction networks from public databases from a list of genes with their annotations and putative functions.
* [rGREAT](https://github.com/jokergoo): Client for GREAT Analysis
* [VisANT](http://visant.bu.edu/): Visual analyses of metabolic networks in cells and ecosystems.

## SOFTWARE AND REFERENCES

[DEPICT](https://data.broadinstitute.org/mpg/depict/) ([GitHub](https://github.com/perslab/depict))

Pers TH et al.(2015) Biological interpretation of genome-wide association studies using predicted gene functions. Nat Commun. 6:5890. doi: 10.1038/ncomms6890.

[MAGENTA](https://software.broadinstitute.org/mpg/magenta/)

Segre AV, et al (2010). Common Inherited Variation in Mitochondrial Genes Is Not Enriched for Associations with Type 2 Diabetes or Related Glycemic Traits. PLoS 
Genet. 12;6(8). pii: e1001058. doi: 10.1371/journal.pgen.1001058

[MAGMA](http://ctg.cncr.nl/software/magma)

de Leeuw C, et al. (2015) MAGMA: Generalized Gene-Set Analysis of GWAS Data. PLoS Comput Biol. 11(4): e1004219. doi:  10.1371/journal.pcbi.1004219

[PASCAL](http://www2.unil.ch/cbg/images/3/3d/PASCAL.zip) ([GitHub](https://github.com/dlampart/Pascal))

Lamparter D, et al. (2016) Fast and Rigorous Computation of Gene and Pathway Scores from SNP-Based Summary Statistics. PLoS Comput Biol. 2016 Jan 25;12(1):e1004714. 
doi: 10.1371/journal.pcbi.1004714

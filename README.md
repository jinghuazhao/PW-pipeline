## PW-pipeline

PathWay pipeline using GWAS summary statistics

### Introduction

Pathway analysis is an important component of work in a GWAS. Broadly, a pathway analysis involves SNP annotation, such as Variant Effect Predictor (VEP), gene 
analysis such as VEGAS2, and gene set analysis. Visualisation of a particular region is typically facilitated with locuszoom, while network(s) built from pathway 
analysis can be visualised via gephi or Cytoscape, which accepts a collection of edges, directed or undirected to build a network. Aspects to consider include part 
or all databases, individual vs summary statistics, computing speed, with and without tissue enrichment. 

### Methods

This repository inovles several software for pathway analysis using GWAS summary statistics, namely,

* Meta-Analysis Gene-set Enrichment of variaNT Associations (MAGENTA)
* Generalized Gene-Set Analysis of GWAS Data (MAGMA)
* Data-Driven Expression Prioritized Integration for Complex Traits (DEPICT)
* Pathway scoring algorithm (PASCAL)

It using several databases that can be supplied to MAGENTA, MAGMA and PASCAL, and also in particular with respect to a DEPICT-to-PASCAL database from PASCAL 
developers which can be used across all software.

Their features are briefly described as follows.

For the discretised DEPICT-to-PASCAL database involving ENSEMBL IDs, the following code readily helps,
```
library(EnsDb.Hsapiens.v86)
chrall <- select(EnsDb.Hsapiens.v86, keys=paste(1:22), keytype="SEQNAME")
chrall_table <- subset(chr22[selcol],!duplicated(chr22[selcol]))
write.table(chrall_table,file="GS.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
```

**MAGENTA**

It first maps SNPs to genes taking 110 Kb upstream and 40 Kb downstream of each gene as extended boundaries to include regulatory regions. Each gene is then 
assigned a genetic set (GS) score, which is the P-value of the most significant SNP within the gene’s extended boundaries, corrected for six potential confounding 
factors of physical and genetic properties of genes through a step-wise multiple linear regression: (i) the physical size of the gene, (ii) number of SNPs per 
kilobase for each gene, (iii) estimated number of independent SNPs per gene, (iv) number of recombination hotspots spanning each gene, (v) genetic distance of the 
gene and (vi) linkage disequilibrium (LD) unit distance per gene.

**MAGMA**

The gene-set analysis is divided into two parts. In the first part a gene analysis is performed to quantify the degree of association each gene has with the 
phenotype. In addition the correlations between genes are estimated. These correlations reflect the LD between genes, and are needed in order to compensate for the 
dependencies between genes during the gene-set analysis. The gene p-values and gene correlation matrix are then used in the second part to perform the actual 
gene-set analysis.

**DEPICT**

DEPICT performs gene set enrichment analyses by testing whether genes in GWAS-associated loci are enriched for reconstituted versions of known molecular pathways 
(jointly referred to as reconstituted gene sets). The reconstitution is accomplished by identifying genes that are co-regulated with other genes in a given gene set 
based on a panel of 77,840 gene expression microarrays. Genes that are found to be transcriptionally co-regulated with genes from the original gene set are added to 
the gene set, which results in the reconstitution. Several types of gene sets were reconstituted in DEPICT: 5,984 protein molecular pathways derived from 169,810 
high-confidence experimentally derived protein-protein interactions, 2,473 phenotypic gene sets derived from 211,882 gene-phenotype pairs from the Mouse Genetics 
Initiative, 737 Reactome database pathways, 184 Kyoto Encyclopedia of Genes and Genomes (KEGG) database pathways and 5,083 Gene Ontology database terms. In total, 
14,461 gene sets were assessed for enrichment in genes in associated regions. DEPICT also facilitates tissue and cell type enrichment analyses by testing whether 
the genes in associated regions are highly expressed in any of the 209 MeSH annotations for 37,427 microarrays on the Affymetrix U133 Plus 2.0 Array platform.

Input to DEPICT was 6696 SNPs achieved 5x10^-8 genomewide significance, from which 111 regions were clumped using PLINK options --clump-kb 500 --clump-p1 5e-08 
--clump-r2 0.1. 209 genes were listed for tissue-specific Z-scores and there were 72 independent loci. 10968 with reconstituted gene set enrichment Z-scores and 205 
prioritised genes.

**PASCAL**

Gene scores are obtained by aggregating SNP p-values from a GWAS meta-analysis while correcting for LD using a reference population via the max and sum of 
chi-squared statistics based on the most significant SNP and the average association signal across the region, respectively. This part by default is done for 
msigBIOCARTA_KEGG_REACTOME.gmt (1,077 pathways) or msigdb.v4.0.entrez.gmt containing all MSigDB genes (10,295 pathways). As of 21 June, 7,949 and 12,198 genes have 
been done, respectively.

Gene sets are based on external databases for reported pathways by combining the scores of genes that belong to the same pathways. Pathway enrichment of 
high-scoring (potentially fused) genes is evaluated using parameter-free procedures (chi-square or empirical score), avoiding any p-value cut-off inherent to 
standard binary enrichment tests.

### Databases

**MAGENTA**

The six databases (\_db) contain a total of 10,327 entries were distributed with the MATLAB implementation: 

Name | Entries
-----|--------
GO_terms_BioProc_MolFunc | 9,433
Ingenuity_pathways | 92
KEGG_pathways | 168
PANTHER_BioProc | 241
PANTHER_MolFunc | 252
PANTHER_pathways | 141

Only 2,529 contain 10 or more genes were used by MAGENTA by default.

**MSigDB**

Gene database | N
--------------|---
c2.all.v6.0.entrez.gmt | 4,731 
msigBIOCARTA_KEGG_REACTOME.gmt | 1077 
msigdb.v4.0.entrez.gmt | 10295 
msigdb.v6.0.entrez.gmt | 17779 

**DEPICT***

In line with the default setup for MAGENTA, pathways with less than 10 genes were excluded leading to 2529 (0.05/2529=1.977066e-05) pathways. Except DEPICT, 
category 2 (c2) or all of pathways in Molecular Signatures Database (MSigDB) v6 is used. The MSigDB is divided into 8 major collections and several sub-collections 
on 17,779 gene sets, c2 containing 4,731 (0.05/4731=1.056859e-05) curated gene sets (from various sources such as online pathway databases, the biomedical 
literature, and knowledge of domain experts. MSigDB/BIOCARTA_KEGG_REACTOME came as default to PASCAL and MSigDB v4.0 is distributed with PASCAL. 
An entry in the MAGENTA pathway database contains a pathway ID, followed by a list of Entrez gene IDs. Although MSigDB has an additional column after the pathway ID 
indicating URLs of the pathway, it would be ignored by MAGMA for instance since these URLs do not match any Entrez gene IDs thus has no effect on the results. This 
feature facilitates comparison of software considerably. Comparative as well as individual results including figures are kept in two excel workbooks called mmp.xlsx 
and xlsx.xlsx, respectively.


### Acknowledgements

The work drives from comparison of their performance using our own GWAS data.

### Software and references

DEPICT (Pers TH et al. Nat Commun. 2015 Jan 19;6:5890. doi: 10.1038/ncomms6890)

MAGENTA (Segre AV, et al (2010). PLoS Genet. 2010 Aug 12;6(8). pii: e1001058. doi: 10.1371/journal.pgen.1001058)7

MAGMA (de Leeuw C, et al. PLoS Comput Biol. 2015 Apr; 11(4): e1004219. doi:  10.1371/journal.pcbi.1004219)

PASCAL (Lamparter D, et al. PLoS Comput Biol. 2016 Jan 25;12(1):e1004714. doi: 10.1371/journal.pcbi.1004714. eCollection 2016 Jan.)

## PW-pipeline

PathWay pipeline using GWAS summary statistics

### Introduction

This a repository for pathway analysis using GWAS summary statistics; the software involved are as follows,

* Meta-Analysis Gene-set Enrichment of variaNT Associations (MAGENTA)
* Generalized Gene-Set Analysis of GWAS Data (MAGMA)
* Data-Driven Expression Prioritized Integration for Complex Traits (DEPICT)
* Pathway scoring algorithm (PASCAL)

It using several databases that can be supplied to MAGENTA, MAGMA and PASCAL, and also in particular with respect to a DEPICT-to-PASCAL database from PASCAL 
developers which can be used across all software.

### Methods

Their features are briefly described as follows.

For the discretised DEPICT-to-PASCAL database involving ENSEMBL IDs, the following code readily helps,
```
library(EnsDb.Hsapiens.v86)
chrall <- select(EnsDb.Hsapiens.v86, keys=paste(1:22), keytype="SEQNAME")
chrall_table <- subset(chr22[selcol],!duplicated(chr22[selcol]))
write.table(chrall_table,file="GS.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
```

### Acknowledgements

The work drives from comparison of their performance using our own GWAS data.

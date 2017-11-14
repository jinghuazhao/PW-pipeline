cat GO_terms_BioProc_MolFunc_db Ingenuity_pathways_db KEGG_pathways_db PANTHER_BioProc_db PANTHER_MolFunc_db PANTHER_pathways_db | \
awk '{$1="MAGENTA"};1' FS="\t" OFS="\t" > magenta.db

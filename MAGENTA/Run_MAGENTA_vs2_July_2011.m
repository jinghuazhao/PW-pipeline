%	Copyright, February 2009, Ayellet Segre, Mark Daly, David Altshuler, Broad Institute, 7 Cambridge Center, Cambridge, MA 02142, USA
%
%	The MAGENTA (Meta-Analysis Gene-set Enrichment of VariaNT Associations) software was written by Ayellet Segre in the labs of David Altshuler and Mark Daly
%	Code written in Matlab version R2009b.
%
%	This code is part of the MAGENTA software package vs1.2, that tests for enrichment of multiple modest genetic effects on a given complex disease or trait,
%	in predefined sets of genes or loci.
%
%	This software accompanies the paper:
%	Ayellet V. Segre, DIAGRAM Consortium, MAGIC investigators, Leif Groop, Vamsi K. Mootha, Mark J. Daly, and David Altshuler. Common Inherited Variation in
%	Mitochondrial Genes is not Enriched for Associations with Type 2 Diabetes or Related Glycemic Traits. PLoS Genetics (in press). June 2010.
%
%	Disclaimer: This software is distributed as is. The authors take no responsibility for any use or misuse.
%
%	If your work benefits from the use of the MAGENTA software package please cite the reference above.
%
%	For questions or comments please contact Ayellet Segre at asegre@broadinstitute.org. You can check for updates at: http://www.broadinstitute.org/mpg/magenta
%
%	Last updated: June 23, 2010
%
%

function Run_MAGENTA_vs2_July_2011(GWAS_SNP_score_file_name,GeneSet_db_file_name,Flag_gene_set_file_name,GWAS_OR_filename,SNP_rs_num_file_name,exp_label,Flag_gs,Remove_gs,num_gs_simul,Gene_boundr_upstr,Gene_boundr_downstr,GSEA_method,top_percen_cutoffs,min_gs_size,max_gs_size,choose_unique_genes,print_trait_gene_names_in_gs,print_gene_scores,print_rs_num,print_best_SNP_OR)

% Sample run: function Run_MAGENTA_vs1_May10_2010('SNPChrNumPosZscorePval_LDL_MetaAnaly_20k_010810','PANTHER_lipid_only_BiolProcesses_GSEA_format_120609','LDL_GeneIDs_Meta20k_010810','LDL_Lipid_Pathways',1,0,10000,110000,40000,1,[5 25],10,50000,1,1,0,0,0);


%  BEFORE BEGINNING MAKE SURE LSF and MATLAB ARE RUNNING (e.g. 'use LSF', and 'use Matlab')

%%%%%  !!!!  PLEASE FILL OUT THE FOLLOWING INPUT FILES AND PARAMETERS:

%%%%%           VARIABLES TO BE DEFINED BY USER             %%%%%%%%%

% REQUIRED INPUT FILES:
%--------------------

% Input file #1:

% Name of file with GWAS results table in the following format (from individual GWAS or GWAS meta-analyses):
% REMOVE SUFFIX FROM FILE NAME
% Each row refers to a single SNP (do not add header)
% columns: (1) SNP chr number, (2) SNP chr position, bp, (3) GWAS z-scores (if available), (4) GWAS p-value (required)

%# GWAS_SNP_score_file_name ='DIAGRAMv3_061909_ChrNumPosZPval_hg19';

% Input file #2:

% Name of file with gene sets to be tested by GSEA (REMOVE SUFFIX FROM FILE NAME):
% File format: Each row refers to a different gene set (do not add header). All entries HAVE to be tab-delimited (very important, not space-delimited)
% Columns: (1) database name [name sould not be separated by spaces], (2) gene set name [name should not be separated by spaces or underscores],
% (columns 3 through n): Entrez IDs of all genes in gene set separated by tabs.

%# GeneSet_db_file_name ='PANTHER_lipid_only_BiolProcesses_GSEA_format_120609'; %'GO_PANTBP_PANTMF_PANT_ING_KEGG2010';
% contains gene sets from several databases: Gene ontology, PANTHER (biological processes, 
% molecular functions, metabolic and signaling pathways), KEGG, Ingenuity pathways
% It takes about 15-20 hours to run MAGENTA on all four databases inputed together; for faster running time each database can be run separately (See README for file names)

% Other gene set files you can use (described in README file) are:
% 'GO_PANTHER_INGENUITY_KEGG_REACTOME_BIOCARTA' 
% 'GO_terms_BioProc_MolFunc_db'
% 'KEGG_pathways_db'
% 'Ingenuity_pathways_db'
% 'PANTHER_BioProc_db'
% 'PANTHER_MolFunc_db'
% 'PANTHER_pathways_db'
% 'Mitochondrial_gene_sets'; % set of nuclear encoded mitochondrial genes from MitoCarta compendium


%OPTIONAL INPUT FILES:
%--------------------

% Input file #3 (recommended):

% File name for list of a set of Gene IDs or Gene Names of genes of interest, such as known disease/trait genes
% (REMOVE SUFFIX FROM FILE NAME)
% This list can be used (1) to label the presence of onr or more of these genes in each of the analyzed gene sets/pathways;
% (2) one can choose to remove these genes from the GSEA analysis.
% Fomart: a single column of gene Entrez IDs or gene symbols - preferably gene IDS (do not add header)
% User can use function: "Extract_genes_around_SNPs.m" to extract all geners within a given physical distance around a given set
% of SNPs. instructions for use is in the function's code/

Flag_gene_set_file_name =''; % 

% State whether to work with Gene ID or Gene Symbol
match_genes_GeneID = 1; % 1=match genes according to gene ID; 0=match genes according to gene HUGO/RefFlat symbols

% Input file #4:

% If available, name of file that contains the odds ratios of all SNPs in first column.
% Optional - to add upper and lower 95% confidence intervals in second and third columns.
% REMOVE HEADER FROM FILE.
% should be same number of rows as GWAS SNP p-value file (GWAS_SNP_score_file_name)
% (REMOVE SUFFIX FROM FILE NAME)

% Used only when printing out gene association scores for all genes in a limited number of gene sets 
% that are listed in the Input file #2.

GWAS_OR_filename =''; 

% Input file #5:

% Name of file that contains in one column the rs numbers of all SNPs without header
% should be same number of rows as GWAS SNP p-value file (GWAS_SNP_score_file_name)
% (REMOVE SUFFIX FROM FILE NAME)

% Used only when printing out gene association scores for all genes in a limited number of gene sets 
% that are listed in the Input file #2.

SNP_rs_num_file_name=''; 

%%%%%%%%%	ENTER INPUT PARAMTERS HERE:     %%%%%%%%%

% Enter here name of experiment (e.g. disease/trait name and gene set database)
exp_label = GeneSet_db_file_name;
code_name = 'Run_MAGENTA_vs2_March_2011.m';

% Please specify the genome build of the SNP positions of your inputed GWAS or meta-analysis SNP association p-values
Genome_build='NCBI37'; % Genome builds NCBI36 (hg18) and NCBI37 (hg19) are currently available

Flag_gs = 0; % 1=flag the input subset of genes (e.g. disease/trait genes)

% 1 = remove a predefined subset of genes listed in 'Flag_gene_set_file_name' from analysis (e.g. genes in validated disease association regions) from GSEA;
% 0 = do not remove predefined set of genes from GSEA
Remove_gs=0;

% 1 = remove genes that fall in HLA region (chr6:29,710,331-33,150,000)
% 0 = keep all genes
Remove_HLA=1; % (default: remove HLA region due to high LD and high gene density in region)

num_gs_simul = 10000; % number of gene set permutations to begin with for GSEA p-value estimation. number will be increased with GSEA p-value<1e-4.
% (i.e. number of times random gene sets of identical size to tested gene-set are randomly sampled from genome for GSEA significance estimation)
num_perm_limit=1000000; % upper limit of gene set permutations (number may increase per gene set if GSEA p<1/num_gs_simul, and will stop at the upper bound set here).
			% RECOMMENDED UPPER BOUND=1e6 PERMUTATIONS

Gene_boundr_upstr = 110000; % added distance (in base pair units) uptream to gene's most extreme start site
Gene_boundr_downstr = 40000; % added distance (in base pair units) downtream to gene's most extreme end site

% Choose "1" if you want to run GSEA, otherwise "0" 
run_GSEA=1; 
% Choose GSEA statistical test: 1=GSEA-cutoff test, 4=Mann-Whitney rank sum test (set "calculate_FDR" option to '0' since option not available yet for rank sum test)
GSEA_method = 1;

% cutoffs to use for GSEA test (top X% of all gene scores ranked in order of significance)
% Simulations suggest 95 percentile cutoff has optimatl power (of 99, 95, 90, 75 and 50 percentiles of all gene scores)
% However, one might want to chose 75 percentile if complex disease or trait is highly polygenic
top_percen_cutoffs = [5 25];  % default: 95 percentile and 75 percentile

% limit of gene set size
min_gs_size = getenv('min_gs_size'); % minimum number of genes in gene set (default>=10 genes)
max_gs_size = 2000; % maximum number of genes in a gene set (default=no limit)

choose_unique_genes = 1; % 1 = correct for physcial proximity of genes in each given gene set, by removing all but one gene of a subset of genes
% that were assigned the same best SNP; choose gene with most significant corrected score from each such subset of genes.
% 0 = do not apply physcial clustering correction.

% Enter 1 to print in the output results file gene names and Entrez IDs of those genes in the input list of genes in file: Flag_gene_set_file_name (e.g. trait/disease genes)
% that belong to each analyzed pathway or gene set.
print_trait_gene_names_in_gs = 0;

% Enter 1 for printing gene association scores for all gene sets present in Input file #2. Each gene set in a separate file 
% (ending with suffix .geneset; ~10-50 kbytes per file on average). 
% To print the gene scores for a limited number of gene sets, include only the gene sets of interst in Input file #2.
% Coming soon - the option to state the names of pathways for which gene score print outs is desired.
print_gene_scores=0; % 1=print list of genes and their scores for one or a few gene sets

% The following two options (print_rs_num, print_best_SNP_OR) are relevant only if printing the gene scores for a given set 
% of gene sets into output files:
% Enter 1 for printing best SNP rs numbers in gene score output file for a given gene set; enter 0 to not print rs numbers
print_rs_num=0;
print_best_SNP_OR=0;

% Enter 1 to calculate a false discovery rate; enter 0 not to calculate FDR
calculate_FDR=1;


%%%%%           END OF VARIABLES DEFINED BY USER        %%%%%%%%%


%%%%%%%%%   ANALYSES BEGIN HERE   %%%%%%%%%%%%%%
tic;
rand('twister',sum(100*clock));

if (Genome_build == 'NCBI36')

	load AllHumanGeneChrPosStrand_18434Genes_RefFlat_111909; % hg18
	load hotspot_boundaries_b36; % hg18
	load CEUPruned_SNP_ChrNumPos_hg18_030811; % hg18

        clear HumanGeneChrPos  	
        HumanGeneChrPos = AllHumanGeneChrPosStrand_18434Genes_RefFlat_111909;

	clear hotspot_boundaries
	hotspot_boundaries=hotspot_boundaries_b36; % chr num, hotspot region start (bp), hotspot region end (bp), hotspot width (kb)

	clear PrunedSNPsChrPos
	PrunedSNPsChrPos=CEUPruned_SNP_ChrNumPos_hg18_030811; % chr num, chr pos (bp)

        % (2) Read in list of all Human genes; column 1: Gene ID, column 2: Gene Symbol (RefFlat)
        clear AllRefFlatGeneID AllRefFlatGeneNames
	[AllRefFlatGeneID,AllRefFlatGeneNames]=textread('RefFlatGeneSymbolGeneID_18434Genes_111909','%n%s'); % hg18

elseif (Genome_build == 'NCBI37')

	load AllHumanGeneChrPosStrand_RefSeq_hg19_072111;
	load hotspot_boundaries_b37_hg19_072111;
	load CEU_HapMap_pruned_SNPs_ChrNumPos_hg19_072111;

        clear HumanGeneChrPos
        HumanGeneChrPos = AllHumanGeneChrPosStrand_RefSeq_hg19_072111;

        clear hotspot_boundaries
        hotspot_boundaries=hotspot_boundaries_b37_hg19_072111; % chr num, hotspot region start (bp), hotspot region end (bp), hotspot width (kb)

        clear PrunedSNPsChrPos
        PrunedSNPsChrPos=CEU_HapMap_pruned_SNPs_ChrNumPos_hg19_072111; % chr num, chr pos (bp)

        % (2) Read in list of all Human genes; column 1: Gene ID, column 2: Gene Symbol (RefFlat)
        clear AllRefFlatGeneID AllRefFlatGeneNames
        [AllRefFlatGeneID,AllRefFlatGeneNames]=textread('AllHumanGeneNames_RefSeq_hg19_072111','%n%s'); % hg19

end


calculate_Scores_confounders=1; % 1=recalcaulted Best-SNP-per-gene score and confounders, 0=do not calculate Best-SNP-per-gene score and confounders
GeneScoreMetric = 1; % 1=Best SNP per gene score; in the future will test other scoring metrics such as an average over a set of SNP in linkage equilibrbium
score_signif_direct = 1;        % 1=lower gene score values are more significant (e.g. p-values); 0=higher gene score values are more significant
                                % (e.g. normalized gene scores, z-scores)

if (Remove_HLA==1)

	st=29710331;
	en=33150000;

	clear a
	a=HumanGeneChrPos;
	find_HLA_reg=find(a(:,1)==6 & ((a(:,2)<st & a(:,3)>st & a(:,3)<=en) | (a(:,2)>=st & a(:,3)<=en) | (a(:,2)>=st & a(:,2)<en & a(:,3)>en) | (a(:,2)<=st & a(:,3)>=en) ) );

        find_not_HLA_reg=setdiff([1:length(a)],find_HLA_reg);
        AllGeneChrPos=a(find_not_HLA_reg,:);
        clear a

else
        AllGeneChrPos=HumanGeneChrPos;
end

% record today's date:
dash=abs('-');
date_string=date;
date_string_in_numbers=abs(date_string);
positions_of_dash=find(date_string_in_numbers==dash);
todays_date = [num2str(date_string((positions_of_dash(1)+1):(positions_of_dash(2)-1))), num2str(date_string(1:(positions_of_dash(1)-1))), '_' , num2str(date_string((positions_of_dash(2)+3):length(date_string)))];


%%%%%%%%%%%% ------  ORIGINAL  -----


% (0) PRINT INTO LOG FILE:

% make directory
Output_dir = [num2str(exp_label) , '_',  num2str(num_gs_simul), 'perm_', num2str(todays_date)];
mkdir_com = ['system(''mkdir ', num2str(Output_dir), ''');'];
eval(mkdir_com);

disp(['All output data is saved in the following directory: ', num2str(Output_dir), '\n']);

GSEA_log_file_name = ([num2str(Output_dir), '/MAGENTA_run_', num2str(exp_label), '_', num2str(num_gs_simul), 'perm_', num2str(todays_date) , '.log']);
FID_log=fopen(GSEA_log_file_name,'w');

fprintf(FID_log, ['MAGENTA written by Ayellet Segre, Altshuler and Daly Labs, Date: ', num2str(todays_date), '\n\n']);
fprintf(FID_log, ['Summary file for running MAGENTA (Gene Set Enrichment Analysis (GSEA) on Genome-wide association study (GWAS) variant results).\n']);
fprintf(FID_log, ['Program run: ', num2str(code_name), '\n']);
fprintf(FID_log, ['GWAS used: ', num2str(exp_label), '\n']);
fprintf(FID_log, ['Number of randomly sampled gene sets for GSEA-GWAS p-value estimation is: ', num2str(num_gs_simul), '.\n']);
fprintf(FID_log, ['Gene boundaries used for mapping SNPs onto genes are: ', num2str(Gene_boundr_upstr/1000), 'kb upstream to most extreme gene transcript start position, and ', num2str(Gene_boundr_downstr/1000)  'kb downstream to most extreme gene transcript end position, taking gene orientation into account.\n']);


% (1) Open GSEA-GWAS results file

if (run_GSEA==1)

GSEA_results_file_name = ([num2str(Output_dir), '/MAGENTA_pval_GeneSetEnrichAnalysis_', num2str(exp_label), '_', num2str(Gene_boundr_upstr/1000), 'kb_upstr_', num2str(Gene_boundr_downstr/1000), 'kb_downstr_', num2str(num_gs_simul), 'perm_', num2str(todays_date) , '.results']);
FID_results=fopen(GSEA_results_file_name,'w');

% Print header

    if (GSEA_method==1)

	if (calculate_FDR==1)
            fprintf(FID_results, ['DB\tGS\tORIG_GS_SIZE\tEFF_GS_SIZE\t#_GENES_ABS_LIST\t#_GENES_WO_SCORE\t#_GENES_REM_PROXIM\tMED_GENE_SIZE_KB\tMEAN_GENE_SIZE_KB\tNOMINAL_GSEA_PVAL_95PERC_CUTOFF\tFDR_95PERC_CUTOFF\tEXP_#_GENES_ABOVE_95PERC_CUTOFF\tOBS_#_GENES_ABOVE_95PERC_CUTOFF\tNOMINAL_GSEA_PVAL_75PERC_CUTOFF\tFDR_75PERC_CUTOFF\tEXP_#_GENES_ABOVE_75PERC_CUTOFF\tOBS_#_GENES_ABOVE_75PERC_CUTOFF\t#_GENES_FLAGGED\tFLAGGED_GENE_NAMES\n']);
	else
	    fprintf(FID_results, ['DB\tGS\tORIG_GS_SIZE\tEFF_GS_SIZE\t#_GENES_ABS_LIST\t#_GENES_WO_SCORE\t#_GENES_REM_PROXIM\tMED_GENE_SIZE_KB\tMEAN_GENE_SIZE_KB\tNOMINAL_GSEA_PVAL_95PERC_CUTOFF\tEXP_#_GENES_ABOVE_95PERC_CUTOFF\tOBS_#_GENES_ABOVE_95PERC_CUTOFF\tNOMINAL_GSEA_PVAL_75PERC_CUTOFF\tEXP_#_GENES_ABOVE_75PERC_CUTOFF\tOBS_#_GENES_ABOVE_75PERC_CUTOFF\t#_GENES_FLAGGED\tFLAGGED_GENE_NAMES\n']);
	end

        % GSEA-wilcoxon rank test
        elseif (GSEA_method==4)
    
                % Print header
                fprintf(FID_results, ['DB\tGS\tORIG_GS_SIZE\tEFF_GS_SIZE\t#_GENES_ABS_LIST\t#_GENES_WO_SCORE\t#_GENES_REM_PROXIM\tMED_GENE_SIZE_KB\tMEAN_GENE_SIZE_KB\tRANKSUM_PVAL_HIGHER_TAIL\t#_GENES_FLAGGED\tFLAGGED_GENE_NAMES\n']);
	end

end % if run GSEA


clear AllGene_Names
if (match_genes_GeneID==1)
        AllGene_Names=AllRefFlatGeneID;
else
        AllGene_Names=AllRefFlatGeneNames;
end

if (Remove_HLA==1)
	
        AllGene_Names_woHLA = AllGene_Names(find_not_HLA_reg);
	clear AllGene_Names
	AllGene_Names=AllGene_Names_woHLA;

	AllRefFlatGeneID_woHLA=AllRefFlatGeneID(find_not_HLA_reg);
	clear AllRefFlatGeneID;
	AllRefFlatGeneID=AllRefFlatGeneID_woHLA;
	clear AllRefFlatGeneID_woHLA;

	AllRefFlatGeneNames_woHLA=AllRefFlatGeneNames(find_not_HLA_reg);
	clear AllRefFlatGeneNames;
	AllRefFlatGeneNames=AllRefFlatGeneNames_woHLA;
	clear AllRefFlatGeneNames_woHLA;
end
 
% (3.1) Read in GWA SNP z-scores and p-values after genomic control

Load_SNP_score_com = ['load ', num2str(GWAS_SNP_score_file_name)];
eval(Load_SNP_score_com);

% 	!!! NEEDS TO BE ADDED: TEST IF FILE NAME HAS SUFFIX. IF SO REMOVE SUFFIX FROM FILE NAME.

% columns: (1) chr num, (2) chr pos, bp, (3) GWAS z-scores, (4) GWAS p-value 
Read_SNP_score_com = ['GWAS_SNP_ChrNumPosZscoresPval = ', num2str(GWAS_SNP_score_file_name)  ,  ';'];
eval(Read_SNP_score_com);

[n_r, n_c]=size(GWAS_SNP_ChrNumPosZscoresPval);

if (n_c<4 & max(GWAS_SNP_ChrNumPosZscoresPval(:,3))<=1)
Input_file_SNPpval_only=1; % =1 if only SNP association p-values are given; =0 if SNP association z-scores are given with or without p-values
else 
Input_file_SNPpval_only=0;
end


% ODDS ratio with lower and upper bound 95% confidence intervals
if (print_best_SNP_OR==1)
 if (GWAS_OR_filename)
  Load_SNP_OR_com = ['load ', num2str(GWAS_OR_filename)];
  eval(Load_SNP_OR_com);

  Read_SNP_OR_com = ['GWAS_SNP_OR_L95_U95 = ', num2str(GWAS_OR_filename)  ,  ';'];
  eval(Read_SNP_OR_com);
 end
end

% Convert GWA SNP p-values to z-scores 

if (Input_file_SNPpval_only==1)
	[GWAS_SNPChrNumPos_ZscoresCalculatedFromPval]=Converting_GWAS_TwoTailedPval_to_Zscores_082809(GWAS_SNP_ChrNumPosZscoresPval);
	clear GWAS_SNP_ChrNumPosZscoresPval
	GWAS_SNP_ChrNumPosZscoresPval=GWAS_SNPChrNumPos_ZscoresCalculatedFromPval;
end

% (3.2) Read in all SNP rs numbers in same order as SNP positions and z-scores/p-values in input table
% one column of all SNP rs numbers

if (SNP_rs_num_file_name)

 if (print_rs_num==1)

    	AllSNP_rs_num=cell([length(GWAS_SNP_ChrNumPosZscoresPval),1]);

	if (SNP_rs_num_file_name)
        	Read_SNP_rs_num_com = ['[AllSNP_rs_num]=textread(''', num2str(SNP_rs_num_file_name)  ,  ''',''%s'');'];
        	eval(Read_SNP_rs_num_com);
	end
 else
        AllSNP_rs_num=cell([length(GWAS_SNP_ChrNumPosZscoresPval),1]);
 end

else
        AllSNP_rs_num=cell([length(GWAS_SNP_ChrNumPosZscoresPval),1]);
end

% (4) Calculate 'gene association score' from local SNP association z-scores

% Best SNP per gene scoring metric

StrandOrient = AllGeneChrPos(:,6);
 
if (GeneScoreMetric==1)

 % AllHumanGeneChrPosStrand_RefFlat_020909 % Columns: (1) chr num, (2) Tx start, (3) Tx end, (4) CDS start, (5) CDS end, (6) strand orientation (1/0)
 % Best SNP per gene test statistic contains original sign
 best_pval_or_z=1; % use p-values to find best SNP per gene; if best_pval_or_z=0 use z-scores to find best SNP per gene
 [BestSNPPerGeneChrNumPosZScores,NumSNPs_per_gene,Best_SNP_rs]=ExtractGeneScoreBestSNP_PvalZscore_NumSNPsPerGene_092909(AllGeneChrPos(:,1:3),GWAS_SNP_ChrNumPosZscoresPval,Gene_boundr_upstr,Gene_boundr_downstr,StrandOrient,best_pval_or_z,AllSNP_rs_num);

 % BestSNPPerGeneChrNumPosZScores: columns: (1) Best SNP per gene Chr num, (2) Best SNP per gene Chr pos, bp, (3) Best SNP per gene GWA z-score
 Uncorr_score = abs(BestSNPPerGeneChrNumPosZScores); % Original GWAS score of best SNP per gene before correction for confounders

% if you want to save Best SNP per gene scores
%  save_com1 = ['save ', num2str(Output_dir) ,'/BestSNPPerGeneScores_uncorr_', num2str(Gene_boundr_upstr/1000), 'kb_upstr_', num2str(Gene_boundr_downstr/1000), 'kb_downstr_', num2str(exp_label), '_', num2str(todays_date), ' BestSNPPerGeneChrNumPosZScores -ascii -tabs'];
% eval(save_com1);

% save_com2 = ['save ', num2str(Output_dir) ,'/NumSNPs_per_gene_', num2str(Gene_boundr_upstr/1000), 'kb_upstr_', num2str(Gene_boundr_downstr/1000), 'kb_downstr_', num2str(exp_label), '_', num2str(todays_date), ' NumSNPs_per_gene -ascii -tabs'];
% eval(save_com2);


 % (5) Correct gene score according to confounders using step-wise multivariate linear regression analysis

 % define confounder matrix (rows: genes, columns: confounders)
 clear confounders interval
 interval = Gene_boundr_upstr+Gene_boundr_downstr; % total distance added to gene upstream and downstream, in base pairs

 % Size of most extreme transcript boundaries for each gene
 AllGeneSizes = AllGeneChrPos(:,3)-AllGeneChrPos(:,2);

 % confounder variables in per kb units (by dividing the variables by gene size plus added boundaries)
 GeneSize_plus_interval_in_kb = (AllGeneSizes + interval)/1000;

% save gene size with added boundaries in kb
% save_com = ['save ', num2str(Output_dir) ,'/GeneSize_plus_interval_in_kb_', num2str(Gene_boundr_upstr/1000), 'kb_upstr_', num2str(Gene_boundr_downstr/1000), 'kb_downstr_', num2str(exp_label), '_', num2str(todays_date), ' GeneSize_plus_interval_in_kb -ascii -tabs'];

% Calculate estimated number of independent SNPs that are in linkage equilibrium per gene region
tic;
[Num_Indep_SNPs_per_gene]=Calculate_NumIndepSNPsPerGene_030911(AllGeneChrPos,PrunedSNPsChrPos,Gene_boundr_upstr,Gene_boundr_downstr);

clear save_com
% save_com = ['save ', num2str(Output_dir) ,'/NumIndepSNPsperGene_', num2str(Gene_boundr_upstr/1000), 'kb_upstr_', num2str(Gene_boundr_downstr/1000), 'kb_downstr_', num2str(exp_label), '_', num2str(todays_date), ' Num_Indep_SNPs_per_gene -ascii -tabs'];
% eval(save_com);

t1=toc;
 
fprintf(FID_log, ['It took ', num2str(t1), ' seconds or ', num2str(t1/60) , ' minutes to calculate number of independent SNPs per gene.\n']);
disp(['It took ', num2str(t1), ' seconds or ', num2str(t1/60) , ' minutes to calculate number of independent SNPs per gene.\n']);


 tic; 
 % Calculate number of hotspots per gene 
 [NumHotspots_per_gene]=Calculate_NumHotspots_per_gene_030911(AllGeneChrPos(:,1:3),hotspot_boundaries,Gene_boundr_upstr,Gene_boundr_downstr,StrandOrient);
  t1=toc;
 disp(['It took ', num2str(t1), ' seconds or ', num2str(t1/60) , ' minutes to calculate number of hotspots per gene.\n']);
 fprintf(FID_log, ['It took ', num2str(t1), ' seconds or ', num2str(t1/60) , ' minutes to calculate number of hotspots per gene.\n']);
 

% four confounders (physical and linkage-related gene properties)
confounders = [GeneSize_plus_interval_in_kb NumSNPs_per_gene./GeneSize_plus_interval_in_kb Num_Indep_SNPs_per_gene./GeneSize_plus_interval_in_kb NumHotspots_per_gene./GeneSize_plus_interval_in_kb];


 tic;
 clear RegressCorrGeneScores_pval
 [RegressCorrGeneScores_pval, Residuals]=GeneScoreRegressCorr_122108(confounders,Uncorr_score(:,3));
 Corr_score = RegressCorrGeneScores_pval;
 t1=toc;
 disp(['It took ', num2str(t1), ' seconds or ', num2str(t1/60) , ' minutes to correct gene scores with stepwise regression.\n']);
 fprintf(FID_log, ['It took ', num2str(t1), ' seconds or ', num2str(t1/60) , ' minutes to correct gene scores with stepwise regression.\n']);

clear save_com   
save_com = ['save ', num2str(Output_dir) ,'/RegressCorrGeneScores_pval_', num2str(Gene_boundr_upstr/1000), 'kb_upstr_', num2str(Gene_boundr_downstr/1000), 'kb_downstr_', num2str(exp_label), '_', num2str(todays_date), ' RegressCorrGeneScores_pval -ascii -tabs'];
eval(save_com);
    
end % for Best-SNP-per-gene scoring metric

% Set SNP score - to be added in the future:
if (GeneScoreMetric==2)
SetScore = [];
Corr_score = SetScore(find_corr_score_isan);
end


% (6) Remove genes without gene scores (NaN), and if requested - remove input subset of genes from gene list and from analyses (e.g. known disease genes)
% Adjust all relevant gene matrices (e.g. Corrected and uncorrected gene scores, Gene chromosome positions, Gene IDs/Names accordingly

clear find_corr_score_isan find_score_isan_not_flag_gs

% indeces only of genes that have a corrected association score
find_corr_score_isan = find(~isnan(Corr_score) & ~isnan(Uncorr_score(:,3)));
AllGene_Names_score_isan = AllGene_Names(find_corr_score_isan);

if (Flag_gene_set_file_name)
if (Remove_gs==1)

         if (match_genes_GeneID==1)     % use Gene IDs

                Load_flag_gs_com = ['load ', num2str(Flag_gene_set_file_name), ';'];
                eval(Load_flag_gs_com);
        
                Assign_flag_gs_com = ['Flag_gs_GeneNames= ', num2str(Flag_gene_set_file_name)  ,  ';'];
                eval(Assign_flag_gs_com);

	else  % use Gene symbols

                Read_flag_gs_names_com =['[Flag_gs_GeneNames]=textread(''', num2str(Flag_gene_set_file_name) , ''',''%s'');'];
                eval(Read_flag_gs_names_com);

	end

	% !!! NOTE: The function 'intersect' and 'unique' order the output vector of numbers or strings in order of alphabet or in numerical order 

	Flag_gs_GeneNames=unique(Flag_gs_GeneNames);

        [AllGenes_wo_flag_gs, AllGenes_wo_flag_gs_ind]=setdiff(AllGene_Names,Flag_gs_GeneNames);
                
        % indeces of all genes that got assigned a corrected score and that do not overlap the input gene subset to remove
           clear find_score_isan_not_flag_gs_unsorted k 
	    find_score_isan_not_flag_gs_unsorted = intersect(AllGenes_wo_flag_gs_ind, find_corr_score_isan);
	    [find_score_isan_not_flag_gs, k]=sort(find_score_isan_not_flag_gs_unsorted);  % indeces of 'AllGene_Names' tat refers to genes that have a score
                                                                                       % and that do not include the flagged genes

else
        find_score_isan_not_flag_gs=find_corr_score_isan;
end 
else
        find_score_isan_not_flag_gs=find_corr_score_isan;
end

% all score and gene matrices adjusted to only genes with a score and
% without genes in flagged gene subset if requested to remove
 Corr_score_isan = Corr_score(find_score_isan_not_flag_gs);
 Uncorr_score_isan = abs(BestSNPPerGeneChrNumPosZScores(find_score_isan_not_flag_gs,:));  % Original GWAS score of best SNP per gene before correction

 AllGene_Names_isan = AllGene_Names(find_score_isan_not_flag_gs); 
 NumSNPs_per_gene_isan = NumSNPs_per_gene(find_score_isan_not_flag_gs);
 AllGeneSizes_isan = AllGeneSizes(find_score_isan_not_flag_gs);

% (7) Read in gene sets

ReadGSFile_com = ['[gs_GeneIDs_cellarray num_genes_gs]=multi_list2cell_GeneSetDB(''' ,  num2str(GeneSet_db_file_name), ''');'];
eval(ReadGSFile_com);

find_gs_min_num_genes = find(num_genes_gs>=min_gs_size & num_genes_gs<=max_gs_size);

% number of gene sets with a minimum number of genes (default=10 genes)
num_gene_sets_min_num_genes =length(find_gs_min_num_genes);

% cell array of gene sets and vector of Gene IDs only for gene sets with minimum number of genes
gs_GeneIDs_cellarray_min_gene_num = {gs_GeneIDs_cellarray{find_gs_min_num_genes}};

gs_GeneIDs_cellarray_min_gene_num{num_gene_sets_min_num_genes+1}{1}='TAIL'; % add an additional element for looping purposes later on

% number of genes in gene sets with minimum number of genes
num_genes_in_gs_min_gene_num = num_genes_gs(find_gs_min_num_genes);

%unique_resources = unique(gs_GeneIDs_cellarray_min_gene_num{:}{1}); % list of resources from which gene sets were compiled
unique_resources = gs_GeneIDs_cellarray_min_gene_num{1}{1};

% (8) PRINT INTO LOG FILE:

fprintf(FID_log, ['The input file name of gene sets analyzed is: ', num2str(GeneSet_db_file_name), '\n']);
fprintf(FID_log, ['Number of gene sets with at least ', num2str(min_gs_size), ' genes is ', num2str(num_gene_sets_min_num_genes), '.\n']);
fprintf(FID_log, ['Gene sets were taken from the following resources:','\n']);
fprintf(FID_log,'%s\t',[num2str(unique_resources)]);
if (print_gene_scores == 1)
	fprintf(FID_log, ['Gene scores for the inputed set of gene-sets were printed into files with the suffix .geneset.' , '\n']);
end

%for i=1:length(unique_resources)
%        fprintf(FID_log,'%s\t',[num2str(unique_resources{i})]);
%end

fprintf(FID_log,'\n'); 

if (run_GSEA==1)

fprintf(FID_log, ['The GSEA results file and other output files are saved in the directory: ', num2str(Output_dir), '\n']);
fprintf(FID_log, ['The GSEA results file name is: ', num2str(GSEA_results_file_name), '\n']);
fprintf(FID_log, ['Here is a key of the header subtitles for each columns: ', '\n']);
fprintf(FID_log, ['DB=databse.', '\n']);
fprintf(FID_log, ['GS=gene set.', '\n']);
fprintf(FID_log, ['ORIG_GS_SIZE=Original number of genes per gene set.', '\n']);
fprintf(FID_log, ['EFF_GS_SIZE=Effective number of genes per gene set analyzed by GSEA, after removing genes that were not assigned a gene score (e.g. no SNPs in their region), or', '\n']);
fprintf(FID_log, ['after adjusting for physcial clustering of genes in a given gene set (removing all but one gene from a subset of genes assigned the same best SNP, ', '\n']);
fprintf(FID_log, ['keeping the gene with the most significant gene score).', '\n']);
fprintf(FID_log, ['#_GENES_ABS_LIST = number of genes absent from the full human genes list used for the analysis.', '\n']);
fprintf(FID_log, ['#_GENES_WO_SCORE = number of genes that did not get assigned a score since they had no SNPs with association statistics in their extended gene region boundaries.', '\n']);
fprintf(FID_log, ['#_GENES_REM_PROXIM = Number of genes removed due to physical proximity adjustment (See above).', '\n']);
fprintf(FID_log, ['MED_GENE_SIZE_KB = median gene size of all genes in gene set in kilobase units.', '\n']);
fprintf(FID_log, ['MEAN_GENE_SIZE_KB = mean gene size of all genes in gene set in kilobase units.', '\n']);
fprintf(FID_log, ['NOMINAL_GSEA_PVAL_95PERC_CUTOFF = GSEA p-value using 95 percentile of all gene scores for the enrichment cutoff.', '\n']);
fprintf(FID_log, ['FDR_95PERC_CUTOFF = Estimated false discovery rate (q-value) using 95 percentile cutoff.', '\n']);
fprintf(FID_log, ['EXP_#_GENES_ABOVE_95PERC_CUTOFF = Expected number of genes with a corrected gene p-value above the 95 percentile enrichment cutoff.', '\n']);
fprintf(FID_log, ['OBS_#_GENES_ABOVE_95PERC_CUTOFF = Observed number of	genes with a corrected gene p-value above the 95 percentile enrichment cutoff.', '\n']);
fprintf(FID_log, ['NOMINAL_GSEA_PVAL_75PERC_CUTOFF = GSEA p-value using 75 percentile of all gene scores for the enrichment cutoff.', '\n']);
fprintf(FID_log, ['FDR_75PERC_CUTOFF =	Estimated false discovery rate (q-value) using 75 percentile cutoff.', '\n']);
fprintf(FID_log, ['EXP_#_GENES_ABOVE_75PERC_CUTOFF = Expected number of genes with a corrected gene p-value above the 75 percentile enrichment cutoff.', '\n']);
fprintf(FID_log, ['OBS_#_GENES_ABOVE_75PERC_CUTOFF = Observed number of genes with a corrected gene p-value above the 75 percentile enrichment cutoff.', '\n']);
fprintf(FID_log, ['#_GENES_FLAGGED = Number of genes in gene set that belong to the subset of genes flagged by the user (e.g. genes near validated disease/trait SNPs).', '\n']);
fprintf(FID_log, ['FLAGGED_GENE_NAMES = Gene symbols of genes in list of flagged genes list that belong to each gene set.', '\n']);

else

fprintf(FID_log, ['GSEA was not run on the input gene sets, as requested by the user.', '\n']);

end

% (9) LOOP OVER ALL GENE SETS

EnrichPval_multipleCutoffs_all_genesets=[]; % matrix that stores GSEA-GWAS p-values for all gene sets

Presence_flag_gs_genes = zeros(1,num_gene_sets_min_num_genes);  % 1=at least one gene of gene subset is present in given gene set; 0=no genes in 
								% gene subset is present in given gene set


Eff_num_genes_gs = zeros(1,num_gene_sets_min_num_genes); 

Record_results_all_genesets = zeros(num_gene_sets_min_num_genes,(length(top_percen_cutoffs)+10));

Rand_gs_fract_p95 = zeros(num_gs_simul,num_gene_sets_min_num_genes); % matrix that records fraction of genes above the 95 percentile cutoff 
								     % for each simulated gene set (row) for each gene set (column)
Rand_gs_fract_p75 = zeros(num_gs_simul,num_gene_sets_min_num_genes); % matrix that records fraction of genes above the 75 percentile cutoff 
                                                                     % for each simulated gene set (row) for each gene set (column)


first_db = gs_GeneIDs_cellarray_min_gene_num{1}{1};
gs_ind = 0;

clear db_name GeneSetLabel Original_gs_size Eff_num_genes_gs Med_genesize_gs Mean_genesize_gs Num_LargeGenes_gs Num_SmallGenes_gs
clear Num_genes_no_match_list num_genes_gs_wo_score num_genes_removed_overlap GeneSetLabel_array GeneSetGeneNames_cell_of_arrays
clear Rand_gs_fract_p95 Rand_gs_fract_p75 obs_fract_above_cutoff_p95 obs_fract_above_cutoff_p75 obs_num_gs_above_cutoff_p95 obs_num_gs_above_cutoff_p75
clear GeneSetGeneNames_cell_of_arrays

for gs=1:num_gene_sets_min_num_genes  % for each gene set with minimum number of genes in gene set
        
	gs_ind=gs_ind+1;
	db_name{gs_ind} = gs_GeneIDs_cellarray_min_gene_num{gs}{1};

	clear GeneSetLabel
	GeneSetLabel = gs_GeneIDs_cellarray_min_gene_num{gs}{2};
        GeneSetLabel_array{gs_ind}=GeneSetLabel;

	clear GeneSetGeneNames
	GeneSetGeneNames = gs_GeneIDs_cellarray_min_gene_num{gs}{3}; % this is a vector of Gene IDs or names		

	%%%%% run GSEA only on selected set of gene sets
	%[Select_DB_Names,Select_GS_Names]=textread('Lipid_Selected_genesets_012510','%s%s');

	%clear ChosenGeneSet
	%[ChosenGeneSet,v11,v22] = intersect([Select_GS_Names],[GeneSetLabel]);

	%if (length(ChosenGeneSet)>0)
	%%%%%

	% (9.1) Make sure gene set ID/name vector is unique

	GeneSetGeneNames=unique(GeneSetGeneNames);
        GeneSetGeneNames_cell_of_arrays{gs_ind}=GeneSetGeneNames; 

	Original_gs_size(gs_ind) = length(GeneSetGeneNames);  % Original number of genes in gene set


	% (9.2) Find genes that have a match with our whole human genes list

        clear GS_GeneNames_match w1 w2
        [GS_GeneNames_match,w1,w2]=intersect(AllGene_Names,GeneSetGeneNames);
	
	% number of genes that do not have a match on our all human genes list
	Num_genes_no_match_list(gs_ind) = Original_gs_size(gs_ind)-length(GS_GeneNames_match);

	% (9.3) Remove genes without a corrected score and without the flagged gene subset if requested
        clear GS_GeneNames_match_isan find_gs_genes_isan_sorted
        [GS_GeneNames_match_isan, find_gs_genes_isan_sorted]=intersect(AllGene_Names_isan,GS_GeneNames_match); 

	% number of genes in gene set without a corrected score
	clear GS_GeneNames_match_w_score find_gs_genes_w_score_sorted 
        [GS_GeneNames_match_w_score, find_gs_genes_w_score_sorted]=intersect(AllGene_Names_score_isan,GS_GeneNames_match);	
	num_genes_gs_wo_score(gs_ind) = length(GS_GeneNames_match)-length(GS_GeneNames_match_w_score);

	[find_gs_genes_isan, sort_find_gs_genes_isan_ind] = sort(find_gs_genes_isan_sorted);

	% uncorrected score for gene set after removing genes without score (columns: (1) Best SNP per gene chr num, (2) Best SNP per gene chr pos, (3) Best SNP per gene z-scores
	clear Uncorr_score_gs_isan Corr_score_gs_isan
	Uncorr_score_gs_isan = Uncorr_score_isan(find_gs_genes_isan,:);
	Corr_score_gs_isan = Corr_score_isan(find_gs_genes_isan,:);
	GS_GeneNames_match_isan = AllGene_Names_isan(find_gs_genes_isan,:);	% SHOULD BE EQUIVALENT TO: GS_GeneNames_match_isan = GS_GeneNames_match_isan(sort_find_gs_genes_isan_ind);


	% (9.4) For subsets of genes that share the same best-SNP-per-gene score, choose the gene with the best score, remove the rest

	if (choose_unique_genes==1)

		clear Uncorr_score_geneset_pval_sorted_ind a
	        [a, Uncorr_score_geneset_pval_sorted_ind]= sort(abs(Corr_score_gs_isan),'ascend'); % sort SNPs according to corrected p-values
            	
	        Uncorr_score_SNP_Chrpos_geneset_sorted = Uncorr_score_gs_isan(Uncorr_score_geneset_pval_sorted_ind,1:2); % Chr num and pos of best SNP per gene scores in order of SNP score significance

    		clear Unique_genes_ind a j 
		[a,Unique_genes_ind,j]=unique(Uncorr_score_SNP_Chrpos_geneset_sorted,'rows','first');

	    	num_genes_removed_overlap(gs_ind) = length(GS_GeneNames_match_isan)-length(Unique_genes_ind); % # genes removed due to best SNP per gene overlap
		Eff_num_genes_gs(gs_ind) = length(Unique_genes_ind); %  Final number of genes in gene set taken forward to GSEA-GWAS analysis

		% Retrieve gene set Gene names/IDs after removing genes due to signal overlap
		clear s1 GS_GeneNames_match_isan_no_overlap
        	s1=GS_GeneNames_match_isan(Uncorr_score_geneset_pval_sorted_ind);
 		GS_GeneNames_match_isan_no_overlap=s1(Unique_genes_ind);

	        clear v find_gene_set_ind_sorted k
        	[v, find_gene_set_ind_sorted, k]=intersect(AllGene_Names_isan,GS_GeneNames_match_isan_no_overlap);

		clear find_gene_set_ind v
		[find_gene_set_ind, v] = sort(find_gene_set_ind_sorted);

		clear GS_GeneNames_final
		GS_GeneNames_final = GS_GeneNames_match_isan_no_overlap;

    	    else
        
        	num_genes_removed_overlap = 0; % # genes removed due to best SNP per gene overlap
		Eff_num_genes_gs(gs_ind) = length(find_gs_genes_isan); %  Final number of genes in gene set taken forward to GSEA-GWAS analysis
  
		clear GS_GeneNames_final
		GS_GeneNames_final = GS_GeneNames_match_isan;
		find_gene_set_ind = find_gs_genes_isan;
    	    end

    	% record gene size of all genes in gene set
        clear GS_genesize_vect
    	GS_genesize_vect = AllGeneSizes_isan(find_gene_set_ind,:); 
    
  	Med_genesize_gs(gs_ind) = median(GS_genesize_vect(find(~isnan(GS_genesize_vect))))/1000; % median gene size of all genes in gene set in kb units        
	Mean_genesize_gs(gs_ind) = mean(GS_genesize_vect(find(~isnan(GS_genesize_vect))))/1000; % mean gene size of all genes in gene set in kb units
        
    	Num_LargeGenes_gs(gs_ind) = length(find(GS_genesize_vect>=100000));  % Number of genes in gene set with size >=100 kb
    	Num_SmallGenes_gs(gs_ind) = length(find(GS_genesize_vect<=10000));  % Number of genes in gene set with size <=10 kb

            
	% (9.5) Record presence of genes from input subset of genes in given gene set (e.g. known disease genes)

if (Flag_gene_set_file_name)

   if (Flag_gs == 1) 

     	 % Flag_gs_GeneNames: list of gene names/IDs
         if (match_genes_GeneID==1)	% use Gene IDs
		Load_flag_gs_com = ['load ', num2str(Flag_gene_set_file_name), ';'];
		eval(Load_flag_gs_com);        

		Assign_flag_gs_com = ['Flag_gs_GeneNames= ', num2str(Flag_gene_set_file_name)  ,  ';'];
		eval(Assign_flag_gs_com);

         else	% use Gene symbols

		Read_flag_gs_names_com =['[Flag_gs_GeneNames]=textread(''', num2str(Flag_gene_set_file_name) , ''',''%s'');'];         
		eval(Read_flag_gs_names_com);

         end

                clear q w x
                [q,w,x]=intersect(GS_GeneNames_match,unique(Flag_gs_GeneNames));

                if (length(q)>0) 
                    Presence_flag_gs_genes(gs_ind)=length(q); 
                else
                    Presence_flag_gs_genes(gs_ind)=0;                     
                end
    else
                    Presence_flag_gs_genes(gs_ind)=0;
    end  % if flag a subset of genes
else
                    Presence_flag_gs_genes(gs_ind)=0;
end


	% (9.6) Print gene scores for each gene set

	% Have option for a small number of gene sets to print the list of genes with their scores and the rs number and chr position of the
	% best SNP per gene

	if (print_gene_scores == 1)

                output_genescore_file_name=([num2str(Output_dir), '/GeneAssocScores_', num2str(exp_label), '_',  num2str(GeneSetLabel) , '_' ,  num2str(Gene_boundr_upstr/1000), 'kb_upstr_' , num2str(Gene_boundr_downstr/1000) , 'kb_downstr_' , num2str(todays_date) , '.geneset']);
		FID_geneset=fopen(output_genescore_file_name,'w');

		fprintf(FID_geneset,['Database\tGene_Set\tGene_Symbol\tEntrez_ID\tGene_p-value\tGene_Chr_Num\tGene_Start_Pos\tGene_End_Pos\tGene_Size_kb\tNum_SNPs_per_Gene\tNum_Indep_HapMap_SNPs_per_Gene\tNum_RecombHotspots_per_Gene\t']);

		if (SNP_rs_num_file_name) 
			if (print_rs_num==1) 	fprintf(FID_geneset,['Best_SNP_rs\t']);	end 
		end

		fprintf(FID_geneset,['Best_SNP_Chr_Num\tBest_SNP_Chr_Pos\tBest_SNP_Z\tBest_SNP_pval\t']);

		if (print_best_SNP_OR==1) 

	                [num_row_OR,num_col_OR]=size(GWAS_SNP_OR_L95_U95);
			if (num_col_OR==1)
				fprintf(FID_geneset,['Best_SNP_ODDS_RATIO_OR_EFFECT\t']);
			elseif (num_col_OR>1)
				fprintf(FID_geneset,['Best_SNP_ODDS_RATIO_OR_EFFECT\tBest_SNP_L95_CI\tBest_SNP_U95_CI\t']);	
			end
		end

		fprintf(FID_geneset,['1=Flagged_gene\n']);

		for i=1:length(GeneSetGeneNames)

			clear find_gene_a
			[MatchedGene, find_gene_a, find_gene_b]=intersect(AllGene_Names,GeneSetGeneNames(i));

			if (length(find_gene_a)==1)  % if given gene exists on our all human gene list

				fprintf(FID_geneset, '%s\t', [num2str(db_name{gs_ind})]);           % Resource name
			        fprintf(FID_geneset, '%s\t', [num2str(GeneSetLabel_array{gs_ind})]);          % Gene set name
			
				fprintf(FID_geneset, '%s\t', num2str(AllRefFlatGeneNames{find_gene_a}) );       % Gene name
				fprintf(FID_geneset, '%1.0f\t', AllRefFlatGeneID(find_gene_a) );		% Gene ID

				fprintf(FID_geneset, '%E\t', Corr_score(find_gene_a));                     % Regression-corrected gene score

				fprintf(FID_geneset, '%1.0f\t%1.0f\t%1.0f\t', AllGeneChrPos(find_gene_a,1:3));		% Gene chr position, bp
				fprintf(FID_geneset, '%1.0f\t', AllGeneSizes(find_gene_a)/1000);          		% Gene size, kb
				fprintf(FID_geneset, '%1.0f\t', NumSNPs_per_gene(find_gene_a));				% Number of SNPs per gene    
                                fprintf(FID_geneset, '%1.0f\t', Num_Indep_SNPs_per_gene(find_gene_a));                  % Estimated number of independent (pruned) SNPs per gene based on HapMap CEU
				fprintf(FID_geneset, '%1.0f\t', NumHotspots_per_gene(find_gene_a)); 			% Number recombination hotspots (from Myers et al, Science 2005) 
															% spanning each gene region (e.g. +110kb, -40kb)

                                % Best SNP per gene rs #
				if (SNP_rs_num_file_name)
					if (print_rs_num==1 & calculate_Scores_confounders==1)
                	                	fprintf(FID_geneset, '%s\t', num2str(Best_SNP_rs{find_gene_a}) );         % Best SNP per gene rs number  
					end
				end

				fprintf(FID_geneset, '%1.0f\t%1.0f\t%5.4f\t%E\t', Uncorr_score(find_gene_a,1:4) );		% Best SNP per gene chr num, chr pos, z-score, p-value

                                if (print_best_SNP_OR==1)
                                        [num_row_OR,num_col_OR]=size(GWAS_SNP_OR_L95_U95);

                                        if (~isnan(Best_SNP_rs{find_gene_a}))
                                                clear find_best_SNP_rs
                                                find_best_SNP_rs=strmatch(Best_SNP_rs{find_gene_a},AllSNP_rs_num,'exact');

                                                for or=1:num_col_OR
                                                        % Best SNP per gene Odds ratio, lower 95% CI, upper 95% CI
                                                        fprintf(FID_geneset, '%5.4f\t', GWAS_SNP_OR_L95_U95(find_best_SNP_rs,or) );
                                                end
                                        else
                                                for or=1:num_col_OR
                                                        % Best SNP per gene Odds ratio, lower 95% CI, upper 95% CI
                                                        fprintf(FID_geneset, 'NaN\t');
                                                end
                                        end
                                end
        

				if (Flag_gene_set_file_name)
	                         if (Flag_gs == 1)
	 			
                                  if (intersect(GeneSetGeneNames(i),Flag_gs_GeneNames))
                                        fprintf(FID_geneset, '%1.0f\n', 1);    % 1=presence of at least one gene in input gene subset in given gene set
                                  else
                                        fprintf(FID_geneset, '%1.0f\n', 0);    % 0=not present in input gene subset
                                  end
				 end
				else
					fprintf(FID_geneset, '\n');		
				end
			end	

		end

	        fclose(FID_geneset);

	end % print gene scores for input gene sets

    % (9.7) Apply GSEA-GWAS to given gene set

    % GSEA-cutoff
    EnrichPval_multipleCutoffs=[];
    
if (run_GSEA==1)
     
    if (GSEA_method==1)

	clear record_rand_gs_fract obs_fract_above_cutoff record_rand_gs_num obs_num_gs_above_cutoff
	
[EnrichPval_multipleCutoffs,obs_num_gs_above_cutoff,obs_fract_above_cutoff,record_rand_gs_num,record_rand_gs_fract]=GSEA_GWAS_080810(Uncorr_score_isan,Corr_score_isan,top_percen_cutoffs,num_gs_simul,find_gene_set_ind,score_signif_direct,choose_unique_genes,num_perm_limit);
     
        [final_num_simul,cut]=size(record_rand_gs_fract);
        Rand_gs_fract_p95(1:final_num_simul,gs_ind) = record_rand_gs_fract(:,1);
        Rand_gs_fract_p75(1:final_num_simul,gs_ind) = record_rand_gs_fract(:,2);
        obs_fract_above_cutoff_p95(gs_ind) = obs_fract_above_cutoff(1);
        obs_fract_above_cutoff_p75(gs_ind) = obs_fract_above_cutoff(2);
	obs_num_gs_above_cutoff_p95(gs_ind) = obs_num_gs_above_cutoff(1);
	obs_num_gs_above_cutoff_p75(gs_ind) = obs_num_gs_above_cutoff(2);

     elseif (GSEA_method==4)	% GSEA-wilcoxon rank test

		[EnrichPval_multipleCutoffs]=GSEA_GWAS_RankSum_092409(Uncorr_score_isan,Corr_score_isan,top_percen_cutoffs,num_gs_simul,find_gene_set_ind,score_signif_direct,choose_unique_genes);
		top_percen_cutoffs=EnrichPval_multipleCutoffs;

     end

	% concatenate the GSEA-GWAS p-values for all gene sets into one matrix
	% rows: gene sets, columns: GSEA p-values for different gene score cutoffs
	EnrichPval_multipleCutoffs_all_genesets = [EnrichPval_multipleCutoffs_all_genesets; EnrichPval_multipleCutoffs];

	Record_results_all_genesets(gs,1) = Original_gs_size(gs_ind);
    	Record_results_all_genesets(gs,2) = Eff_num_genes_gs(gs_ind);
    	Record_results_all_genesets(gs,3) = Num_genes_no_match_list(gs_ind); 
	Record_results_all_genesets(gs,4) = num_genes_gs_wo_score(gs_ind);
	Record_results_all_genesets(gs,5) = num_genes_removed_overlap(gs_ind);
	Record_results_all_genesets(gs,6) = Presence_flag_gs_genes(gs_ind);
    	Record_results_all_genesets(gs,7) = Med_genesize_gs(gs_ind);
   	Record_results_all_genesets(gs,8) = Mean_genesize_gs(gs_ind);
    	Record_results_all_genesets(gs,9) = Num_LargeGenes_gs(gs_ind);
    	Record_results_all_genesets(gs,10) = Num_SmallGenes_gs(gs_ind);

        if (GSEA_method==1)  % cutoff-based GSEA
	    	if (EnrichPval_multipleCutoffs)
    			Record_results_all_genesets(gs,11:(length(top_percen_cutoffs)+10)) = EnrichPval_multipleCutoffs;
    		end
        elseif (GSEA_method==4) % rank-sum-based GSEA
                if (EnrichPval_multipleCutoffs)
                        Record_results_all_genesets(gs,11) = EnrichPval_multipleCutoffs(1);
                end 
	end
                
	% Calculate FDR when finish iterating over all gene-sets of a given database; a minimum of 10 gene sets required
                        
	if (calculate_FDR==1)
           
	        % if database changes between this gene set and the following gene
        	% set or if this is the last gene set in the input list calculate
	        % FDR for all gene sets analyzed until now.

        	if ( (strcmp(gs_GeneIDs_cellarray_min_gene_num{gs}{1},gs_GeneIDs_cellarray_min_gene_num{gs+1}{1})==0) | gs==num_gene_sets_min_num_genes) 

			% Calculate FDR q-value for each gene set, a p-value that corrects for multiple hypothesis testing
			score_cutoffs=prctile(Corr_score_isan,top_percen_cutoffs);

			clear Norm_rand_gs_fract_p95 Norm_rand_gs_fract_p75 Norm_obs_fract_above_cutoff_p95 Norm_obs_fract_above_cutoff_p75             
			Norm_rand_gs_fract_p95=zeros(num_gs_simul,gs_ind);
			Norm_rand_gs_fract_p75=zeros(num_gs_simul,gs_ind);
			Norm_obs_fract_above_cutoff_p95 = zeros(1,gs_ind);
			Norm_obs_fract_above_cutoff_p75 = zeros(1,gs_ind);
         
			clear FDR_gs_p95 FDR_gs_p75 FWER_p95 FWER_p75 
			FDR_gs_p95 = zeros(1,gs_ind); % false discovery rate
			FDR_gs_p75 =zeros(1,gs_ind); % family-wise error rate
%			FWER_p95 =zeros(1,gs_ind);
%			FWER_p75 =zeros(1,gs_ind);

			for each_gs=1:gs_ind

			        Norm_rand_gs_fract_p95(:,each_gs) = (Rand_gs_fract_p95(:,each_gs) - mean(Rand_gs_fract_p95(:,each_gs))) / (std(Rand_gs_fract_p95(:,each_gs)));
			        Norm_obs_fract_above_cutoff_p95(each_gs) = (obs_fract_above_cutoff_p95(each_gs) - mean(Rand_gs_fract_p95(:,each_gs))) / (std(Rand_gs_fract_p95(:,each_gs)));
        
			        Norm_rand_gs_fract_p75(:,each_gs) = (Rand_gs_fract_p75(:,each_gs) - mean(Rand_gs_fract_p75(:,each_gs))) / (std(Rand_gs_fract_p75(:,each_gs)));
			        Norm_obs_fract_above_cutoff_p75(each_gs) = (obs_fract_above_cutoff_p75(each_gs) - mean(Rand_gs_fract_p75(:,each_gs))) / (std(Rand_gs_fract_p75(:,each_gs)));
                                                                                        
			end

			
			% Print FDR q-values into output file:                       
			[n_simul,n_gs]=size(Norm_rand_gs_fract_p95);
			num_rand_gs = n_simul*n_gs;
        
			for each_gs=1:gs_ind

				% Calculate FDR per gene set        
        			FDR_gs_p95(each_gs) = (length(find(Norm_rand_gs_fract_p95>=Norm_obs_fract_above_cutoff_p95(each_gs)))/num_rand_gs) / (length(find(Norm_obs_fract_above_cutoff_p95>=Norm_obs_fract_above_cutoff_p95(each_gs))) / length(Norm_obs_fract_above_cutoff_p95) );
			        FDR_gs_p75(each_gs) = (length(find(Norm_rand_gs_fract_p75>=Norm_obs_fract_above_cutoff_p75(each_gs)))/num_rand_gs) / (length(find(Norm_obs_fract_above_cutoff_p75>=Norm_obs_fract_above_cutoff_p75(each_gs))) / length(Norm_obs_fract_above_cutoff_p75) );
  
				% Calculated a family-wise error rate (FWER) for each gene-set
%       			FWER_p95(each_gs) = length(find(Norm_rand_gs_fract_p95>=Norm_obs_fract_above_cutoff_p95(each_gs)))/num_rand_gs;
%			        FWER_p75(each_gs) = length(find(Norm_rand_gs_fract_p75>=Norm_obs_fract_above_cutoff_p75(each_gs)))/num_rand_gs;
        
%        			if (FWER_p95(each_gs)==0)
%			                FWER_p95(each_gs) =0.99/num_rand_gs; % lower bound
%			        end
%			        if (FWER_p75(each_gs)==0)
%			                FWER_p75(each_gs) =0.99/num_rand_gs; % lower bound
%			        end


				if (FDR_gs_p95(each_gs) > 1) 
					FDR_gs_p95(each_gs) =1; 
				end

                                if (FDR_gs_p75(each_gs)	> 1) 
                                        FDR_gs_p75(each_gs) =1;    
                                end

			        % Print Gene Set Enrichment results for each gene set into output file
                                
        			fprintf(FID_results, '%s\t', [num2str(db_name{each_gs})]);               % Resource name
			        fprintf(FID_results, '%s\t', [num2str(GeneSetLabel_array{each_gs})]);          % Gene set name
			        fprintf(FID_results, '%1.0f\t', Original_gs_size(each_gs));          % Original gene set size
			        fprintf(FID_results, '%1.0f\t', Eff_num_genes_gs(each_gs));              % Effective gene set size

			        fprintf(FID_results, '%1.0f\t', Num_genes_no_match_list(each_gs));       % # genes not found on our all human gene list
			        fprintf(FID_results, '%1.0f\t', num_genes_gs_wo_score(each_gs));         % # of genes that matched our gene list that did not get a corrected score
			        fprintf(FID_results, '%1.0f\t', num_genes_removed_overlap(each_gs));     % # genes removed due to best SNP per gene overlap
        
			        fprintf(FID_results, '%1.0f\t', Med_genesize_gs(each_gs));       % median gene size of all genes in gene set in kb units
			        fprintf(FID_results, '%1.0f', Mean_genesize_gs(each_gs));      % mean gene size of all genes in gene set in kb units  


			    if (find(EnrichPval_multipleCutoffs))
    
   			      if (GSEA_method==1)
        
			        fprintf(FID_results,'\t%E', EnrichPval_multipleCutoffs_all_genesets(each_gs,1));   % GSEA-GWAS p-value, 95 percentile cutoff
                                fprintf(FID_results,'\t%E', FDR_gs_p95(each_gs)); % FDR q-value, p-value correcting for multiple hypothesis testing
        	                fprintf(FID_results,'\t%1.0f', round((top_percen_cutoffs(1)/100)*Eff_num_genes_gs(each_gs)));   % expected number of genes with score above percentile cutoffs
 	        	        fprintf(FID_results,'\t%1.0f', obs_num_gs_above_cutoff_p95(each_gs));   % observed number of genes with score above percentile cutoffs
%			        fprintf(FID_results,'\t%E', FWER_p95(each_gs)); % FDR q-value, p-value correcting for multiple hypothesis testing
                
			        fprintf(FID_results,'\t%E', EnrichPval_multipleCutoffs_all_genesets(each_gs,2));   % GSEA-GWAS p-value, 75 percentile cutoff
                                fprintf(FID_results,'\t%E', FDR_gs_p75(each_gs)); % FDR q-value, p-value correcting for multiple hypothesis testing
				fprintf(FID_results,'\t%1.0f', round((top_percen_cutoffs(2)/100)*Eff_num_genes_gs(each_gs)));   % expected number of genes with score above percentile cutoffs
	                        fprintf(FID_results,'\t%1.0f', obs_num_gs_above_cutoff_p75(each_gs));   % observed number of genes with score above percentile cutoffs
%				fprintf(FID_results,'\t%E', FWER_p75(each_gs)); % FDR q-value, p-value correcting for multiple hypothesis testing
   
		               elseif (GSEA_method==4)
                		fprintf(FID_results,'\t%E', EnrichPval_multipleCutoffs_all_genesets(each_gs,1) );   % GSEA rank sum p-value
        		       end

     			     end

                            fprintf(FID_results, '\t%1.0f', Presence_flag_gs_genes(each_gs));        % number of flagged genes in given gene set
	
			% PRINT GENE NAMES OF ALL DISEASE GENES (FLAGGED GENES) IN A GIVEN GENE SET
			if (Flag_gene_set_file_name)

			if (Flag_gs == 1)
			  if (print_trait_gene_names_in_gs==1)
			     if (intersect(GeneSetGeneNames_cell_of_arrays{each_gs},Flag_gs_GeneNames))
        
			                clear Flagged_GeneNames_in_GeneSet
			                [Flagged_GeneNames_in_GeneSet,n1,n2]=intersect(GeneSetGeneNames_cell_of_arrays{each_gs},Flag_gs_GeneNames);
     
					fprintf(FID_results, '\t');
			                for u=1:(length(Flagged_GeneNames_in_GeneSet)-1)
			                        clear find_flag_gene_name Flag_gene_name
			                        [m1,find_flag_gene_name,m2]=intersect(AllRefFlatGeneID, Flagged_GeneNames_in_GeneSet(u));

 		                         if (find_flag_gene_name) % if gene ID appears on human gene list
			                        Flag_gene_name = AllRefFlatGeneNames{find_flag_gene_name(1)};
			                        fprintf(FID_results, '%s', num2str(Flag_gene_name));  % gene name
					        fprintf(FID_results, ', ');
			      %                  fprintf(FID_results, '|');
			      %                  fprintf(FID_results, '%1.0f', Flagged_GeneNames_in_GeneSet(u));  % gene ID
					 end	
			                end
					% for last gene
                                        clear find_flag_gene_name Flag_gene_name
                                        [m1,find_flag_gene_name,m2]=intersect(AllRefFlatGeneID, Flagged_GeneNames_in_GeneSet(length(Flagged_GeneNames_in_GeneSet)));
                                        if (find_flag_gene_name) % if gene ID appears on human gene list
                                                Flag_gene_name = AllRefFlatGeneNames{find_flag_gene_name(1)};
                                                fprintf(FID_results, '%s', num2str(Flag_gene_name));  % gene name
                              %                  fprintf(FID_results, '|');
                              %                  fprintf(FID_results, '%1.0f', Flagged_GeneNames_in_GeneSet(length(Flagged_GeneNames_in_GeneSet) ));  % gene ID
                                         end


			      end
			  end   
			end
			end
			fprintf(FID_results,'\n');

		   end % for each gene set within a given database
			                                                                                        
       
                        gs_ind = 0; % start counting gene sets from next database starting with one

			clear db_name Original_gs_size Eff_num_genes_gs Med_genesize_gs Mean_genesize_gs Num_LargeGenes_gs Num_SmallGenes_gs
			clear Num_genes_no_match_list num_genes_gs_wo_score Presence_flag_gs_genes num_genes_removed_overlap GeneSetLabel_array GeneSetGeneNames_cell_of_arrays
			clear Rand_gs_fract_p95 Rand_gs_fract_p75 obs_fract_above_cutoff_p95 obs_fract_above_cutoff_p75 obs_num_gs_above_cutoff_p95 obs_num_gs_above_cutoff_p75
			clear GeneSetGeneNames_cell_of_arrays

                        clear Norm_rand_gs_fract_p95 Norm_rand_gs_fract_p75 Norm_obs_fract_above_cutoff_p95 Norm_obs_fract_above_cutoff_p75

			EnrichPval_multipleCutoffs_all_genesets=[]; % matrix that stores GSEA-GWAS p-values for all gene sets


        end % calculate FDR for each database separately

    else % do not calculate FDR, print just nominal GSEA p-values
     
	% (9.8) Print Gene Set Enrichment results for each gene set into output file

    	fprintf(FID_results, '%s\t', [num2str(db_name{gs_ind})]); 		% Resource name
	fprintf(FID_results, '%s\t', [num2str(GeneSetLabel_array{gs_ind})]);		% Gene set name
	fprintf(FID_results, '%1.0f\t', Original_gs_size(gs_ind));		% Original gene set size
    	fprintf(FID_results, '%1.0f\t', Eff_num_genes_gs(gs_ind));   		% Effective gene set size

    	fprintf(FID_results, '%1.0f\t', Num_genes_no_match_list(gs_ind)); 	% # genes not found on our all human gene list
    	fprintf(FID_results, '%1.0f\t', num_genes_gs_wo_score(gs_ind)); 	% # of genes that matched our gene list that did not get a corrected score
    	fprintf(FID_results, '%1.0f\t', num_genes_removed_overlap(gs_ind));	% # genes removed due to best SNP per gene overlap        

    	fprintf(FID_results, '%1.0f\t', Med_genesize_gs(gs_ind));	% median gene size of all genes in gene set in kb units  
    	fprintf(FID_results, '%1.0f', Mean_genesize_gs(gs_ind));	% mean gene size of all genes in gene set in kb units
%    	fprintf(FID_results, '%1.0f\t', Num_LargeGenes_gs(gs_ind));	% Number of genes in gene set with size >=100 kb
%    	fprintf(FID_results, '%1.0f', Num_SmallGenes_gs(gs_ind));	% Number of genes in gene set with size <=10 kb
                
     	if (find(EnrichPval_multipleCutoffs))

         if (GSEA_method==1)

                for y=1:length(EnrichPval_multipleCutoffs)
                        fprintf(FID_results,'\t%E', EnrichPval_multipleCutoffs(y));   % GSEA-GWAS p-value for different percentile cutoffs
			fprintf(FID_results,'\t%1.0f', round((top_percen_cutoffs(y)/100)*Eff_num_genes_gs(gs_ind)));   % expected number of genes with score above percentile cutoffs
			fprintf(FID_results,'\t%1.0f', obs_num_gs_above_cutoff(y));   % observed number of genes with score above percentile cutoffs
                end
       
         elseif (GSEA_method==4)
                fprintf(FID_results,'\t%E', EnrichPval_multipleCutoffs(1,1));   % GSEA rank sum p-value
         end
     	end

        fprintf(FID_results, '\t%1.0f', Presence_flag_gs_genes(gs_ind));        % number of flagged genes in given gene set

                   % PRINT GENE NAMES OF ALL DISEASE GENES (FLAGGED GENES) IN A GIVEN GENE SET
                        if (Flag_gene_set_file_name)

                        if (Flag_gs == 1)
			  if (print_trait_gene_names_in_gs==1)
                             if (intersect(GeneSetGeneNames_cell_of_arrays{gs_ind},Flag_gs_GeneNames))
                        
                                        clear Flagged_GeneNames_in_GeneSet
                                        [Flagged_GeneNames_in_GeneSet,n1,n2]=intersect(GeneSetGeneNames_cell_of_arrays{gs_ind},Flag_gs_GeneNames);
                                    	fprintf(FID_results, '\t');
                                        for u=1:(length(Flagged_GeneNames_in_GeneSet)-1)
                                                clear find_flag_gene_name Flag_gene_name
                                                [m1,find_flag_gene_name,m2]=intersect(AllRefFlatGeneID, Flagged_GeneNames_in_GeneSet(u));

                                         if (find_flag_gene_name) % if gene ID appears on human gene list
                                                Flag_gene_name = AllRefFlatGeneNames{find_flag_gene_name(1)};
                                                fprintf(FID_results, '%s', num2str(Flag_gene_name));  % gene name
                                                fprintf(FID_results, ', ');
                              %                  fprintf(FID_results, '|');
                              %                  fprintf(FID_results, '%1.0f', Flagged_GeneNames_in_GeneSet(u));  % gene ID
                                         end
                                        end
                                        % for last gene
                                        clear find_flag_gene_name Flag_gene_name
					[m1,find_flag_gene_name,m2]=intersect(AllRefFlatGeneID, Flagged_GeneNames_in_GeneSet(length(Flagged_GeneNames_in_GeneSet)));

                                        if (find_flag_gene_name) % if gene ID appears on human gene list
                                                Flag_gene_name = AllRefFlatGeneNames{find_flag_gene_name(1)};
                                                fprintf(FID_results, '%s', num2str(Flag_gene_name));  % gene name
                              %                  fprintf(FID_results, '|');
                              %                  fprintf(FID_results, '%1.0f', Flagged_GeneNames_in_GeneSet(length(Flagged_GeneNames_in_GeneSet) ));  % gene ID
                                         end


                              end
                          end
			end
			end
                        fprintf(FID_results,'\n');

%      end % for each selected gene set

    end % if calculate FDR

end % if calculate gene set enrichment

end  % for each gene set

if (run_GSEA==1)

 clear save_com
 save_com = ['save ', num2str(Output_dir) ,'/GeneSetEnrichPval_AllGeneSets_', num2str(exp_label) , '_',  num2str(Gene_boundr_upstr/1000), 'kb_upstr_', num2str(Gene_boundr_downstr/1000), 'kb_downstr_', num2str(num_gs_simul), 'perm_', num2str(todays_date) , ' EnrichPval_multipleCutoffs_all_genesets'];
 eval(save_com);

% clear save_com
%  save_com = ['save ', num2str(Output_dir) ,'/Record_results_all_genesets_', num2str(exp_label) , '_',  num2str(Gene_boundr_upstr/1000), 'kb_upstr_', num2str(Gene_boundr_downstr/1000), 'kb_downstr_', num2str(num_gs_simul), 'perm_', num2str(todays_date) , ' Record_results_all_genesets -ascii -tabs'];
% eval(save_com);

 fclose(FID_results);

end % if calculate gene set enrichment

% (10) Print into LOG FILE time it took to run

Run_time=toc

fprintf(FID_log, ['Time it took the program to run: ', num2str(Run_time/60), ' minutes = ', num2str(Run_time/3600), ' hours = ', num2str(Run_time/(3600*24)), ' days.\n']);

fclose(FID_log);

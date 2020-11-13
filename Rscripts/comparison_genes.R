### Load packages ###
library(cmapR)


#### Define variables #### 
place = "home" #home or work
type_genes = 'curated' # curated or all

if (place=="work"){
  main_directory = "/home/quim/PHD/Projects/camda"
  references_directory = "/home/quim/PHD/References"
  bigmem_directory = "/sbi/users/interchange/emre/quim/camda"
} else {
  main_directory = "/Users/quim/Dropbox/UPF/PHD/Projects/camda"
  references_directory = "/Users/quim/Dropbox/UPF/PHD/References"
  bigmem_directory = "/Users/quim/Documents/DATA/camda"
}


#### Define directories and files #### 
# Data files
drugs_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-pert_iname.rda", sep="/")
drug_info_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-info.rda", sep="/")
dilirank_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-dilirank.v2.rda", sep="/")
expression_file <- paste(bigmem_directory, "CAMDA_l1000_1314compounds-GSE92742_Level5_gct.rda", sep="/")
gene_info_file <- paste(main_directory, "additional_data/GSE92742_Broad_LINCS_gene_info.txt", sep="/")
# Specific data files
wilcox_file <- paste(main_directory, "results/reverse_engineering/reverse_signature_phh_notcorrected.txt", sep="/")
phenotype2gene_file <- paste(main_directory, "guildify_data/phenotype2gene.tsv", sep="/")
disease2gene_file <- paste(main_directory, "guildify_data/disease2gene_disgenet_guildify.tsv", sep="/")
redundantphenotypes_file <- paste(main_directory, "guildify_data/redundant_phenotypes.tsv", sep="/")
guildify_dir <- paste(main_directory, "guildify_data/guildify_results/", sep="/")
peng_file <- paste(references_directory, "Peng_ToxLetters19_supplementary/TableS5.txt", sep="/")
# Output files
output_peng_number_shared_genes_table <- paste(main_directory, "outputs/files/Peng_number_of_shared_genes.tsv", sep="/")
output_peng_list_shared_genes_table <- paste(main_directory, "outputs/files/Peng_list_of_shared_genes.tsv", sep="/")
output_phenotype_to_genes_table <- paste(main_directory, "outputs/files/curated_phenotypes_to_genes.tsv", sep="/")


### Load files ###
#source(functions_file)
#load(drugs_file)
#load(drug_info_file)
#load(dilirank_file)
#load(expression_file) # Requires cmapR


#### Subset drugs ####
drug.dataset <- subset.drug.dataset(drank.sel, outliers=outliers, remove.outliers=remove.outliers)


#### Read Peng et al. table#### 
peng_df <- read.csv(peng_file, header=TRUE, sep="\t", stringsAsFactors = F)


#### Prepare output table ####
cols <- c("Dataset", "Num. genes", "Num. shared genes with Peng dataset")
peng_shared_df <- data.frame(matrix(ncol = length(cols), nrow=0))
colnames(peng_shared_df) <- cols


#### Get landmark genes ####
gene_info_df <- read.csv(gene_info_file, header=TRUE, sep="\t", stringsAsFactors = F)
landmark_genes <- gene_info_df$pr_gene_id[gene_info_df$pr_is_lm==1]
landmark_genes_info <- gene_info_df[gene_info_df$pr_is_lm==1,]
shared_genes_landmark_peng <- landmark_genes_info[landmark_genes_info$pr_gene_symbol %in% peng_df$Hepatotoxicity,]
print(nrow(shared_genes_landmark_peng))
print(shared_genes_landmark_peng$pr_gene_symbol[order(shared_genes_landmark_peng$pr_gene_symbol)])
peng_shared_df[nrow(peng_shared_df)+1,] <- c( "Landmark", length(landmark_genes), length(shared_genes_landmark_peng$pr_gene_symbol) )


#### Read wilcoxon signature table and get overlap with Peng genes#### 
wilcox_df <- read.csv(wilcox_file, header=FALSE, sep="\t", stringsAsFactors = F)
wilcox_genes <- unique(wilcox_df[,1])
wilcox_genes_info <- gene_info_df[gene_info_df$pr_gene_id %in% wilcox_genes,]
shared_genes_wilcox_peng <- wilcox_genes_info[wilcox_genes_info$pr_gene_symbol %in% peng_df$Hepatotoxicity,]
print(nrow(shared_genes_wilcox_peng))
print(shared_genes_wilcox_peng$pr_gene_symbol[order(shared_genes_wilcox_peng$pr_gene_symbol)])
peng_shared_df[nrow(peng_shared_df)+1,] <- c( "DILI Landmark", length(wilcox_genes_info$pr_gene_symbol), length(shared_genes_wilcox_peng$pr_gene_symbol) )


#### Read genes associated to DisGeNET/GUILDify phenotypes and get overlap with Peng genes ####

# Read phenotypes to gene data
disgenet2gene <- read.csv(disease2gene_file, header=TRUE, sep="\t", stringsAsFactors = F)
redundantphenotypes_df <- read.csv(redundantphenotypes_file, header=TRUE, sep="\t", stringsAsFactors = F)
categories <- c(paste("disgenet.",type_genes, sep = ""), paste("guildify.",type_genes, sep = ""))
redundantphenotypes_df <- redundantphenotypes_df[redundantphenotypes_df$category %in% categories,]
disgenet2gene <- disgenet2gene[!(disgenet2gene$diseaseid %in% redundantphenotypes_df$redundant.phenotype),]
phenotypes <- unique(disgenet2gene[,c('diseaseid', 'diseaseterm')])
phenotypes <- phenotypes[order(phenotypes$diseaseterm),]

# Get genes from disgenet and guildify and calculate the overlap for each phenotype
cols <- c("DILI Phenotypes", "UMLS", "Number of genes associated in DisGeNET", "Number of genes associated in GUILDify")
phenotypes_to_genes_df <- data.frame(matrix(ncol = length(cols), nrow=0))
colnames(phenotypes_to_genes_df) <- cols
all_disgenet_genes <- c()
all_guildify_genes <- c()
disgenet_phenotype_to_genes <- list()
guildify_phenotype_to_genes <- list()
for (i in 1:nrow(phenotypes)){
  diseaseid <- phenotypes[i,c("diseaseid")]
  diseaseterm <- phenotypes[i,c("diseaseterm")]
  disgenet2gene_df <- disgenet2gene[c("geneid", "diseaseid", "diseaseterm", paste("disgenet", type_genes, sep='.'))]
  colnames(disgenet2gene_df) <- c("geneid", "diseaseid", "diseaseterm", "type_analysis")
  disgenet_genes <- unique(disgenet2gene_df$geneid[disgenet2gene_df$type_analysis==1 & disgenet2gene_df$diseaseid==diseaseid])
  guildify2gene_df <- disgenet2gene[c("geneid", "diseaseid", "diseaseterm", paste("guildify", type_genes, sep='.'))]
  colnames(guildify2gene_df) <- c("geneid", "diseaseid", "diseaseterm", "type_analysis")
  guildify_genes <- unique(guildify2gene_df$geneid[guildify2gene_df$type_analysis==1 & guildify2gene_df$diseaseid==diseaseid])
  if (length(disgenet_genes)>=10) {
    # Annotate genes of DisGeNET and GUILDify
    phenotypes_to_genes_df[nrow(phenotypes_to_genes_df)+1,] <- c( diseaseterm, diseaseid, length(disgenet_genes), length(guildify_genes) )
    all_disgenet_genes <- unique(c(all_disgenet_genes, disgenet_genes))
    all_guildify_genes <- unique(c(all_guildify_genes, guildify_genes))
    #### Get genes from DisGeNET shared with Peng ####
    disgenet_genes_info_df <- gene_info_df[gene_info_df$pr_gene_id %in% disgenet_genes,]
    shared_genes_disgenet_peng_df <- disgenet_genes_info_df[disgenet_genes_info_df$pr_gene_symbol %in% peng_df$Hepatotoxicity,]
    peng_shared_df[nrow(peng_shared_df)+1,] <- c( paste(diseaseterm, '(DisGeNET)', sep=" "), length(disgenet_genes), length(shared_genes_disgenet_peng_df$pr_gene_symbol) )
    disgenet_phenotype_to_genes[[diseaseterm]] <- list(disgenet_genes)
    #### Get genes from GUILDify shared with Peng ####
    if (length(guildify_genes) < length(disgenet_genes)){
      next
    }
    guildify_genes_info_df <- gene_info_df[gene_info_df$pr_gene_id %in% guildify_genes,]
    shared_genes_guildify_peng_df <- guildify_genes_info_df[guildify_genes_info_df$pr_gene_symbol %in% peng_df$Hepatotoxicity,]
    peng_shared_df[nrow(peng_shared_df)+1,] <- c( paste(diseaseterm, '(GUILDify)', sep=" "), length(guildify_genes), length(shared_genes_guildify_peng_df$pr_gene_symbol) )
    guildify_phenotype_to_genes[[diseaseterm]] <- list(guildify_genes)
  }
}
phenotypes_to_genes_df[nrow(phenotypes_to_genes_df)+1,] <- c("Number of different genes associated to DILI phenotypes", "-", length(all_disgenet_genes), length(all_guildify_genes))


#### Make Peng's list of genes ####
cols <- c("Gene Symbol", "Gene ID", "Included in the following datasets")
peng_shared_list_df <- data.frame(matrix(ncol = length(cols), nrow=0))
colnames(peng_shared_list_df) <- cols
peng_genes <- unique(peng_df$Hepatotoxicity[order(peng_df$Hepatotoxicity)])
for (peng_gene in peng_genes){
  gene_id <- "-"
  datasets <- c()
  if(peng_gene %in% gene_info_df$pr_gene_symbol){
    gene_id <- gene_info_df$pr_gene_id[gene_info_df$pr_gene_symbol == peng_gene]
  }
  if(peng_gene %in% shared_genes_landmark_peng$pr_gene_symbol){
    datasets <- c(datasets, "Landmark")
  }
  if(peng_gene %in% shared_genes_wilcox_peng$pr_gene_symbol){
    datasets <- c(datasets, "DILI Landmark")
  }
  for(phenotype in names(disgenet_phenotype_to_genes)){
    if(gene_id %in% disgenet_phenotype_to_genes[[phenotype]][[1]]){
      datasets <- c(datasets, paste(phenotype, " (DisGeNET)"))
    }
    if(phenotype %in% names(guildify_phenotype_to_genes)){
      if(gene_id %in% guildify_phenotype_to_genes[[phenotype]][[1]]){
        datasets <- c(datasets, paste(phenotype, " (GUILDify)"))
      }
    }
  }
  if(length(datasets) == 0){
    datasets <- "-"
  }
  peng_shared_list_df[nrow(peng_shared_list_df)+1,] <- c(peng_gene, gene_id, paste(datasets, collapse = "; "))
}


#### Write output tables ####
write.table(peng_shared_df, file = output_peng_number_shared_genes_table, row.names=FALSE, na="-", col.names=TRUE, sep="\t", quote=FALSE)
write.table(peng_shared_list_df, file = output_peng_list_shared_genes_table, row.names=FALSE, na="-", col.names=TRUE, sep="\t", quote=FALSE)
write.table(phenotypes_to_genes_df, file = output_phenotype_to_genes_table, row.names=FALSE, na="-", col.names=TRUE, sep="\t", quote=FALSE)


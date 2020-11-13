### Load packages ###
library(cmapR)
library(ggplot2)
library(caret)


### Define variables ###
remove.outliers = FALSE
outliers = c('daunorubicin', 'vorinostat')
number.cv = 10
number.repetitions = 10
fraction_train = 0.7
place = "home" #home or work
if (place=="work"){
  main_directory = "/home/quim/PHD/Projects/camda"
  bigmem_directory = "/sbi/users/interchange/emre/quim/camda"
} else {
  main_directory = "/Users/quim/Dropbox/UPF/PHD/Projects/camda"
  bigmem_directory = "/Users/quim/Documents/DATA/camda"
}


### Define files ###
# Data files
drugs_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-pert_iname.rda", sep="/")
drug_info_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-info.rda", sep="/")
dilirank_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-dilirank.v2.rda", sep="/")
expression_file <- paste(bigmem_directory, "CAMDA_l1000_1314compounds-GSE92742_Level5_gct.rda", sep="/")
landmark_genes_file <- paste(main_directory, "additional_data/GSE92742_Broad_LINCS_gene_info_delta_landmark.txt", sep="/")
gene_info_file <- paste(main_directory, "additional_data/GSE92742_Broad_LINCS_gene_info.txt", sep="/")
cell_info_file <- paste(main_directory, "additional_data/GSE92742_Broad_LINCS_cell_info.txt", sep="/")
functions_file <- paste(main_directory, "Rscripts/camda_functions.R", sep="/")
# Specific data files
phenotype2gene_file <- paste(main_directory, "guildify_data/phenotype2gene.tsv", sep="/")
guildify_dir <- paste(main_directory, "guildify_data/guildify_results/", sep="/")
# Output files
disease2gene_file <- paste(main_directory, "guildify_data/disease2gene_disgenet_guildify.tsv", sep="/")
redundantphenotypes_file <- paste(main_directory, "guildify_data/redundant_phenotypes.tsv", sep="/")


### Get phenotypes from DisGeNET ###
phenotype2gene <- read.csv(phenotype2gene_file, header=TRUE, sep="\t", stringsAsFactors = F)
phenotypes <- unique(phenotype2gene$diseaseid)
all_phenotypes <- unique(phenotype2gene$diseaseid[phenotype2gene$source == "ALL"])
curated_phenotypes <- unique(phenotype2gene$diseaseid[phenotype2gene$source == "CURATED"])
#curated_phenotypes <- c("C0023890", "C0239946")


### Get genes from DisGeNET and GUILDify for each phenotype ###
cols <- c("geneid", "diseaseid", "diseaseterm", "disgenet.curated", "disgenet.all", "guildify.curated", "guildify.all")
diseasegene_df <- data.frame(matrix(ncol = length(cols), nrow=0))
colnames(diseasegene_df) <- cols
for (phenotype in all_phenotypes){
  # Get genes from DisGeNET (ALL)
  diseaseterm <- unique(phenotype2gene$name[phenotype2gene$diseaseid==phenotype])
  disgenet_info <- phenotype2gene[c("geneid")][phenotype2gene$diseaseid==phenotype & phenotype2gene$source=="ALL",]

  # Get genes from GUILDify (ALL)
  guildify_filename = paste(phenotype, "_ALL.tsv", sep="")
  guildify_file <- paste(guildify_dir, guildify_filename, sep="")
  guildify_result <- read.csv(guildify_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  guildify_info_raw <- as.character(guildify_result$gene.id)
  guildify_info <- as.numeric(unlist(strsplit(guildify_info_raw, split = ","))) # Split multiple entries in "gene.id" column
  # Check if the expansion is correct
  if (length(guildify_info) <= 1){
    guildify_info <- c()
  }
  
  # Get genes from DisGeNET (CURATED)
  disgenet_curated <- phenotype2gene[c("geneid")][phenotype2gene$diseaseid==phenotype & phenotype2gene$source=="CURATED",]
  
  # Get genes from GUILDify (CURATED)
  guildify_curated <- c()
  guildify_cur_filename = paste(phenotype, "_CURATED.tsv", sep="")
  guildify_cur_file <- paste(guildify_dir, guildify_cur_filename, sep="")
  if (file.exists(guildify_cur_file)){
    guildify_cur_result <- read.csv(guildify_cur_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    guildify_cur_raw <- as.character(guildify_cur_result$gene.id)
    guildify_curated <- as.numeric(unlist(strsplit(guildify_cur_raw, split = ","))) # Split multiple entries in "gene.id" column
    # Check if the expansion is correct
    if (length(guildify_curated) <= 1){
      guildify_curated <- c()
    }
  }

  # Create final table
  all_genes <- union(union(union(disgenet_info, guildify_info), disgenet_curated), guildify_curated)
  disease_df <- data.frame(geneid=all_genes, diseaseid=phenotype, diseaseterm=diseaseterm)
  disease_df$disgenet.curated <- ifelse(disease_df$geneid %in% disgenet_curated, 1, 0)
  disease_df$disgenet.all <- ifelse(disease_df$geneid %in% disgenet_info, 1, 0)
  disease_df$guildify.curated <- ifelse(disease_df$geneid %in% guildify_curated, 1, 0)
  disease_df$guildify.all <- ifelse(disease_df$geneid %in% guildify_info, 1, 0)
  diseasegene_df <- rbind(diseasegene_df, disease_df)
}

# Write output table
write.table(diseasegene_df, file = disease2gene_file,row.names=FALSE, na="-",col.names=TRUE, sep="\t")


### Find redundant phenotypes ###

#### Get redundant phenotypes (with same genes). ####
get.redundant.phenotypes<-function(diseasegene_table, category, phenotypes, output_table) {
  diseasegene_cat <- diseasegene_table[c("geneid", "diseaseid", category)]
  colnames(diseasegene_cat) <- c("geneid", "diseaseid", "category")
  for (i in 1:length(phenotypes)){
    phenotype1 <- phenotypes[i]
    genes1 <- diseasegene_cat$geneid[diseasegene_cat$category==1 & diseasegene_cat$diseaseid==phenotype1]
    for (j in i:length(phenotypes)){
      phenotype2 <- phenotypes[j]
      if (phenotype1 != phenotype2){
        genes2 <- diseasegene_cat$geneid[diseasegene_cat$category==1 & diseasegene_cat$diseaseid==phenotype2]
        if (!phenotype2 %in% output_table$redundant.phenotype){
          if (setequal(genes1, genes2)){
            output_table[nrow(output_table)+1,] <- c(category, phenotype1, phenotype2)
          }
        }
      }
    }
  }
return(output_table);
}

# Define output table of redundant phenotypes
cols <- c("category", "main.phenotype", "redundant.phenotype")
redundantphenotypes_df <- data.frame(matrix(ncol = length(cols), nrow=0))
colnames(redundantphenotypes_df) <- cols
# Get redundant phenotypes for all categories
redundantphenotypes_df <- get.redundant.phenotypes(diseasegene_df, category = "disgenet.all", phenotypes = all_phenotypes, output_table = redundantphenotypes_df)
redundantphenotypes_df <- get.redundant.phenotypes(diseasegene_df, category = "disgenet.curated", phenotypes = curated_phenotypes, output_table = redundantphenotypes_df)
redundantphenotypes_df <- get.redundant.phenotypes(diseasegene_df, category = "guildify.all", phenotypes = all_phenotypes, output_table = redundantphenotypes_df)
redundantphenotypes_df <- get.redundant.phenotypes(diseasegene_df, category = "guildify.curated", phenotypes = curated_phenotypes, output_table = redundantphenotypes_df)
# Write output table
write.table(redundantphenotypes_df, file = redundantphenotypes_file, row.names=FALSE, na="-",col.names=TRUE, sep="\t")



### Load packages ###
library(cmapR)
library(ggplot2)
library(caret)


### Define variables ###
place = "home" #home or work
remove.outliers = FALSE
outliers = c('daunorubicin', 'vorinostat')
number.cv = 10
number.repetitions = 10
fraction_train = 0.7

if (place=="work"){
  main_directory = "/home/quim/PHD/Projects/camda"
  bigmem_directory = "/sbi/users/interchange/emre/quim/camda"
} else {
  main_directory = "/Users/quim/Dropbox/UPF/PhD/Projects/camda"
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
redundantphenotypes_file <- paste(main_directory, "guildify_data/redundant_phenotypes.tsv", sep="/")
wilcox_file <- paste(main_directory, "results/reverse_engineering/reverse_signature_phh_notcorrected.txt", sep="/")
tanimoto_file <- paste(main_directory, "/additional_data/tanimoto_smiles.tsv", sep="/")
targets_file <- paste(main_directory, "additional_data/targets/targets_dgidb_hitpick_sea.tsv", sep="/")
# Output files
output.cv.rf <- paste(main_directory, "results/crossvalidation/cv_landmark_rf.txt", sep="/")
output.cv.gbm <- paste(main_directory, "results/crossvalidation/cv_landmark_gbm.txt", sep="/")
output.ind.rf = paste(main_directory, "camda_data/independent_validation/JAguirre_predictions_landmark_rf.txt", sep="/")
output.ind.gbm = paste(main_directory, "camda_data/independent_validation/JAguirre_predictions_landmark_gbm.txt", sep="/")


### Load files ###
source(functions_file)
load(drugs_file)
load(drug_info_file)
load(dilirank_file)
load(expression_file) # Requires cmapR


#### Subset drugs ####
drug.dataset <- subset.drug.dataset(drank.sel, outliers=outliers, remove.outliers=remove.outliers)

#### Get landmark genes ####
gene_info_df <- read.csv(gene_info_file, header=TRUE, sep="\t")
landmark_genes <- gene_info_df$pr_gene_id[gene_info_df$pr_is_lm==1]

### Get genes associated to phenotypes from DisGeNET ###
phenotype2gene <- read.csv(phenotype2gene_file, header=TRUE, sep="\t")
phenotypes <- unique(phenotype2gene$diseaseid)
curated_phenotypes <- unique(phenotype2gene$diseaseid[phenotype2gene$source == "CURATED"])

### Prepare SMILES data ###
tanimoto_df <- read.csv(tanimoto_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tanimoto_df$DILIConcern[tanimoto_df$DILIConcern == "No-DILI-concern"] <- "No-DILI-Concern" # Correct the 4 entries with lowercase concern
tanimoto_df$pert_iname = rownames(tanimoto_df)
colnames(tanimoto_df)[match("DILIConcern", colnames(tanimoto_df))] <- "dilirank"
tanimoto_df$vDILIConcern <- NULL
tanimoto_ind_df <- tanimoto_df[colnames(tanimoto_df) %in% drug.dataset$independent_drugs, rownames(tanimoto_df) %in% drug.dataset$drugs]
tanimoto_df<-tanimoto_df[colnames(tanimoto_df) %in% drug.dataset$drugs, rownames(tanimoto_df) %in% drug.dataset$drugs]

### Prepare targets info ###
# Read targets and map them to dilirank info (the dilirank info on targets is wrong!)
targets_df <- read.csv(targets_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(targets_df)[match("drug", colnames(targets_df))] <- "pert_iname" # Change name to pert_iname
targets_ind_df <- targets_df[targets_df$pert_iname %in% drug.dataset$independent_drugs,] # Get independent set
targets_df <- merge(x = targets_df, y = drug.dataset$dilirank_df[c("pert_iname", "DILIConcern", "Severity.Class")], by = "pert_iname")
targets_df$DILI <- NULL
targets_df$severity <- NULL
targets_ind_df$DILI <- NULL
targets_ind_df$severity <- NULL
colnames(targets_df)[match("DILIConcern", colnames(targets_df))] <- "dilirank"
colnames(targets_df)[match("Severity.Class", colnames(targets_df))] <- "severity"
targets_df <- targets_df[targets_df$pert_iname %in% drug.dataset$drugs,]

### Load wilcoxon test analysis ###
wilcox_df = read.csv(wilcox_file, header = FALSE, sep = "\t")
#selected_genes <- unique(wilcox_df$gene_id[wilcox_df$p.value<0.05 & wilcox_df$landmark==TRUE])
selected_genes <- unique(wilcox_df[,1])

### Subset GCT object by landmark genes ###
# Subset the GCT object by cell ID PHH, 10 µM, 24 h
expression_df <- subset.expression(gct, landmark_genes, drug.dataset$drugs, drug.dataset$dilirank_df, cell_id="PHH", pert_idose="10 µM", pert_itime="24 h", merge_samples = TRUE)
# Subset gene expression for independent drugs as well
expression_ind_df <- subset.expression(gct, landmark_genes, drug.dataset$independent_drugs, drug.dataset$independent_df, cell_id="PHH", pert_idose="10 µM", pert_itime="24 h", merge_samples = TRUE)

### Subset GCT object by wilcoxon genes ###
# Subset the GCT object by cell ID PHH, 10 µM, 24 h
expression_wilcox_df <- subset.expression(gct, selected_genes, drug.dataset$drugs, drug.dataset$dilirank_df, cell_id="PHH", pert_idose="10 µM", pert_itime="24 h", merge_samples = TRUE)
# Subset gene expression for independent drugs as well
expression_wilcox_ind_df <- subset.expression(gct, selected_genes, drug.dataset$independent_drugs, drug.dataset$independent_df, cell_id="PHH", pert_idose="10 µM", pert_itime="24 h", merge_samples = TRUE)

### Prepare balanced machine learning datasets ###
datasets.list.landmark <- prepare.balanced.datasets(expression_df, number.repetitions, drug.dataset$most_concern_drugs, drug.dataset$less_concern_drugs, drug.dataset$no_concern_drugs, type_analysis = "discrete", fraction_train=fraction_train)
datasets.list.wilcoxon <- prepare.balanced.datasets(expression_wilcox_df, number.repetitions, drug.dataset$most_concern_drugs, drug.dataset$less_concern_drugs, drug.dataset$no_concern_drugs, type_analysis = "discrete", fraction_train=fraction_train)
datasets.list.smiles <- prepare.balanced.datasets(tanimoto_df, number.repetitions, drug.dataset$most_concern_drugs, drug.dataset$less_concern_drugs, drug.dataset$no_concern_drugs, type_analysis = "discrete", fraction_train=fraction_train)
datasets.list.targets <- prepare.balanced.datasets(targets_df, number.repetitions, drug.dataset$most_concern_drugs, drug.dataset$less_concern_drugs, drug.dataset$no_concern_drugs, type_analysis = "discrete", fraction_train=fraction_train)

### Define balanced drugs for machine learning datasets (disgenet/guildify) ###
balanced.drugs.list <- prepare.balanced.drugs(number.repetitions, drug.dataset$most_concern_drugs, drug.dataset$less_concern_drugs, drug.dataset$no_concern_drugs, fraction_train=fraction_train)

### Count the number of different drugs in each dataset ###
# Create new dataframe to store gene expression values for pca
cols <- c("Type of drug", "Num. DILIrank drugs", "Num. Most drugs", "Num. Less drugs", "Num. No drugs", "Num. independent drugs", "Num. train drugs", "Num. Most train drugs", "Num. Less train drugs", "Num. No train drugs", "Num. test drugs", "Num. Most test drugs", "Num. Less test drugs", "Num. No test drugs")
compounds_df <- data.frame(matrix(ncol = length(cols), nrow=0))
colnames(compounds_df) <- cols
compounds_df[nrow(compounds_df)+1,] <- c( "Landmark", length(expression_df$pert_iname), length(expression_df[expression_df$dilirank == "Most-DILI-Concern",]$pert_iname), length(expression_df[expression_df$dilirank == "Less-DILI-Concern",]$pert_iname), length(expression_df[expression_df$dilirank == "No-DILI-Concern",]$pert_iname), length(expression_ind_df$pert_iname),
                                          length(datasets.list.landmark$training_datasets_with_names[[1]]$pert_iname), length(datasets.list.landmark$training_datasets_with_names[[1]][datasets.list.landmark$training_datasets_with_names[[1]]$dilirank == "Most-DILI-Concern",]$dilirank), length(datasets.list.landmark$training_datasets_with_names[[1]][datasets.list.landmark$training_datasets_with_names[[1]]$dilirank == "Less-DILI-Concern",]$dilirank), length(datasets.list.landmark$training_datasets_with_names[[1]][datasets.list.landmark$training_datasets_with_names[[1]]$dilirank == "No-DILI-Concern",]$dilirank), 
                                          length(datasets.list.landmark$testing_with_names$pert_iname), length(datasets.list.landmark$testing_with_names$dilirank[datasets.list.landmark$testing_with_names$dilirank == "Most-DILI-Concern"]), length(datasets.list.landmark$testing_with_names$dilirank[datasets.list.landmark$testing_with_names$dilirank == "Less-DILI-Concern"]), length(datasets.list.landmark$testing_with_names$dilirank[datasets.list.landmark$testing_with_names$dilirank == "No-DILI-Concern"]) )
compounds_df[nrow(compounds_df)+1,] <- c( "DILI Landmark", length(expression_wilcox_df$pert_iname), length(expression_wilcox_df[expression_wilcox_df$dilirank == "Most-DILI-Concern",]$pert_iname), length(expression_wilcox_df[expression_wilcox_df$dilirank == "Less-DILI-Concern",]$pert_iname), length(expression_wilcox_df[expression_wilcox_df$dilirank == "No-DILI-Concern",]$pert_iname), length(expression_wilcox_ind_df$pert_iname),
                                          length(datasets.list.wilcoxon$training_datasets_with_names[[1]]$pert_iname), length(datasets.list.wilcoxon$training_datasets_with_names[[1]][datasets.list.wilcoxon$training_datasets_with_names[[1]]$dilirank == "Most-DILI-Concern",]$dilirank), length(datasets.list.wilcoxon$training_datasets_with_names[[1]][datasets.list.wilcoxon$training_datasets_with_names[[1]]$dilirank == "Less-DILI-Concern",]$dilirank), length(datasets.list.wilcoxon$training_datasets_with_names[[1]][datasets.list.wilcoxon$training_datasets_with_names[[1]]$dilirank == "No-DILI-Concern",]$dilirank), 
                                          length(datasets.list.wilcoxon$testing_with_names$pert_iname), length(datasets.list.wilcoxon$testing_with_names$dilirank[datasets.list.wilcoxon$testing_with_names$dilirank == "Most-DILI-Concern"]), length(datasets.list.wilcoxon$testing_with_names$dilirank[datasets.list.wilcoxon$testing_with_names$dilirank == "Less-DILI-Concern"]), length(datasets.list.wilcoxon$testing_with_names$dilirank[datasets.list.wilcoxon$testing_with_names$dilirank == "No-DILI-Concern"]) )
compounds_df[nrow(compounds_df)+1,] <- c( "DisGeNET", length(drug.dataset$drugs), length(drug.dataset$most_concern_drugs), length(drug.dataset$less_concern_drugs), length(drug.dataset$no_concern_drugs), length(drug.dataset$independent_drugs),
                                          length(balanced.drugs.list$drugs_training_list[[1]]), length(balanced.drugs.list$drugs_training_list[[1]][balanced.drugs.list$drugs_training_list[[1]] %in% drug.dataset$most_concern_drugs]), length(balanced.drugs.list$drugs_training_list[[1]][balanced.drugs.list$drugs_training_list[[1]] %in% drug.dataset$less_concern_drugs]), length(balanced.drugs.list$drugs_training_list[[1]][balanced.drugs.list$drugs_training_list[[1]] %in% drug.dataset$no_concern_drugs]),
                                          length(balanced.drugs.list$drugs_testing), length(balanced.drugs.list$drugs_testing[balanced.drugs.list$drugs_testing %in% drug.dataset$most_concern_drugs]), length(balanced.drugs.list$drugs_testing[balanced.drugs.list$drugs_testing %in% drug.dataset$less_concern_drugs]), length(balanced.drugs.list$drugs_testing[balanced.drugs.list$drugs_testing %in% drug.dataset$no_concern_drugs]) )
compounds_df[nrow(compounds_df)+1,] <- c( "GUILDify", length(drug.dataset$drugs), length(drug.dataset$most_concern_drugs), length(drug.dataset$less_concern_drugs), length(drug.dataset$no_concern_drugs), length(drug.dataset$independent_drugs),
                                          length(balanced.drugs.list$drugs_training_list[[1]]), length(balanced.drugs.list$drugs_training_list[[1]][balanced.drugs.list$drugs_training_list[[1]] %in% drug.dataset$most_concern_drugs]), length(balanced.drugs.list$drugs_training_list[[1]][balanced.drugs.list$drugs_training_list[[1]] %in% drug.dataset$less_concern_drugs]), length(balanced.drugs.list$drugs_training_list[[1]][balanced.drugs.list$drugs_training_list[[1]] %in% drug.dataset$no_concern_drugs]),
                                          length(balanced.drugs.list$drugs_testing), length(balanced.drugs.list$drugs_testing[balanced.drugs.list$drugs_testing %in% drug.dataset$most_concern_drugs]), length(balanced.drugs.list$drugs_testing[balanced.drugs.list$drugs_testing %in% drug.dataset$less_concern_drugs]), length(balanced.drugs.list$drugs_testing[balanced.drugs.list$drugs_testing %in% drug.dataset$no_concern_drugs]) )
compounds_df[nrow(compounds_df)+1,] <- c( "SMILES", length(tanimoto_df$pert_iname), length(tanimoto_df[tanimoto_df$dilirank == "Most-DILI-Concern",]$pert_iname), length(tanimoto_df[tanimoto_df$dilirank == "Less-DILI-Concern",]$pert_iname), length(tanimoto_df[tanimoto_df$dilirank == "No-DILI-Concern",]$pert_iname), length(tanimoto_ind_df$pert_iname),
                                          length(datasets.list.smiles$training_datasets_with_names[[1]]$pert_iname), length(datasets.list.smiles$training_datasets_with_names[[1]][datasets.list.smiles$training_datasets_with_names[[1]]$dilirank == "Most-DILI-Concern",]$dilirank), length(datasets.list.smiles$training_datasets_with_names[[1]][datasets.list.smiles$training_datasets_with_names[[1]]$dilirank == "Less-DILI-Concern",]$dilirank), length(datasets.list.smiles$training_datasets_with_names[[1]][datasets.list.smiles$training_datasets_with_names[[1]]$dilirank == "No-DILI-Concern",]$dilirank), 
                                          length(datasets.list.smiles$testing_with_names$pert_iname), length(datasets.list.smiles$testing_with_names$dilirank[datasets.list.smiles$testing_with_names$dilirank == "Most-DILI-Concern"]), length(datasets.list.smiles$testing_with_names$dilirank[datasets.list.smiles$testing_with_names$dilirank == "Less-DILI-Concern"]), length(datasets.list.smiles$testing_with_names$dilirank[datasets.list.smiles$testing_with_names$dilirank == "No-DILI-Concern"]) )
compounds_df[nrow(compounds_df)+1,] <- c( "Targets", length(targets_df$pert_iname), length(targets_df[targets_df$dilirank == "Most-DILI-Concern",]$pert_iname), length(targets_df[targets_df$dilirank == "Less-DILI-Concern",]$pert_iname), length(targets_df[targets_df$dilirank == "No-DILI-Concern",]$pert_iname), length(targets_ind_df$pert_iname),
                                          length(datasets.list.targets$training_datasets_with_names[[1]]$pert_iname), length(datasets.list.targets$training_datasets_with_names[[1]][datasets.list.targets$training_datasets_with_names[[1]]$dilirank == "Most-DILI-Concern",]$dilirank), length(datasets.list.targets$training_datasets_with_names[[1]][datasets.list.targets$training_datasets_with_names[[1]]$dilirank == "Less-DILI-Concern",]$dilirank), length(datasets.list.targets$training_datasets_with_names[[1]][datasets.list.targets$training_datasets_with_names[[1]]$dilirank == "No-DILI-Concern",]$dilirank), 
                                          length(datasets.list.targets$testing_with_names$pert_iname), length(datasets.list.targets$testing_with_names$dilirank[datasets.list.targets$testing_with_names$dilirank == "Most-DILI-Concern"]), length(datasets.list.targets$testing_with_names$dilirank[datasets.list.targets$testing_with_names$dilirank == "Less-DILI-Concern"]), length(datasets.list.targets$testing_with_names$dilirank[datasets.list.targets$testing_with_names$dilirank == "No-DILI-Concern"]) )
compounds_df


### Calculate the drugs without targets ###
drugs_no_target <- drug.dataset$drugs[!(drug.dataset$drugs %in% targets_df$pert_iname)]
independent_drugs_no_target <- drug.dataset$independent_drugs[!(drug.dataset$independent_drugs %in% targets_ind_df$pert_iname)]
drugs_no_target
independent_drugs_no_target
drug.dataset$most_concern_drugs[!(drug.dataset$most_concern_drugs %in% targets_df$pert_iname)]
drug.dataset$less_concern_drugs[!(drug.dataset$less_concern_drugs %in% targets_df$pert_iname)]
drug.dataset$no_concern_drugs[!(drug.dataset$no_concern_drugs %in% targets_df$pert_iname)]



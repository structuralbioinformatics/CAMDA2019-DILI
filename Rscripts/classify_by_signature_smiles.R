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
  main_directory = "/Users/quim/Dropbox/UPF/PHD/Projects/camda"
  bigmem_directory = "/Users/quim/Documents/DATA/camda"
}


### Define files ###
drugs_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-pert_iname.rda", sep="/")
drug_info_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-info.rda", sep="/")
dilirank_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-dilirank.v2.rda", sep="/")
expression_file <- paste(bigmem_directory, "CAMDA_l1000_1314compounds-GSE92742_Level5_gct.rda", sep="/")
landmark_genes_file <- paste(main_directory, "additional_data/GSE92742_Broad_LINCS_gene_info_delta_landmark.txt", sep="/")
gene_info_file <- paste(main_directory, "additional_data/GSE92742_Broad_LINCS_gene_info.txt", sep="/")
cell_info_file <- paste(main_directory, "additional_data/GSE92742_Broad_LINCS_cell_info.txt", sep="/")
functions_file <- paste(main_directory, "Rscripts/camda_functions.R", sep="/")
sig_info_file <- paste(main_directory, "camda_data/GSE92742_Broad_LINCS_sig_info.txt", sep="/")
sig_metrics_file <- paste(main_directory, "camda_data/GSE92742_Broad_LINCS_sig_metrics.txt", sep="/")
# Specific data files
#wilcox_file <- paste(main_directory, "results/gene_test_landmark_phh_10_24_mixlessmost.tsv", sep="/")
wilcox_file <- paste(main_directory, "results/reverse_engineering/reverse_signature_phh_notcorrected.txt", sep="/")
#wilcox_file <- paste(main_directory, "results/reverse_engineering/reverse_signature_phh_noout_notcorrected.txt", sep="/")
tanimoto_file <- paste(main_directory, "/additional_data/tanimoto_smiles.tsv", sep="/")
# Output files
signature_smiles_data_file <- paste(main_directory, "/additional_data/dili_landmark_smiles.tsv", sep="/")
output.cv.rf <- paste(main_directory, "results/crossvalidation/cv_signature_smiles_rf.txt", sep="/")
output.cv.gbm <- paste(main_directory, "results/crossvalidation/cv_signature_smiles_gbm.txt", sep="/")
output.cv.glm <- paste(main_directory, "results/crossvalidation/cv_signature_smiles_glm.txt", sep="/")
output.ind.rf = paste(main_directory, "camda_data/independent_validation/JAguirre_predictions_signature_smiles_rf.txt", sep="/")
output.ind.gbm = paste(main_directory, "camda_data/independent_validation/JAguirre_predictions_signature_smiles_gbm.txt", sep="/")
output.ind.glm = paste(main_directory, "camda_data/independent_validation/JAguirre_predictions_signature_smiles_glm.txt", sep="/")


### Load files ###
source(functions_file)
load(drugs_file)
load(drug_info_file)
load(dilirank_file)
load(expression_file) # Requires cmapR


#### Subset drugs ####
drug.dataset <- subset.drug.dataset(drank.sel, outliers=outliers, remove.outliers=remove.outliers)


### Prepare SMILES data ###
tanimoto_df <- read.csv(tanimoto_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tanimoto_df$DILIConcern[tanimoto_df$DILIConcern == "No-DILI-concern"] <- "No-DILI-Concern" # Correct the 4 entries with lowercase concern
tanimoto_df$pert_iname = rownames(tanimoto_df)
colnames(tanimoto_df)[match("DILIConcern", colnames(tanimoto_df))] <- "dilirank"
tanimoto_df$vDILIConcern <- NULL
tanimoto_ind_df <- tanimoto_df[colnames(tanimoto_df) %in% drug.dataset$independent_drugs, rownames(tanimoto_df) %in% drug.dataset$drugs]
tanimoto_df<-tanimoto_df[colnames(tanimoto_df) %in% drug.dataset$drugs, rownames(tanimoto_df) %in% drug.dataset$drugs]


### Load wilcoxon test analysis ###
wilcox_df = read.csv(wilcox_file, header = FALSE, sep = "\t")
#selected_genes <- unique(wilcox_df$gene_id[wilcox_df$p.value<0.05 & wilcox_df$landmark==TRUE])
selected_genes <- unique(wilcox_df[,1])


### Subset GCT object ###
# Subset the GCT object by cell ID PHH, 10 µM, 24 h
expression_wilcox_df <- subset.expression(gct, selected_genes, drug.dataset$drugs, drug.dataset$dilirank_df, cell_id="PHH", pert_idose="10 µM", pert_itime="24 h", merge_samples = TRUE)
# Subset gene expression for independent drugs as well
expression_wilcox_ind_df <- subset.expression(gct, selected_genes, drug.dataset$independent_drugs, drug.dataset$independent_df, cell_id="PHH", pert_idose="10 µM", pert_itime="24 h", merge_samples = TRUE)


### Combine signature and SMILES data ###
# Combine the main datasets
tanimoto_df$dilirank <- NULL
tanimoto_df$severity <- NULL
combined_data <- merge(x = expression_wilcox_df, y = tanimoto_df, by = "pert_iname")
# Combine the independent datasets
tanimoto_ind_df$dilirank <- NULL
tanimoto_ind_df$severity <- NULL
combined_data_ind <- merge(x = expression_wilcox_ind_df, y = tanimoto_ind_df, by = "pert_iname")
# Combine everything and save it as a separated data file
signature_smiles_df <- rbind(combined_data, combined_data_ind)
signature_smiles_df$severity <- NULL
signature_smiles_df$dilirank[signature_smiles_df$pert_iname %in% combined_data_ind$pert_iname] <- "Ambiguous DILI-concern"
wilcox_genes <- substring(names(expression_wilcox_df)[c(4:length(names(expression_wilcox_df)))], 2)
col_names <- c('DrugName', 'DILIrank', wilcox_genes, names(combined_data)[c((length(wilcox_genes)+4):length(combined_data))])
names(signature_smiles_df) <- col_names
write.table(signature_smiles_df, file = signature_smiles_data_file,row.names=FALSE, na="-",col.names=TRUE, sep="\t")


### Prepare balanced machine learning datasets ###
set.seed(21)
datasets.list <- prepare.balanced.datasets(combined_data, number.repetitions, drug.dataset$most_concern_drugs, drug.dataset$less_concern_drugs, drug.dataset$no_concern_drugs, type_analysis = "discrete", fraction_train=fraction_train)


### Train the machine learning classifiers ###
# Train models by RF
signature.rf.results <- train.combine.classifiers(datasets.list, ml.method="rf", type_analysis="discrete", number.cv=number.cv, number.repetitions=number.repetitions)
write.machine.learning.results(signature.rf.results, output.cv.rf)
# Test RF model on independent dataset
signature.ind.rf <- validate.classifier(combined_data_ind, 
                                        output.ind.rf, 
                                        signature.rf.results$modFit.comb, 
                                        signature.rf.results$train.list,
                                        type_analysis="discrete")
# Train models by GBM
signature.gbm.results <- train.combine.classifiers(datasets.list, ml.method="gbm", type_analysis="discrete", number.cv=number.cv, number.repetitions=number.repetitions)
write.machine.learning.results(signature.gbm.results, output.cv.gbm)
# Test GBM model on independent dataset
signature.ind.gbm <- validate.classifier(combined_data_ind, 
                                         output.ind.gbm, 
                                         signature.gbm.results$modFit.comb, 
                                         signature.gbm.results$train.list,
                                         type_analysis="discrete")
# Train models by GLM
signature.glm.results <- train.combine.classifiers(datasets.list, ml.method="glm", type_analysis="discrete", number.cv=number.cv, number.repetitions=number.repetitions)
write.machine.learning.results(signature.glm.results, output.cv.glm)
# Test GLM model on independent dataset
signature.ind.glm <- validate.classifier(combined_data_ind, 
                                         output.ind.glm, 
                                         signature.glm.results$modFit.comb, 
                                         signature.glm.results$train.list,
                                         type_analysis="discrete")


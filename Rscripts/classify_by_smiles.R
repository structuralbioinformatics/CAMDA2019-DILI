### Load packages ###
library(cmapR)
library(ggplot2)
library(caret)


### Define variables ###
place = "home" #home or work
remove.outliers = TRUE
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
sig_info_file <- paste(main_directory, "camda_data/GSE92742_Broad_LINCS_sig_info.txt", sep="/")
sig_metrics_file <- paste(main_directory, "camda_data/GSE92742_Broad_LINCS_sig_metrics.txt", sep="/")
# Specific data files
tanimoto_file <- paste(main_directory, "/additional_data/tanimoto_smiles.tsv", sep="/")
# Output files
smiles_data_file <- paste(main_directory, "/additional_data/smiles.tsv", sep="/")
output.cv.rf <- paste(main_directory, "results/crossvalidation/cv_smiles_rf.txt", sep="/")
output.cv.gbm <- paste(main_directory, "results/crossvalidation/cv_smiles_gbm.txt", sep="/")
output.cv.glm <- paste(main_directory, "results/crossvalidation/cv_smiles_glm.txt", sep="/")
output.ind.rf = paste(main_directory, "camda_data/independent_validation/JAguirre_predictions_smiles_rf.txt", sep="/")
output.ind.gbm = paste(main_directory, "camda_data/independent_validation/JAguirre_predictions_smiles_gbm.txt", sep="/")
output.ind.glm = paste(main_directory, "camda_data/independent_validation/JAguirre_predictions_smiles_glm.txt", sep="/")


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
# Save data as a separated file
smiles_data_df1 <- rbind(tanimoto_df, tanimoto_ind_df)
smiles_data_df1$dilirank[smiles_data_df1$pert_iname %in% tanimoto_ind_df$pert_iname] <- "Ambiguous DILI-concern"
smiles_data_df1$severity <- NULL
smiles_df <- drank.sel[c("pert_iname", "SMILES")][drank.sel$pert_iname %in% signature_smiles_df1$pert_iname,]
smiles_data_df2 <- merge(x = smiles_data_df1, y = smiles_df, by = "pert_iname") # Include the SMILES
smiles_data_df2 <- smiles_data_df2[match(smiles_data_df1$pert_iname, smiles_data_df2$pert_iname),] # To maintain the order of the drugs (first the training set, then the hold-out set)
smiles_data_final_df <- data.frame(DrugName = smiles_data_df2$pert_iname, DILIrank = smiles_data_df2$dilirank, SMILES = smiles_data_df2$SMILES, smiles_data_df2[!(names(smiles_data_df2) %in% c("pert_iname", "dilirank", "SMILES"))])
write.table(smiles_data_final_df, file = smiles_data_file,row.names=FALSE, na="-",col.names=TRUE, sep="\t")


### Prepare balanced machine learning datasets ###
datasets.list <- prepare.balanced.datasets(tanimoto_df, number.repetitions, drug.dataset$most_concern_drugs, drug.dataset$less_concern_drugs, drug.dataset$no_concern_drugs, type_analysis = "discrete", fraction_train=fraction_train)


### Train the machine learning classifiers ###
# Train models by RF
smiles.rf.results <- train.combine.classifiers(datasets.list, ml.method="rf", type_analysis="discrete", number.cv=number.cv, number.repetitions=number.repetitions)
write.machine.learning.results(smiles.rf.results, output.cv.rf)
# Test RF model on independent dataset
smiles.ind.rf <- validate.classifier(tanimoto_ind_df, 
                                        output.ind.rf, 
                                        smiles.rf.results$modFit.comb, 
                                        smiles.rf.results$train.list,
                                        type_analysis="discrete")
# Train models by GBM
smiles.gbm.results <- train.combine.classifiers(datasets.list, ml.method="gbm", type_analysis="discrete", number.cv=number.cv, number.repetitions=number.repetitions)
write.machine.learning.results(smiles.gbm.results, output.cv.gbm)
# Test GBM model on independent dataset
smiles.ind.gbm <- validate.classifier(tanimoto_ind_df,
                                      output.ind.gbm, 
                                      smiles.gbm.results$modFit.comb, 
                                      smiles.gbm.results$train.list,
                                      type_analysis="discrete")
# Train models by GLM
smiles.glm.results <- train.combine.classifiers(datasets.list, ml.method="glm", type_analysis="discrete", number.cv=number.cv, number.repetitions=number.repetitions)
write.machine.learning.results(smiles.glm.results, output.cv.glm)
# Test GLM model on independent dataset
smiles.ind.glm <- validate.classifier(tanimoto_ind_df, 
                                       output.ind.glm, 
                                       smiles.glm.results$modFit.comb, 
                                       smiles.glm.results$train.list,
                                       type_analysis="discrete")


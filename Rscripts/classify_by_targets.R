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
# Specific data files
targets_file <- paste(main_directory, "additional_data/targets/targets_dgidb_hitpick_sea.tsv", sep="/")
# Output files
targets_data_file <- paste(main_directory, "/additional_data/targets.tsv", sep="/")
output.cv.rf <- paste(main_directory, "results/crossvalidation/cv_targets_rf.txt", sep="/")
output.cv.gbm <- paste(main_directory, "results/crossvalidation/cv_targets_gbm.txt", sep="/")
output.cv.glm <- paste(main_directory, "results/crossvalidation/cv_targets_glm.txt", sep="/")
output.ind.rf = paste(main_directory, "camda_data/independent_validation/JAguirre_predictions_targets_rf.txt", sep="/")
output.ind.gbm = paste(main_directory, "camda_data/independent_validation/JAguirre_predictions_targets_gbm.txt", sep="/")
output.ind.glm = paste(main_directory, "camda_data/independent_validation/JAguirre_predictions_targets_glm.txt", sep="/")


### Load files ###
source(functions_file)
load(drugs_file)
load(drug_info_file)
load(dilirank_file)
load(expression_file) # Requires cmapR


#### Subset drugs ####
drug.dataset <- subset.drug.dataset(drank.sel, outliers=outliers, remove.outliers=remove.outliers)


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
# Save it as a separated data file
targets_data_df <- targets_df
targets_data_df$severity <- NULL
targets_ind_data_df <- targets_ind_df[ order(targets_ind_df$pert_iname), ]
targets_ind_data_df$dilirank <- "Ambiguous DILI-concern"
targets_combined_df <- rbind(targets_data_df, targets_ind_data_df)
targets_combined_final_df <- data.frame(DrugName = targets_combined_df$pert_iname, DILIrank = targets_combined_df$dilirank, targets_combined_df[!(names(targets_combined_df) %in% c("pert_iname", "dilirank"))])
write.table(targets_combined_final_df, file = targets_data_file,row.names=FALSE, na="-",col.names=TRUE, sep="\t")


### Prepare balanced machine learning datasets ###
datasets.list <- prepare.balanced.datasets(targets_df, number.repetitions, drug.dataset$most_concern_drugs, drug.dataset$less_concern_drugs, drug.dataset$no_concern_drugs, type_analysis = "discrete", fraction_train=fraction_train)


### Train the machine learning classifiers ###
# Train models by RF
targets.rf.results <- train.combine.classifiers(datasets.list, ml.method="rf", type_analysis="discrete", number.cv=number.cv, number.repetitions=number.repetitions)
write.machine.learning.results(targets.rf.results, output.cv.rf)
# Test RF model on independent dataset
targets.ind.rf <- validate.classifier(targets_ind_df, 
                                     output.ind.rf, 
                                     targets.rf.results$modFit.comb, 
                                     targets.rf.results$train.list,
                                     type_analysis="discrete")
# Train models by GBM
targets.gbm.results <- train.combine.classifiers(datasets.list, ml.method="gbm", type_analysis="discrete", number.cv=number.cv, number.repetitions=number.repetitions)
write.machine.learning.results(targets.gbm.results, output.cv.gbm)
# Test GBM model on independent dataset
targets.ind.gbm <- validate.classifier(targets_ind_df,
                                      output.ind.gbm, 
                                      targets.gbm.results$modFit.comb, 
                                      targets.gbm.results$train.list,
                                      type_analysis="discrete")
# Train models by GLM
targets.glm.results <- train.combine.classifiers(datasets.list, ml.method="glm", type_analysis="discrete", number.cv=number.cv, number.repetitions=number.repetitions)
write.machine.learning.results(targets.glm.results, output.cv.glm)
# Test GLM model on independent dataset
targets.ind.glm <- validate.classifier(targets_ind_df, 
                                         output.ind.glm, 
                                         targets.glm.results$modFit.comb, 
                                         targets.glm.results$train.list,
                                         type_analysis="discrete")


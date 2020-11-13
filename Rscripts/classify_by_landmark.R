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


#### Get landmark genes ####
gene_info_df <- read.csv(gene_info_file, header=TRUE, sep="\t")
landmark_genes <- gene_info_df$pr_gene_id[gene_info_df$pr_is_lm==1]


#### Subset drugs ####
drug.dataset <- subset.drug.dataset(drank.sel, outliers=outliers, remove.outliers=remove.outliers)


### Subset GCT object ###
# Subset the GCT object by cell ID PHH, 10 µM, 24 h
expression_df <- subset.expression(gct, landmark_genes, drug.dataset$drugs, drug.dataset$dilirank_df, cell_id="PHH", pert_idose="10 µM", pert_itime="24 h", merge_samples = TRUE)
# Subset gene expression for independent drugs as well
expression_ind_df <- subset.expression(gct, landmark_genes, drug.dataset$independent_drugs, drug.dataset$independent_df, cell_id="PHH", pert_idose="10 µM", pert_itime="24 h", merge_samples = TRUE)


### Prepare balanced machine learning datasets ###
datasets.list <- prepare.balanced.datasets(expression_df, number.repetitions, drug.dataset$most_concern_drugs, drug.dataset$less_concern_drugs, drug.dataset$no_concern_drugs, type_analysis = "discrete", fraction_train=fraction_train)


### Train the machine learning classifiers ###
# Train models by RF
landmark.rf.results <- train.combine.classifiers(datasets.list, ml.method="rf", type_analysis="discrete", number.cv=number.cv, number.repetitions=number.repetitions)
write.machine.learning.results(landmark.rf.results, output.cv.rf)
# Test RF model on independent dataset
landmark.ind.rf <- validate.classifier(expression_ind_df, 
                                       output.ind.rf, 
                                       landmark.rf.results$modFit.comb, 
                                       landmark.rf.results$train.list,
                                       type_analysis="discrete")
# Train models by GBM
landmark.gbm.results <- train.combine.classifiers(datasets.list, ml.method="gbm", type_analysis="discrete", number.cv=number.cv, number.repetitions=number.repetitions)
write.machine.learning.results(landmark.gbm.results, output.cv.gbm)
# Test GBM model on independent dataset
landmark.ind.gbm <- validate.classifier(expression_ind_df, 
                                        output.ind.gbm, 
                                        landmark.gbm.results$modFit.comb, 
                                        landmark.gbm.results$train.list,
                                        type_analysis="discrete")


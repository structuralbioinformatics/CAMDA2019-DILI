### Load packages ###
library(cmapR)
library(ggplot2)
library(caret)


### Define variables ###
remove.outliers = TRUE
outliers = c('daunorubicin', 'vorinostat')
number.cv = 10
number.repetitions = 10
fraction_train = 0.7
num_top <- 2
corr_threshold <- 0.5


### Define files ###
# Data files
drugs_file <- "/home/quim/PHD/Projects/camda/camda_data/CAMDA_l1000_1314compounds-pert_iname.rda"
drug_info_file <- "/sbi/users/interchange/emre/quim/camda/CAMDA_l1000_1314compounds-info.rda"
dilirank_file <- "/home/quim/PHD/Projects/camda/camda_data/CAMDA_l1000_1314compounds-dilirank.v2.rda"
expression_file <- "/sbi/users/interchange/emre/quim/camda/CAMDA_l1000_1314compounds-GSE92742_Level5_gct.rda"
gene_info_file <- "/home/quim/PHD/Projects/camda/additional_data/GSE92742_Broad_LINCS_gene_info.txt"
cell_info_file <- "/home/quim/PHD/Projects/camda/additional_data/GSE92742_Broad_LINCS_cell_info.txt"
functions_file <- "/home/quim/PHD/Projects/camda/Rscripts/camda_functions.R"
# Specific data files
reverse_file = sprintf("/home/quim/PHD/Projects/camda/results/reverse_engineering/gene_test_top%s_%s_correlated_samples.tsv", num_top, corr_threshold)
input_corr_samples = sprintf("/home/quim/PHD/Projects/camda/results/top%s_%s_correlated_samples.txt", num_top, corr_threshold)
input.list <- sprintf("/home/quim/PHD/Projects/camda/results/reverse_engineering/reverse_signature_top%s_%s_correlated_fdr.txt", num_top, corr_threshold)
tanimoto_file <- "/home/quim/PHD/Projects/camda/additional_data/tanimoto_smiles.tsv"
# Output files
output.cv.rf <- sprintf("/home/quim/PHD/Projects/camda/results/crossvalidation/cv_signature_smiles_top%s_%s_correlated_rf.txt", num_top, corr_threshold)
output.cv.gbm <- sprintf("/home/quim/PHD/Projects/camda/results/crossvalidation/cv_signature_smiles_top%s_%s_correlated_gbm.txt", num_top, corr_threshold)
output.ind.rf = sprintf("/home/quim/PHD/Projects/camda/camda_data/independent_validation/JAguirre_predictions_signature_smiles_top%s_%s_correlated_rf.txt", num_top, corr_threshold)
output.ind.gbm = sprintf("/home/quim/PHD/Projects/camda/camda_data/independent_validation/JAguirre_predictions_signature_smiles_top%s_%s_correlated_gbm.txt", num_top, corr_threshold)


### Load files ###
source(functions_file)
load(drugs_file)
load(drug_info_file)
load(dilirank_file)
load(expression_file) # Requires cmapR


### Load wilcoxon test analysis by reverse engineering ###
selected_genes_fdr = read.csv(input.list, header = FALSE, stringsAsFactors = FALSE)[,1]
selected_genes <- selected_genes_fdr


### Get correlated samples and subset drugs ###
# Load correlated samples
corr_samples = read.csv(input_corr_samples, header = FALSE, stringsAsFactors = FALSE)[,1]
# Obtain drug dataset and subset it by drugs from correlated samples
drug.dataset <- subset.drug.dataset(drank.sel, outliers=outliers, remove.outliers=remove.outliers)


### Prepare SMILES data ###
tanimoto_df <- read.csv(tanimoto_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tanimoto_df$DILIConcern[tanimoto_df$DILIConcern == "No-DILI-concern"] <- "No-DILI-Concern" # Correct the 4 entries with lowercase concern
tanimoto_df$pert_iname = rownames(tanimoto_df)
colnames(tanimoto_df)[match("DILIConcern", colnames(tanimoto_df))] <- "dilirank"
tanimoto_df$vDILIConcern <- NULL
tanimoto_ind_df <- tanimoto_df[colnames(tanimoto_df) %in% drug.dataset$independent_drugs, rownames(tanimoto_df) %in% drug.dataset$drugs]
tanimoto_df<-tanimoto_df[colnames(tanimoto_df) %in% drug.dataset$drugs, rownames(tanimoto_df) %in% drug.dataset$drugs]


### Subset GCT object ###
# Subset the GCT object by correlated samples
expression_reverse_df <- subset.expression.by.samples(gct, selected_genes, corr_samples, drug.dataset$drugs, drug.dataset$dilirank_df, merge_samples=TRUE, merge_criteria="median")
# Subset gene expression for independent drugs as well
expression_reverse_ind_df <- subset.expression.by.samples(gct, selected_genes, corr_samples, drug.dataset$independent_drugs, drug.dataset$independent_df, merge_samples = TRUE, merge_criteria="median")


### Combine signature and SMILES data ###
# Combine the main datasets
tanimoto_df$dilirank <- NULL
tanimoto_df$severity <- NULL
combined_data <- merge(x = expression_reverse_df, y = tanimoto_df, by = "pert_iname")
# Combine the independent datasets
tanimoto_ind_df$dilirank <- NULL
tanimoto_ind_df$severity <- NULL
combined_data_ind <- merge(x = expression_reverse_ind_df, y = tanimoto_ind_df, by = "pert_iname")


### Prepare balanced machine learning datasets ###
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


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
phenotype2gene_file <- paste(main_directory, "guildify_data/phenotype2gene.tsv", sep="/")
# Output files
output.cv.rf <- paste(main_directory, "results/crossvalidation/disgenet/cv_disgenet_rf", sep="/")
output.cv.gbm <- paste(main_directory, "results/crossvalidation/disgenet/cv_disgenet_gbm", sep="/")
output.ind.rf = paste(main_directory, "camda_data/independent_validation/disgenet/JAguirre_predictions_disgenet_rf", sep="/")
output.ind.gbm = paste(main_directory, "camda_data/independent_validation/disgenet/JAguirre_predictions_disgenet_gbm", sep="/")
output.cv.comb.rf <- paste(main_directory, "results/crossvalidation/cv_disgenet_comb_rf.txt", sep="/")
output.cv.comb.gbm <- paste(main_directory, "results/crossvalidation/cv_disgenet_comb_gbm.txt", sep="/")
output.single.ind.rf = paste(main_directory, "camda_data/independent_validation/JAguirre_predictions_disgenet_single_rf.txt", sep="/")
output.single.ind.gbm = paste(main_directory, "camda_data/independent_validation/JAguirre_predictions_disgenet_single_gbm.txt", sep="/")
output.comb.ind.rf = paste(main_directory, "camda_data/independent_validation/JAguirre_predictions_disgenet_comb_rf.txt", sep="/")
output.comb.ind.gbm = paste(main_directory, "camda_data/independent_validation/JAguirre_predictions_disgenet_comb_gbm.txt", sep="/")


### Load files ###
source(functions_file)
load(drugs_file)
load(drug_info_file)
load(dilirank_file)
load(expression_file) # Requires cmapR


#### Subset drugs ####
drug.dataset <- subset.drug.dataset(drank.sel, outliers=outliers, remove.outliers=remove.outliers)


### Get genes associated to phenotypes from DisGeNET ###
phenotype2gene <- read.csv(phenotype2gene_file, header=TRUE, sep="\t")
phenotypes <- unique(phenotype2gene$diseaseid)
curated_phenotypes <- unique(phenotype2gene$diseaseid[phenotype2gene$source == "CURATED"])
#curated_phenotypes <- c("C0023890", "C0239946")


### Define balanced drugs for machine learning datasets ###
balanced.drugs.list <- prepare.balanced.drugs(number.repetitions, drug.dataset$most_concern_drugs, drug.dataset$less_concern_drugs, drug.dataset$no_concern_drugs, fraction_train=fraction_train)


### Create models for each phenotype ###
disgenet.datasets.list <- list()
disgenet.rf.list <- list()
disgenet.gbm.list <- list()
independent.datasets.list <- list()
for (phenotype in curated_phenotypes){
  curated_gene_ids <-unique(phenotype2gene$geneid[phenotype2gene$diseaseid==phenotype & phenotype2gene$source=="CURATED"])
  name <- unique(phenotype2gene$name[phenotype2gene$diseaseid==phenotype])
  name_and_extension <- paste(phenotype, '.txt', sep="")
  if (length(curated_gene_ids)>=10) {
    ### Subset GCT object ###
    # Subset the GCT object by cell ID PHH, 10 µM, 24 h
    expression_df <- subset.expression(gct, curated_gene_ids, drug.dataset$drugs, drug.dataset$dilirank_df, cell_id="PHH", pert_idose="10 µM", pert_itime="24 h", merge_samples = TRUE)
    # Subset gene expression for independent drugs as well
    expression_ind_df <- subset.expression(gct, curated_gene_ids, drug.dataset$independent_drugs, drug.dataset$independent_df, cell_id="PHH", pert_idose="10 µM", pert_itime="24 h", merge_samples = TRUE)
    independent.datasets.list[[length(independent.datasets.list)+1]] <- expression_ind_df

    ### Prepare balanced machine learning datasets ###
    datasets.list <- prepare.datasets.from.balanced.drugs(expression_df, balanced.drugs.list$drugs_testing, balanced.drugs.list$drugs_training_list, type_analysis="discrete")
    disgenet.datasets.list[[length(disgenet.datasets.list)+1]] <- datasets.list
    
    ### Train the machine learning classifiers ###
    # Train models by RF
    disgenet.rf.results <- train.combine.classifiers(datasets.list, ml.method="rf", type_analysis="discrete", number.cv=number.cv, number.repetitions=number.repetitions)
    disgenet.rf.list[[length(disgenet.rf.list)+1]] <- disgenet.rf.results
    output.cv.rf.phenotype <- paste(output.cv.rf, name_and_extension, sep = "_")
    write.machine.learning.results(disgenet.rf.results, output.cv.rf.phenotype)
    # Test RF model on independent dataset
    output.ind.rf.phenotype <- paste(output.ind.rf, name_and_extension, sep = "_")
    disgenet.ind.rf <- validate.classifier(expression_ind_df, 
                                           output.ind.rf.phenotype, 
                                           disgenet.rf.results$modFit.comb, 
                                           disgenet.rf.results$train.list,
                                           type_analysis="discrete")
    # Train models by GBM
    disgenet.gbm.results <- train.combine.classifiers(datasets.list, ml.method="gbm", type_analysis="discrete", number.cv=number.cv, number.repetitions=number.repetitions)
    disgenet.gbm.list[[length(disgenet.gbm.list)+1]] <- disgenet.gbm.results
    output.cv.gbm.phenotype <- paste(output.cv.gbm, name_and_extension, sep = "_")
    write.machine.learning.results(disgenet.gbm.results, output.cv.gbm.phenotype)
    # Test GBM model on independent dataset
    output.ind.gbm.phenotype <- paste(output.ind.gbm, name_and_extension, sep = "_")
    disgenet.ind.gbm <- validate.classifier(expression_ind_df, 
                                            output.ind.gbm.phenotype, 
                                            disgenet.gbm.results$modFit.comb, 
                                            disgenet.gbm.results$train.list,
                                            type_analysis="discrete")
  }
}


### Combine all models ###
# Combine models of rf
combination.results.rf <- combine.disgenet.models(disgenet.datasets.list, disgenet.rf.list, output.cv.comb.rf)
# Test rf models on independent dataset
disgenet.ind.rf <- validate.disgenet.combined.models(independent.datasets.list,
                                                     output.single.ind.rf,
                                                     output.comb.ind.rf,
                                                     combination.results.rf$modFit.single,
                                                     combination.results.rf$modFit.combined,
                                                     disgenet.rf.list)
# Combine models of gbm
combination.results.gbm <- combine.disgenet.models(disgenet.datasets.list, disgenet.gbm.list, output.cv.comb.gbm)
# Test gbm models on independent dataset
disgenet.ind.gbm <- validate.disgenet.combined.models(independent.datasets.list, 
                                                      output.single.ind.gbm,
                                                      output.comb.ind.gbm,
                                                      combination.results.gbm$modFit.single,
                                                      combination.results.gbm$modFit.combined,
                                                      disgenet.gbm.list)



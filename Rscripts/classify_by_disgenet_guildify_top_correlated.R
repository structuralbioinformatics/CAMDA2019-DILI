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
database = 'guildify' # disgenet or guildify
type_genes = 'curated' # curated or all
num_top <- 2
corr_threshold <- 0.5
minimum_genes <- 10



### Define files ###
# Data files
drugs_file <- "/home/quim/PHD/Projects/camda/camda_data/CAMDA_l1000_1314compounds-pert_iname.rda"
drug_info_file <- "/sbi/users/interchange/emre/quim/camda/CAMDA_l1000_1314compounds-info.rda"
dilirank_file <- "/home/quim/PHD/Projects/camda/camda_data/CAMDA_l1000_1314compounds-dilirank.v2.rda"
expression_file <- "/sbi/users/interchange/emre/quim/camda/CAMDA_l1000_1314compounds-GSE92742_Level5_gct.rda"
landmark_genes_file <- "/home/quim/PHD/Projects/camda/additional_data/GSE92742_Broad_LINCS_gene_info_delta_landmark.txt"
gene_info_file <- "/home/quim/PHD/Projects/camda/additional_data/GSE92742_Broad_LINCS_gene_info.txt"
cell_info_file <- "/home/quim/PHD/Projects/camda/additional_data/GSE92742_Broad_LINCS_cell_info.txt"
functions_file <- "/home/quim/PHD/Projects/camda/Rscripts/camda_functions.R"
# Specific data files
phenotype2gene_file <- "/home/quim/PHD/Projects/camda/guildify_data/phenotype2gene.tsv"
disease2gene_file <- "/home/quim/PHD/Projects/camda/guildify_data/disease2gene_disgenet_guildify.tsv"
redundantphenotypes_file <- "/home/quim/PHD/Projects/camda/guildify_data/redundant_phenotypes.tsv"
# Output files
output.cv.rf <- sprintf("/home/quim/PHD/Projects/camda/results/crossvalidation/correlated_samples/cv_%s_%s_top%s_%s_correlated_rf", database, type_genes, num_top, corr_threshold)
output.cv.gbm <- sprintf("/home/quim/PHD/Projects/camda/results/crossvalidation/correlated_samples/cv_%s_%s_top%s_%s_correlated_gbm", database, type_genes, num_top, corr_threshold)
output.ind.rf = sprintf("/home/quim/PHD/Projects/camda/camda_data/independent_validation/correlated_samples/JAguirre_predictions_%s_%s_top%s_%s_correlated_rf", database, type_genes, num_top, corr_threshold)
output.ind.gbm = sprintf("/home/quim/PHD/Projects/camda/camda_data/independent_validation/correlated_samples/JAguirre_predictions_%s_%s_top%s_%s_correlated_gbm", database, type_genes, num_top, corr_threshold)
output.cv.comb.rf <- sprintf("/home/quim/PHD/Projects/camda/results/crossvalidation/cv_%s_%s_comb_top%s_%s_correlated_rf", database, type_genes, num_top, corr_threshold)
output.cv.comb.gbm <- sprintf("/home/quim/PHD/Projects/camda/results/crossvalidation/cv_%s_%s_comb_top%s_%s_correlated_gbm", database, type_genes, num_top, corr_threshold)
output.single.ind.rf = sprintf("/home/quim/PHD/Projects/camda/camda_data/independent_validation/JAguirre_predictions_%s_%s_single_top%s_%s_correlated_rf", database, type_genes, num_top, corr_threshold)
output.single.ind.gbm = sprintf("/home/quim/PHD/Projects/camda/camda_data/independent_validation/JAguirre_predictions_%s_%s_single_top%s_%s_correlated_gbm", database, type_genes, num_top, corr_threshold)
output.comb.ind.rf = sprintf("/home/quim/PHD/Projects/camda/camda_data/independent_validation/JAguirre_predictions_%s_%s_comb_top%s_%s_correlated_rf", database, type_genes, num_top, corr_threshold)
output.comb.ind.gbm = sprintf("/home/quim/PHD/Projects/camda/camda_data/independent_validation/JAguirre_predictions_%s_%s_comb_top%s_%s_correlated_gbm", database, type_genes, num_top, corr_threshold)


### Load files ###
source(functions_file)
load(drugs_file)
load(drug_info_file)
load(dilirank_file)
load(expression_file) # Requires cmapR


### Get correlated samples ###
corr_samples = read.csv(input_corr_samples, header = FALSE, stringsAsFactors = FALSE)[,1]


#### Subset drugs ####
drug.dataset <- subset.drug.dataset(drank.sel, outliers=outliers, remove.outliers=remove.outliers)


### Get genes associated to phenotypes from DisGeNET ###
type_analysis <- paste(database, type_genes, sep='.')
disgenet2gene <- read.csv(disease2gene_file, header=TRUE, sep="\t", stringsAsFactors = F)
redundantphenotypes_df <- read.csv(redundantphenotypes_file, header=TRUE, sep="\t", stringsAsFactors = F)
disgenet2gene_unique <- disgenet2gene[!disgenet2gene$diseaseid %in% redundantphenotypes_df$redundant.phenotype,]
disgenet2gene_df <- disgenet2gene_unique[c("geneid", "diseaseid", "diseaseterm", type_analysis)]
colnames(disgenet2gene_df) <- c("geneid", "diseaseid", "diseaseterm", "type_analysis")
disgenet2gene_df <- disgenet2gene_df[disgenet2gene_df$type_analysis==1,]
phenotypes <- unique(disgenet2gene_df$diseaseid)


### Define balanced drugs for machine learning datasets ###
balanced.drugs.list <- prepare.balanced.drugs(number.repetitions, drug.dataset$most_concern_drugs, drug.dataset$less_concern_drugs, drug.dataset$no_concern_drugs, fraction_train=fraction_train)


### Create models for each phenotype ###
disgenet.datasets.list <- list()
disgenet.rf.list <- list()
disgenet.gbm.list <- list()
independent.datasets.list <- list()
for (phenotype in phenotypes){
  gene_ids <- unique(disgenet2gene_df$geneid[disgenet2gene_df$diseaseid == phenotype])
  name_and_extension <- paste(phenotype, '.txt', sep="")
  if (length(gene_ids)>=minimum_genes) {
    ### Subset GCT object ###
    # Subset the GCT object by correlated samples
    expression_df <- subset.expression.by.samples(gct, gene_ids, corr_samples, drug.dataset$drugs, drug.dataset$dilirank_df, merge_samples=TRUE, merge_criteria="median")
    # Subset gene expression for independent drugs as well
    expression_ind_df <- subset.expression.by.samples(gct, gene_ids, corr_samples, drug.dataset$independent_drugs, drug.dataset$independent_df, merge_samples = TRUE, merge_criteria="median")
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



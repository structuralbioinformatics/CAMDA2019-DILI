#### Load packages ####
# For cmap analysis
library(cmapR)
library(ggplot2)


### Define variables ###
place = "home" #home or work
remove.outliers = TRUE
outliers = c('daunorubicin', 'vorinostat')
num_top <- 2
corr_threshold <- 0.5

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
gene_info_file <- paste(main_directory, "additional_data/GSE92742_Broad_LINCS_gene_info.txt", sep="/")
cell_info_file <- paste(main_directory, "additional_data/GSE92742_Broad_LINCS_cell_info.txt", sep="/")
functions_file <- paste(main_directory, "Rscripts/camda_functions.R", sep="/")
# Output files
output_corr_samples <- paste(main_directory, sprintf("results/top%s_%s_correlated_samples.txt", num_top, corr_threshold), sep="/")
output_table <- paste(main_directory, sprintf("results/reverse_engineering/gene_test_top%s_%s_correlated_samples.tsv", num_top, corr_threshold), sep="/")
output.list.pval <- paste(main_directory, sprintf("results/reverse_engineering/reverse_signature_top%s_%s_correlated_notcorrected.txt", num_top, corr_threshold), sep="/")
output.list.fdr <- paste(main_directory, sprintf("results/reverse_engineering/reverse_signature_top%s_%s_correlated_fdr.txt", num_top, corr_threshold), sep="/")
output.list.bf <- paste(main_directory, sprintf("results/reverse_engineering/reverse_signature_top%s_%s_correlated_bf.txt", num_top, corr_threshold), sep="/")


### Load files ###
source(functions_file)
load(drugs_file)
load(drug_info_file)
load(dilirank_file)
load(expression_file) # Requires cmapR


### Get landmark genes ###
gene_info_df <- read.csv(gene_info_file, header=TRUE, sep="\t")
landmark_genes <- gene_info_df$pr_gene_id[gene_info_df$pr_is_lm==1]


### Subset drugs ###
drug.dataset <- subset.drug.dataset(drank.sel, outliers=outliers, remove.outliers=remove.outliers)


### Do not subset ###
#id_subset <- which(gct@cdesc$pert_idose == "10 µM" & gct@cdesc$pert_itime == "24 h")
#gct_subset <- subset.gct(gct, cid=id_subset)
gct_subset <- gct


### Function to get the n largest values and positions
### of a matrix.
### https://stackoverflow.com/questions/32544566/find-the-largest-values-on-a-matrix-in-r
nlargest <- function(m, n, sim = TRUE) {
  mult <- 1;
  if (sim) mult <- 2;
  res <- order(m, decreasing=TRUE)[seq_len(n) * mult];
  pos <- arrayInd(res, dim(m), useNames = TRUE);
  list(values = m[res],
       position = pos)
}

### Function obtain the top correlated samples of a list of samples.
check_top_samples <- function(nlarg_result, top_samples_to_check, top_samples, num_top){
  top_results <- data.frame(matrix(nrow=0, ncol=2))
  colnames(top_results) <- c("top_sample", "value")
  for (sample_to_check in top_samples_to_check){
    topres <- nlarg_result$position[sample_to_check == nlarg_result$position[,c("row")] | sample_to_check == nlarg_result$position[,c("col")],]
    pair_topres <- as.vector(topres[!(topres[,c("row")] %in% top_samples & topres[,c("col")] %in% top_samples),][1,])
    pos_topres <- which(pair_topres[1] == nlarg_result$position[,c("row")] & pair_topres[2] == nlarg_result$position[,c("col")])
    value_pair <- nlarg_result$values[pos_topres]
    sample_topres <- pair_topres[pair_topres != sample_to_check]
    top_results[nrow(top_results)+1,] <- c(sample_topres, value_pair)
  }
  top_results <- top_results[order(-top_results$value),]
  new_samples <- c()
  values_new_samples <- c()
  top_resulting_samples <- top_results$top_sample
  top_resulting_values <- top_results$value
  for (i in 1:nrow(top_results)){
    top_result <- top_resulting_samples[i]
    top_value <- top_resulting_values[i]
    if (length(top_samples) < num_top){
      if (!top_result %in% top_samples){
        top_samples <- append(top_samples, top_result)
        new_samples <- append(new_samples, top_result)
        values_new_samples <- append(values_new_samples, top_value)
      }
    }
  }
  return(list(top_samples=top_samples, new_samples=new_samples, values_new_samples=values_new_samples));
}

### Function to obtain all the top correlated samples of a drug.
### Needs the function "check_top_samples()"
get_top_correlated_samples <- function(nlarg_result, num_top){
  pair <- nlarg_result$position[1,] # Get the top pair
  value <- nlarg_result$values[1] # Get the top corr value
  sample1 <- as.numeric(pair[1]) # Get first top sample
  sample2 <- as.numeric(pair[2]) # Get second top sample
  top_samples <- c(sample1, sample2)
  top_samples_to_check <- c(sample1, sample2)
  repeat{
    checked_top_samples <- check_top_samples(nlarg_result, top_samples_to_check, top_samples, num_top)
    top_samples <- checked_top_samples$top_samples
    if (length(top_samples) >= num_top){
      break
    }
    print(checked_top_samples$values_new_samples)
    top_samples_to_check <- checked_top_samples$new_samples
  }
  return(top_samples);
}

### Function to obtain from a list of samples, either the top N samples correlated with them or below a threshold of correlation.
check_samples_by_threshold_and_top <- function(nlarg_result, top_samples_to_check, top_samples, num_top, corr_threshold=1){
  top_results <- data.frame(matrix(nrow=0, ncol=2))
  colnames(top_results) <- c("top_sample", "value")
  for (sample_to_check in top_samples_to_check){
    topres <- nlarg_result$position[sample_to_check == nlarg_result$position[,c("row")] | sample_to_check == nlarg_result$position[,c("col")],] #Get the top pairs containing the sample that we are checking
    rownames(topres) <- 1:nrow(topres) #Give them a row index
    pos_topres <- as.numeric(rownames(topres[!(topres[,c("row")] %in% top_samples & topres[,c("col")] %in% top_samples),])) #Get the position of only the ones that do not include any of the top_samples
    pairs_topres <- topres[!(topres[,c("row")] %in% top_samples & topres[,c("col")] %in% top_samples),] #Get the pairs that do not include any of the top_samples
    value_pairs <- nlarg_result$values[pos_topres] # Get the correlation value of the ones selected before
    samples_topres <- as.vector(unlist(apply(pairs_topres, 1, function(x){return(x[x!=sample_to_check])}))) #Get the sample number
    top_results_sample <- data.frame(top_sample=samples_topres, value=value_pairs) # Include the sample number and the correlation value in a table in order
    top_results <- rbind(top_results, top_results_sample) # Include the results in the table top_results
  }
  top_results <- top_results[order(-top_results$value),]
  new_samples <- c()
  values_new_samples <- c()
  top_resulting_samples <- top_results$top_sample
  top_resulting_values <- top_results$value
  # Get top N correlated samples OR samples below a corr_threshold!
  for (i in 1:nrow(top_results)){
    top_result <- top_resulting_samples[i]
    top_value <- top_resulting_values[i]
    if ((length(top_samples) < num_top) | (top_value > corr_threshold)){
      if (!top_result %in% top_samples){
        top_samples <- append(top_samples, top_result)
        new_samples <- append(new_samples, top_result)
        values_new_samples <- append(values_new_samples, top_value)
      }
    }
  }
  return(list(top_samples=top_samples, new_samples=new_samples, values_new_samples=values_new_samples));
}

### Function to obtain all the correlated samples of a drug by n_top and threshold.
### Needs the function "check_samples_by_threshold_and_top()"
get_correlated_samples_by_threshold_and_top <- function(nlarg_result, num_top, corr_threshold=1){
  pair <- nlarg_result$position[1,] # Get the top pair
  value <- nlarg_result$values[1] # Get the top corr value
  sample1 <- as.numeric(pair[1]) # Get first top sample
  sample2 <- as.numeric(pair[2]) # Get second top sample
  top_samples <- c(sample1, sample2)
  top_samples_to_check <- c(sample1, sample2)
  repeat{
    checked_top_samples <- check_samples_by_threshold_and_top(nlarg_result, top_samples_to_check, top_samples, num_top, corr_threshold)
    top_samples <- checked_top_samples$top_samples
    #print(checked_top_samples$values_new_samples)
    if (length(top_samples) >= num_top){
      break
    }
    top_samples_to_check <- checked_top_samples$new_samples
  }
  return(top_samples);
}


### Check n top correlated samples for each drug ###
all_drugs <- union(drug.dataset$drugs, drug.dataset$independent_drugs)
corr_samples <- c()
for (drug in all_drugs){
  # Subset by drug
  id_drug <- which(gct_subset@cdesc$pert_iname == drug)
  gct_drug <- subset.gct(gct_subset, cid=id_drug)
  # Calculate correlation between samples of the same drug
  drug_cor <- cor(gct_drug@mat)
  diag(drug_cor) <- 0 # Replace the values of the diagonal to 0
  n_samples <- nrow(drug_cor)
  n_pairs <- round(n_samples**2/2) # Number of pairs of samples
  nlarg_result <- nlargest(drug_cor, n_pairs) # Get the "n_pairs" of largest values and positions of a matrix ("drug_cor")
  top_samples <- get_correlated_samples_by_threshold_and_top(nlarg_result, num_top, corr_threshold=corr_threshold)
  top_sample_names <- rownames(drug_cor)[top_samples]
  corr_samples <- append(corr_samples, top_sample_names)
}
corr_samples <- unique(corr_samples)
# Write the resulting correlated samples
write.table(corr_samples, file = output_corr_samples,row.names=FALSE, na="-",col.names=FALSE, sep="\t")


### Subset by correlated samples ###
id_corr <- which(gct_subset@cid %in% corr_samples)
gct_corr <- subset.gct(gct_subset, cid=id_corr)


### Subset by DILI/NO-DILI drugs ###
# Subset DILI drugs
dili_concern_drugs <- union(drug.dataset$most_concern_drugs, drug.dataset$less_concern_drugs)
id_corr_dili <- which(gct_corr@cdesc$pert_iname %in% dili_concern_drugs)
gct_corr_dili <- subset.gct(gct_corr, cid=id_corr_dili)
# Subset NO-DILI drugs
id_corr_nodili <- which(gct_corr@cdesc$pert_iname %in% drug.dataset$no_concern_drugs)
gct_corr_nodili <- subset.gct(gct_corr, cid=id_corr_nodili)


### Find a gene signature by reverse engineering ###
# Define the new table
cols <- c("gene_id", "statistic", "p.value")
reverse_df <- data.frame(matrix(ncol = length(cols), nrow=0))
colnames(reverse_df) <- cols
# Calculate the wilcoxon for each landmark gene and store the results
for (gene_id in landmark_genes) {
  dili_expression <- gct_corr_dili@mat[gct_corr_dili@rid == gene_id]
  nodili_expression <- gct_corr_nodili@mat[gct_corr_nodili@rid == gene_id]
  res <- wilcox.test(dili_expression, nodili_expression, alternative = "two.sided", paired = FALSE)
  reverse_df[nrow(reverse_df)+1,] <- c(gene_id, res$statistic, res$p.value)
}
reverse_df$p.value <- as.numeric(reverse_df$p.value)
reverse_df$p.adjust.bf <- p.adjust(reverse_df$p.value, "bonferroni")
reverse_df$p.adjust.fdr <- p.adjust(reverse_df$p.value, "fdr")
write.table(reverse_df, file = output_table,row.names=FALSE, na="-",col.names=TRUE, sep="\t")


### Write lists of genes ###
# Not corrected
selected_genes_pval <- unique(reverse_df$gene_id[reverse_df$p.value<0.05])
write.table(selected_genes_pval, file = output.list.pval,row.names=FALSE, na="-",col.names=FALSE, sep="\t")
# Corrected by FDR
selected_genes_fdr <- unique(reverse_df$gene_id[reverse_df$p.adjust.fdr<0.05])
write.table(selected_genes_fdr, file = output.list.fdr,row.names=FALSE, na="-",col.names=FALSE, sep="\t")
# Corrected by Bonferroni
selected_genes_bf <- unique(reverse_df$gene_id[reverse_df$p.adjust.bf<0.05])
write.table(selected_genes_bf, file = output.list.bf,row.names=FALSE, na="-",col.names=FALSE, sep="\t")


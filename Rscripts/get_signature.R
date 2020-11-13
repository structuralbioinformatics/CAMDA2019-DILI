#### Load packages ####
# For cmap analysis
library(cmapR)
library(ggplot2)


### Define variables ###
place = "home" #home or work
remove.outliers = FALSE
outliers = c('daunorubicin', 'vorinostat')
number.cv = 10
number.repetitions = 10
fraction_train = 0.7
cell_id = "PHH"
pert_idose="10 µM"
pert_itime="24 h"


#### Define files ####
# Data files
if (place=="work"){
  main_directory = "/home/quim/PHD/Projects/camda"
  bigmem_directory = "/sbi/users/interchange/emre/quim/camda"
} else {
  main_directory = "/Users/quim/Dropbox/UPF/PHD/Projects/camda"
  bigmem_directory = "/Users/quim/Documents/DATA/camda"
}
drugs_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-pert_iname.rda", sep="/")
dilirank_file <- paste(main_directory, "camda_data/CAMDA_l1000_1314compounds-dilirank.v2.rda", sep="/")
expression_file <- paste(bigmem_directory, "CAMDA_l1000_1314compounds-GSE92742_Level5_gct.rda", sep="/")
gene_info_file <- paste(main_directory, "additional_data/GSE92742_Broad_LINCS_gene_info.txt", sep="/")
cell_info_file <- paste(main_directory, "additional_data/GSE92742_Broad_LINCS_cell_info.txt", sep="/")
functions_file <- paste(main_directory, "Rscripts/camda_functions.R", sep="/")
# Output files
if (remove.outliers==TRUE){
  output_table <- paste(main_directory, "results/reverse_engineering/gene_test_phh_noout.tsv", sep="/")
  output.list.pval <- paste(main_directory, "results/reverse_engineering/reverse_signature_phh_noout_notcorrected.txt", sep="/")
  output.list.pval.info <- paste(main_directory, "results/reverse_engineering/reverse_signature_phh_noout_notcorrected_info.txt", sep="/")
  output.list.fdr <- paste(main_directory, "results/reverse_engineering/reverse_signature_phh_noout_fdr.txt", sep="/")
  output.list.bf <- paste(main_directory, "results/reverse_engineering/reverse_signature_phh_noout_bf.txt", sep="/")
} else {
  output_table <- paste(main_directory, "results/reverse_engineering/gene_test_phh.tsv", sep="/")
  output.list.pval <- paste(main_directory, "results/reverse_engineering/reverse_signature_phh_notcorrected.txt", sep="/")
  output.list.pval.info <- paste(main_directory, "results/reverse_engineering/reverse_signature_phh_notcorrected_info.txt", sep="/")
  output.list.fdr <- paste(main_directory, "results/reverse_engineering/reverse_signature_phh_fdr.txt", sep="/")
  output.list.bf <- paste(main_directory, "results/reverse_engineering/reverse_signature_phh_bf.txt", sep="/")
}


### Load files ###
source(functions_file)
load(drugs_file)
load(dilirank_file)
load(expression_file) # Requires cmapR


### Get landmark genes ###
gene_info_df <- read.csv(gene_info_file, header=TRUE, sep="\t")
landmark_genes <- gene_info_df$pr_gene_id[gene_info_df$pr_is_lm==1]


### Subset drugs ###
drug.dataset <- subset.drug.dataset(drank.sel, outliers=outliers, remove.outliers=remove.outliers)


### Subset by landmark genes, cell ID, dose and time ###
id_selected_samples <- which(gct@cdesc$cell_id == cell_id & gct@cdesc$pert_idose == pert_idose & gct@cdesc$pert_itime == pert_itime)
id_selected_genes <- which(gct@rid %in% landmark_genes)
gct_subset <- subset.gct(gct, cid=id_selected_samples, rid=id_selected_genes)


### Subset by DILI/NO-DILI drugs ###
# Subset DILI drugs
dili_concern_drugs <- union(drug.dataset$most_concern_drugs, drug.dataset$less_concern_drugs)
id_subset_dili <- which(gct_subset@cdesc$pert_iname %in% dili_concern_drugs)
gct_subset_dili <- subset.gct(gct_subset, cid=id_subset_dili)
# Subset NO-DILI drugs
id_subset_nodili <- which(gct_subset@cdesc$pert_iname %in% drug.dataset$no_concern_drugs)
gct_subset_nodili <- subset.gct(gct_subset, cid=id_subset_nodili)


### Find a gene signature by reverse engineering ###
# Define the new table
cols <- c("gene_id", "statistic", "p.value")
reverse_df <- data.frame(matrix(ncol = length(cols), nrow=0))
colnames(reverse_df) <- cols
# Calculate the wilcoxon for each landmark gene and store the results
for (gene_id in landmark_genes) {
  dili_expression <- gct_subset_dili@mat[gct_subset_dili@rid == gene_id]
  nodili_expression <- gct_subset_nodili@mat[gct_subset_nodili@rid == gene_id]
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
selected_genes_pval_df <- gene_info_df[gene_info_df$pr_gene_id %in% selected_genes_pval,]
selected_genes_pval_df <- selected_genes_pval_df[c("pr_gene_id", "pr_gene_symbol", "pr_gene_title")]
write.table(selected_genes_pval_df, file = output.list.pval.info,row.names=FALSE, na="-",col.names=TRUE, sep="\t")
# Corrected by FDR
selected_genes_fdr <- unique(reverse_df$gene_id[reverse_df$p.adjust.fdr<0.05])
write.table(selected_genes_fdr, file = output.list.fdr,row.names=FALSE, na="-",col.names=FALSE, sep="\t")
# Corrected by Bonferroni
selected_genes_bf <- unique(reverse_df$gene_id[reverse_df$p.adjust.bf<0.05])
write.table(selected_genes_bf, file = output.list.bf,row.names=FALSE, na="-",col.names=FALSE, sep="\t")


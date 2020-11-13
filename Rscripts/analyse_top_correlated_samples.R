#### Load packages ####
# For cmap analysis
library(cmapR)
library(ggplot2)


### Define variables ###
remove.outliers = TRUE
outliers = c('daunorubicin', 'vorinostat')


#### Define files ####
# Data files
drugs_file <- "/home/quim/PHD/Projects/camda/camda_data/CAMDA_l1000_1314compounds-pert_iname.rda"
drug_info_file <- "/sbi/users/interchange/emre/quim/camda/CAMDA_l1000_1314compounds-info.rda"
dilirank_file <- "/home/quim/PHD/Projects/camda/camda_data/CAMDA_l1000_1314compounds-dilirank.v2.rda"
expression_file <- "/sbi/users/interchange/emre/quim/camda/CAMDA_l1000_1314compounds-GSE92742_Level5_gct.rda"
gene_info_file <- "/home/quim/PHD/Projects/camda/additional_data/GSE92742_Broad_LINCS_gene_info.txt"
cell_info_file <- "/home/quim/PHD/Projects/camda/additional_data/GSE92742_Broad_LINCS_cell_info.txt"
functions_file <- "/home/quim/PHD/Projects/camda/Rscripts/camda_functions.R"
# Output files
output_plot_pdf <- "/home/quim/PHD/Projects/camda/results/plots/analysis_top_correlated_samples.pdf"


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
nlargest <- function(m, n, sim = TRUE) {
  mult <- 1;
  if (sim) mult <- 2;
  res <- order(m, decreasing=TRUE)[seq_len(n) * mult];
  pos <- arrayInd(res, dim(m), useNames = TRUE);
  list(values = m[res],
       position = pos)
}

### Check n top correlated samples for each drug ###
all_drugs <- union(drug.dataset$drugs, drug.dataset$independent_drugs)
top_thresholds <- seq(2, 10, by=1)
table_top_corr <- data.frame(matrix(ncol=3, nrow=0))
colnames(table_top_corr) <- c("pert_iname", "correlation", "num_top")
for (drug in all_drugs){
  # Subset by drug
  id_drug <- which(gct_subset@cdesc$pert_iname == drug)
  gct_drug <- subset.gct(gct_subset, cid=id_drug)
  drug_cor <- cor(gct_drug@mat)
  diag(drug_cor) <- 0 # Replace the values of the diagonal to 0
  # Check the largest correlation values for each top threshold
  for (top_threshold in top_thresholds){
    nlarg_result <- nlargest(drug_cor, top_threshold) # Get the largest n values
    result_table <- data.frame(pert_iname=drug, correlation=nlarg_result$values, num_top=top_threshold)
    table_top_corr <- rbind(table_top_corr, result_table)
  }
}

### Get the number of drugs for each correlation threshold ###
corr_thresholds <- seq(0.1, 0.9, by=0.05)
table_num_drugs <- data.frame(matrix(ncol=3, nrow=0))
colnames(table_num_drugs) <- c("num_drugs", "correlation", "num_top")
for (corr_threshold in corr_thresholds){
  for (top_threshold in top_thresholds){
    num_drugs <- length(unique(table_top_corr[table_top_corr$num_top==top_threshold & table_top_corr$correlation>corr_threshold,]$pert_iname))
    result_table <- data.frame(num_drugs=num_drugs, correlation=corr_threshold, num_top=top_threshold)
    table_num_drugs <- rbind(table_num_drugs, result_table)
  }
}

table_num_drugs$num_top <- as.character(table_num_drugs$num_top)

### Graphic of number of different types of drugs depending on threshold ###
pdf(output_plot_pdf)

ggplot(table_num_drugs, aes(x=correlation, y=num_drugs, group=num_top)) +
  geom_line(aes(linetype=num_top))

dev.off()

